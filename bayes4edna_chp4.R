## - Bayesian models for eDNA (2023 samples)


## created 5 July 2024
## by Molly M Kressler

########
## Load data  
########
## cleaned in pipeline documents

pacman::p_load(sf,tidyverse,dplyr,ggplot2, patchwork, cowplot,lubridate,flextable, rnaturalearth, lme4, modelsummary, readr, readxl, brms, bayesplot, tidybayes, modelr, ggdist, nimble, MCMCvis)

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/')

## qPCR results, all species, not standards, test assay positive nor negative test controls 
fieldsamples <- read.csv('EDNA/data_edna/qPCRresults/processedQPCRresults_cornwallspecies_june2024.csv')%>%
	filter(NTC.Amplified == 0)%>%
	dplyr::select(-NTC.Amplified)%>%
		as_tibble() # st_write converted the TRUE/FALSE to 0s and 1s. 0 = TRUE and 1 = FALSE. This removes any plate where the NTC amplified before the LOQ. 


	## prepare data 
	data <- fieldsamples %>% 
		mutate(taxa = as.character(case_when(Target.Name == 'Engraulis' ~ 'Fish', Target.Name == 'Scomber' ~ 'Fish',Target.Name == 'Prionace' ~ 'Shark', Target.Name == 'Alopias' ~ 'Shark')))%>%
		mutate(log10.copies = log(copies.sampavg+0.1),
			round.copies = as.numeric(round(copies.sampavg,0)),
			re.level = as.numeric(as.numeric(factor(eventID))),
			methodID = as.numeric(as.numeric(factor(method.type))),
			taxaID = as.numeric(as.numeric(factor(taxa))))%>%
		dplyr::select('eventID', 're.level', 'Sample.Name', 'Target.Name','taxa', 'method.type', 'methodID', 'taxaID', 'copies','copies.techrepavg','log10.copies', 'round.copies')%>%
		group_by(Sample.Name, Target.Name)%>%
		slice(1) # nrow = 373
	data # waterbottles = 2, metaprobes = 1,fish = 1, shark = 2

	soaks <- fieldsamples %>%
		as_tibble%>%
		filter(method.type=='metaprobe')%>%
		dplyr::select(Sample.Name, Target.Name, eventID, copies, copies.techrepavg,method.type, time.timeIN, timeOUT)%>%
		mutate(timeOUT = as.character(paste0(timeOUT,':00')), time.timeIN = as.character(paste0(time.timeIN,':00')))%>%
		mutate(timeOUT2 = hms::as_hms(strptime(timeOUT, "%H:%M:%S")), timeIN2 = hms::as_hms(strptime(time.timeIN, "%H:%M:%S")), soaktime = as.duration(timeOUT2-timeIN2))%>%
		mutate(soaktime.min = as.numeric(soaktime)/60)%>%
		mutate(soaktime.hr = as.numeric(soaktime.min)/60)%>%
		mutate(re.level = as.numeric(as.numeric(factor(eventID))))



########
## Data exploration 
########

	hist<- ggplot(data = fieldsamples, aes(x =log(copies.sampavg+1), fill = Target.Name))+
		geom_histogram(position = 'dodge', bins = 5)+
	  	scale_fill_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus'))+
		theme_bw()+
		labs(y = NULL, x = 'DNA yield (log10)')

	ggsave(hist, file = 'EDNA/data_edna/brms/figures/hist_bySpecies.png', device = 'png', units = 'in', height = 6, width = 6, dpi = 800)


########
## BRMs models - normal distribution with log10 transformation 
########

	## export for dave Hudson
		#write.csv(data,'EDNA/data_edna/data_for_bayes4edna.csv')

	######
	### Fit Models
	######		

	## fit1 = gaussian with log10 data
	fit1 <- brm(log10.copies ~ method.type * taxa + (1|eventID), 
			data = data, 
			family = gaussian(), 
			file = 'EDNA/data_edna/brms/model1_gaussian_log10data.RDS') # auto-saves the model object as an RDS file

	fit1 <- readRDS('EDNA/data_edna/brms/model1_gaussian_log10data.RDS')

	## fit2 = poisson - hates this. wont sample. 
	fit2 <- brm(round.copies ~ method.type * taxa + (1|eventID), 
			data = data, 
			family = poisson(link = 'log'), 
			file = 'EDNA/data_edna/brms/model2_Poisson_notransformation.RDS') # auto-saves the model object as an RDS file

	# fit3, soaktime 
	soaks <- soaks %>% mutate(round.copies = round(copies.techrepavg, 0))
	fit3 <- brm(round.copies ~ soaktime.min + (1|eventID),
			data = soaks,
			family = poisson(link = 'log'),
			file = 'EDNA/data_edna/brms/model3_soaktime_poisson.RDS')


	######
	## Evaluate 
	######
		f3 <- readRDS('EDNA/data_edna/brms/model3_soaktime_poisson.RDS')

		summary_fit3 <- MCMCsummary(f3,round=3,pg0=TRUE,prob=c(0.05,0.95)) %>%
		      tibble::rownames_to_column()%>%
		      rename_with(str_to_title)%>%
		      rename(Parameter = Rowname)%>%
		      filter(!str_detect(Parameter,'-R'))%>%
		      filter(!str_detect(Parameter,'lp'))%>%
		      mutate(Parameter = case_when(
		     	Parameter == 'Intercept' ~ Parameter,
		     	str_detect(Parameter, 'soak') ~ 'Soak Time (min)',
		     	str_detect(Parameter, 'sd_eventID') ~ 'Random Effect \n\ of Event'
		      	))%>%
		      rename('pg0'='P>0')%>%
		      mutate(pg00 = case_when(Mean >= 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
		      rename('Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
		      mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
		      dplyr::select(-lower,-upper,-Sd, -pg0)%>%
		      flextable()%>%
		      theme_zebra()%>%
		      set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
		      align(align = 'center', part = 'all')%>%
		      font(fontname = 'Arial', part = 'all')%>%
		      color(color='black',part='all')%>%
		      fontsize(size = 10, part = 'all')%>%
		      autofit()
	    summary_fit3

		# Gather variable indices into a separate column (one row per draw, one column per parameter)
	
		draws <- fit3 %>%
			gather_draws(`b_Intercept`, `b_soaktime.min`)%>%
			median_qi(.width = c(.95, .5)) # calculates the median and the 95% quantile interval for the parameter

		caterpillars <- draws %>%
			dplyr::filter(.variable != 'r_eventID') %>%
			mutate(key = case_when(
					str_detect(.variable, 'Inte') ~ 'Intercept', 
					str_detect(.variable, 'soak') ~ 'Soak Time (min)'
				))%>%
			ggplot(aes(y = key, x = .value, xmin = .lower, xmax = .upper))+
			geom_pointinterval()+      
			geom_vline(xintercept=0,linetype=3)+
			theme_bw()+
			labs(y = NULL, x = 'Estimate & CI')
		caterpillars

		ggsave(caterpillars, file = 'EDNA/data_edna/brms/figures/caterpillars_brmsfit3_soaktimemin_july2024.png', device = 'png', units = 'in', height = 5, width = 5, dpi = 850)

		 save_as_image(summary_fit3, path = 'EDNA/data_edna/brms/figures/mcmcsummary_brmsfit3_soaktimemin_july2024.png', webshot = 'webshot2')
		 save_as_docx(summary_fit3, path = 'EDNA/data_edna/brms/figures/mcmcsummary_brmsfit3_soaktimemin_july2024.docx')

		## Diagnostics 

		## Overdispersion 
		n <- nrow(soaks)
		## fit3 
		res_fit3 <- resid(fit3)
		overdis_fit3 = sum(res_fit3^2)/(n-2)
		overdis_fit3 # 11081.38, overdispersed (Gaussian on log10 copies). needs a zlognormal distb I bet. leave for noe (23 July) and return in a day or two. 



########
## Nimble  models - poissons 
########
	# nimble model zlogNormal, method & taxa on copies, Dave Hudson. dependent on file: 'ZlogNormal_dist.R' written by Dave Hudson. 
	
		source("EDNA/eDNA_in_Cornish_blue_repo/ZlogNormal_dist.R")
		set.seed(11)
	
		modelcode_zlog <- nimbleCode({
		  # priors #
		  for(i in 1:E){
		    e[i] ~ dnorm(0, tau.re)
		    }
		  beta0 ~ dnorm(0, 1)
		  beta1 ~ dnorm(0, 1)
		  beta2 ~ dnorm(0, 1)
		  beta3 ~ dnorm(0, 1)
		  sdL ~ dexp(1)
		  pz ~ dunif(0, 1)
		  tau.re <- pow(sigma.re, -2)
		  sigma.re ~ dunif(0, 10)
		 
		  # likelihood #
		  for(i in 1:N){
		    Y[i] ~ dZlogNormal(pz, muL[i], sdL)
		    muL[i] <- beta0 + beta1 * methodID[i] + 
		                          beta2 * taxaID[i] + 
		                          beta3 * methodID[i] * taxaID[i] + 
		                          e[re[i]]
		 }
		 # derived parameters #
             b0 <- exp(beta0) # method 0, taxa 0
             b1 <- exp(beta0 + beta1) # method 1, taxa 0
             b2 <- exp(beta0 + beta2) # method 0, taxa 1
             b3 <- exp(beta0 + beta1 + beta2 + beta3) # method 1, taxa 1
		})
		
		fitZ.constants <- list(
			E = length(unique(data$re.level)), 
			N = nrow(data),
			methodID = data$methodID - 1, 
			taxaID = data$taxaID - 1,
			re = data$re.level)

		fitZ.data <- list(Y = data$copies.techrepavg)

		fitZ.init <- list(
			e = rnorm(15,0,1),
			beta0 = rnorm(1, 0, 1),
			beta1 = rnorm(1, 0, 1),
			beta2 = rnorm(1, 0, 1),
			beta3 = rnorm(1, 0, 1),
			b0 = rnorm(1, 0, 0.01),
			b1 = rgamma(1, 1,1),
			b2 = rgamma(1, 1,1),
			b3 = rgamma(1, 1,1),
			sdL = rexp(1, 1),
			sigma.re = runif(1, 0, 10),
			pz = runif(1, 0, 1)
			)

		### compile
		model_fitZ<-nimbleModel(code=modelcode_zlog, 
		                          name="model_fitZ",
		                          data=fitZ.data, 
		                          constants = fitZ.constants,
		                          inits=fitZ.init) #define the model
		  
		#model_fit2$initializeInfo()  

		CmfZ<-compileNimble(model_fitZ) 
		confZ <- configureMCMC(model_fitZ, monitors = c('sigma.re', 'beta1','beta2', 'beta3', 'beta0', 'pz', 'sdL'))#, onlySLICE = FALSE)
		MCMC_model_fitZ <- buildMCMC(confZ, na.rm = TRUE)
		ccMCMCZ <- compileNimble(MCMC_model_fitZ, project = model_fitZ)

		samplesZ <- runMCMC(ccMCMCZ,niter = 50000, nburnin = 12000,nchains = 3, samplesAsCodaMCMC = TRUE)

		MCMCsummary(samplesZ)		

		saveRDS(samplesZ, 'EDNA/data_edna/nimble_outputs_modelobjects/ZinflatedPoissonLogNormal_eDNA_50kiter_12kburn_3c_july2024.RDS')
	
	# summary and caterpillars

		samplesListZ <- readRDS('EDNA/data_edna/nimble_outputs_modelobjects/ZinflatedPoissonLogNormal_eDNA_50kiter_12kburn_3c_july2024.RDS')

		draws <- samplesListZ %>%
			gather_draws(`beta0`,`beta1`, `beta2`, `beta3`)%>%
			median_hdi(.width = c(.95, .5)) # calculates the median and the 95% quantile interval for the parameter

		spreads <- samplesListZ %>%
			spread_draws(`beta0`,`beta1`, `beta2`, `beta3`) # leaveunti Alice sends you the code to try to piece this together. 

		caterpillars <- draws %>%
			mutate(key = case_when(.variable == 'beta0' ~ 'Metaprobe, Fish', .variable == 'beta1' ~ 'Water Bottle, Fish', .variable == 'beta2' ~ 'Metaprobe, Shark', .variable == 'beta3' ~ 'Water Bottle, Shark'))%>%
			ggplot(aes(y = key, x = .value, xmin = .lower, xmax = .upper, col = key))+
			geom_pointinterval()+
		    scale_color_manual(values=c('#DAA507','#DAA507','#8EC7D2','#8EC7D2'))+
			geom_vline(xintercept=0,linetype=3)+
			theme_bw()+
			guides(col = 'none')+
			labs(x = 'Estimate', y = NULL)
		caterpillars

		ggsave(caterpillars, file = 'EDNA/data_edna/figures_and_tables/nimble_outputs/caterpillars_modelZ_july2024.png', device = 'png', units = 'in', height = 4, width = 5, dpi = 1080)
 
		summary_fitZ <- MCMCsummary(samplesListZ,round=3,pg0=TRUE,prob=c(0.05,0.95)) %>%
		      tibble::rownames_to_column()%>%
		      rename_with(str_to_title)%>%
		      rename(Parameter = Rowname)%>%
		      rename('% of posterior with \n\ same sign as estimate' = 'P>0', Estimate = 'Mean','lower'='5%',upper='95%')%>%
		      mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
		      dplyr::select(-lower,-upper,-Sd)%>% 
		      flextable()%>%
		      theme_zebra()%>%
		      set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
		      align(align = 'center', part = 'all')%>%
		      font(fontname = 'Arial', part = 'all')%>%
		      color(color='black',part='all')%>%
		      fontsize(size = 10, part = 'all')%>%
		      autofit()
		    summary_fitZ

		 save_as_image(summary_fitZ, path = 'EDNA/data_edna/figures_and_tables/nimble_outputs/mcmcsummary_modelZ_50kiter_12kburn_3c_july2024.png', webshot = 'webshot2')
		 save_as_docx(summary_fitZ, path = 'EDNA/data_edna/figures_and_tables/nimble_outputs/mcmcsummary_modelZ_50kiter_12kburn_3c_july2024.docx')

		maintext_fitZ <- MCMCsummary(samplesListZ,round=3,pg0=TRUE,prob=c(0.05,0.95))%>% 
		      tibble::rownames_to_column()%>%
		      rename_with(str_to_title)%>%
		      rename(Parameter = Rowname)%>%
			mutate(Case = case_when(Parameter == 'beta0' ~ 'Metaprobe, Fish', Parameter == 'beta1' ~ 'Water Bottle, Fish', Parameter == 'beta2' ~ 'Metaprobe, Shark', Parameter == 'beta3' ~ 'Water Bottle, Shark'), .before = 'Parameter')%>%
			filter(!is.na(Case))%>%
			dplyr::select(-Parameter)%>%
		      rename('pg0'='P>0')%>%
		      mutate(pg00 = case_when(Mean >= 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
		      rename('Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
		      mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
		      dplyr::select(-lower,-upper,-Sd, -pg0)%>%
		      flextable()%>%
		      theme_zebra()%>%
		      align(align = 'center', part = 'all')%>%
		      font(fontname = 'Arial', part = 'all')%>%
		      color(color='black',part='all')%>%
		      fontsize(size = 10, part = 'all')%>%
		      autofit()
		   maintext_fitZ
		 save_as_image(maintext_fitZ, path = 'EDNA/data_edna/figures_and_tables/nimble_outputs/plainlanguageParameter_mcmcsummary_modelZ_50kiter_12kburn_3c_july2024.png', webshot = 'webshot2')
		 save_as_docx(maintext_fitZ, path = 'EDNA/data_edna/figures_and_tables/nimble_outputs/plainlanguageParameter_mcmcsummary_modelZ_50kiter_12kburn_3c_july2024.docx')

		
	
##########################
##########################
		 
	##########
	# nimble model for poisson fit2 
		  modelcode_fit2 <- nimbleCode({
             # priors #
                 beta1 ~ dnorm(0, 0.01)
                 beta2 ~ dnorm(0, 0.01)
                 beta3 ~ dnorm(0, 0.01)
           
             # priors for random intercept, informative with Half-Cauchy priors for sigma
             for(i in 1:E){
                 e[i] ~ dnorm(0, tau.re)
             }
                 #beta0 ~ dnorm(0,.01) # acts as an intercept for tracing 
                 num ~ dnorm(0, 0.0016)
                 denom ~ dnorm(0,1)
                 sigma.re <- abs(num/denom)
                 tau.re <- 1/(sigma.re * sigma.re)
             
             # likelihood #
             for(i in 1:N){
                 Y[i] ~ dpois(mu[i])
 
                 log(mu[i]) <- e[re[i]] + beta1*methodID[i] + beta2*taxaID[i] +beta3*methodID[i]*taxaID[i] 
             }
             # derived parameters
             b0 <- exp(beta0) # method 0, taxa 0
             b1 <- exp(beta0 + beta1) # method 1, taxa 0
             b2 <- exp(beta0 + beta2) # method 0, taxa 1
             b3 <- exp(beta0 + beta1 + beta2 + beta3) # method 1, taxa 1
         })

		
		### no random effect
		
		modelcode_fit2_noRE <- nimbleCode({
             # priors #
                 beta1 ~ dnorm(0, 0.01)
                 beta2 ~ dnorm(0, 0.01)
                 beta3 ~ dnorm(0, 0.01)
           
             # priors for random intercept, informative with Half-Cauchy priors for sigma
                 beta0 ~ dnorm(0,.01) # acts as an intercept for tracing 
                 num ~ dnorm(0, 0.0016)
                 denom ~ dnorm(0,1)
                 sigma.re <- abs(num/denom)
                 tau.re <- 1/(sigma.re * sigma.re)
             
             # likelihood #
             for(i in 1:N){
                 Y[i] ~ dpois(mu[i])
 
                 log(mu[i]) <- beta0 + beta1*methodID[i] + beta2*taxaID[i] +beta3*methodID[i]*taxaID[i] 
             }
             # derived parameters
             b0 <- exp(beta0) # method 0, taxa 0
             b1 <- exp(beta0 + beta1) # method 1, taxa 0
             b2 <- exp(beta0 + beta2) # method 0, taxa 1
             b3 <- exp(beta0 + beta1 + beta2 + beta3) # method 1, taxa 1
         })


		### bits and bobs for compilation
		fit2.constants <- list(
				K = max(data$methodID),
				E = length(unique(data$re.level)), 
				N = nrow(data),
				methodID = data$methodID, 
				taxaID = data$taxaID,
				re = data$re.level)

		fit2.data <- list(Y = round(data$copies.techrepavg,0))

		fit2.init <- list(
			e = rnorm(15,0,1),
			beta0 = rnorm(1, 0, 0.01),
			beta1 = rnorm(1, 0, 0.01),
			beta2 = rnorm(1, 0, 0.01),
			beta3 = rnorm(1, 0, 0.01),
			b0 = rnorm(1, 0, 0.01),
			b1 = rgamma(1, 1,1),
			b2 = rgamma(1, 1,1),
			b3 = rgamma(1, 1,1),
			tau.re = rgamma(1,1,1),
			num = rnorm(1, 0, 0.0016),
			denom = rnorm(1, 0, 1)
			)

		### compile
	    model_fit2<-nimbleModel(code=modelcode_fit2, name="model_fit2",data=fit2.data, constants = fit2.constants,inits=fit2.init) #define the model

	    	model_fit2$calculate() # if = NA, indicates missing or invalid initial values, and you have to fix the model until it is numeric.
	    	
	    	Cmf2<-compileNimble(model_fit2) 
	    	conf2 <- configureMCMC(model_fit2, monitors = c('tau.re', 'beta0','beta1','beta2','beta3', 'b0', 'b1', 'b2', 'b3'), onlySLICE = FALSE)
		    MCMC_model_fit2 <- buildMCMC(conf2, na.rm = TRUE)
		    ccMCMC2 <- compileNimble(MCMC_model_fit2, project = model_fit2)
		    samples2b <- runMCMC(ccMCMC2, niter = 10000, nburnin = 5000, nchains = 3, samplesAsCodaMCMC = TRUE)

	    MCMCsummary(samples2b) 

	    saveRDS(samples2, 'EDNA/data_edna/nimble_outputs/model_fit2_2000iter_500burn_3ch_worksButDoesntConverge_july2024.RDS')

	    s2 <- readRDS('EDNA/data_edna/nimble_outputs/model_fit2_2000iter_500burn_3ch_worksButDoesntConverge_july2024.RDS')

	######
	### Exploration
	######	
	get_variables(fit1) # raw list of variable names. b_Intercept is the global mean, and the r_eventID[level, Intercept] are offsets of that mean for each condition. 

	## flextable of summary
	fit1table <- modelsummary(fit1, fmt=3,estimate='estimate', statistic='conf.int',conf_level=0.95,output='flextable')%>%
		theme_zebra()%>%
		set_header_labels('Model 1' = 'Method & Target\n\ Taxa GLMM')%>%
	    align(align = 'center', part = 'all')%>%
	    font(fontname = 'Arial', part = 'all')%>%
	    fontsize(size = 10, part = 'all')%>%
	    autofit()

	    ## fit 2
	    	
	    s2 <- readRDS('EDNA/data_edna/nimble_outputs/model_fit2_2000iter_500burn_3ch_worksButDoesntConverge_july2024.RDS')

		summary_fit2 <- MCMCsummary(s2,round=4,pg0=TRUE,prob=c(0.05,0.95)) %>%
		      tibble::rownames_to_column()%>%
		      rename_with(str_to_title)%>%
		      rename(Parameter = Rowname)%>%
		      rename('% of posterior with \n\ same sign as estimate' = 'P>0', Estimate = 'Mean','lower'='5%',upper='95%')%>%
		      mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
		      dplyr::select(-lower,-upper,-Sd)%>% 
		      flextable()%>%
		      theme_zebra()%>%
		      set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
		      align(align = 'center', part = 'all')%>%
		      font(fontname = 'Arial', part = 'all')%>%
		      color(color='black',part='all')%>%
		      fontsize(size = 10, part = 'all')%>%
		      autofit()
		    summary_fit2


	# Gather variable indices into a separate column (one row per draw, one column per parameter)
	## fit 1
		draws <- fit1 %>%
			gather_draws(`b_Intercept`, b_method.typewaterbottle, b_taxaShark,`b_method.typewaterbottle:taxaShark`, r_eventID[condition, ])%>%
			median_qi(.width = c(.95, .5)) # calculates the median and the 95% quantile interval for the parameter

		caterpillars <- draws %>%
			dplyr::filter(.variable != 'r_eventID') %>%
			ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper))+
			geom_pointinterval()+      
			geom_vline(xintercept=0,linetype=3)+
			theme_bw()
		caterpillars

		ggsave(caterpillars, file = 'EDNA/data_edna/brms/figures/caterpillars_model1_gaussian_log10data.png', device = 'png', unit = 'in', height = 4, dpi = 850)	

	## fit 2
		draws <- s2 %>%
			gather_draws(`b0`, `b1`, `b2`, `b3`, tau.re)%>%
			median_qi(.width = c(.95, .5))# calculates the median and the 95% quantile interval for the parameter

		caterpillars <- draws %>%
			ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper))+
			geom_pointinterval()+      
			geom_vline(xintercept=0,linetype=3)+
			theme_bw()
		caterpillars

		ggsave(caterpillars, file = 'EDNA/data_edna/brms/figures/caterpillars_model1_gaussian_log10data.png', device = 'png', unit = 'in', height = 4, dpi = 850)	



	# Add draws from posterior fit and (then) residuals, using epred_draws

	postfit <- data %>%
		dplyr::select('eventID', 'Sample.Name', 'taxa', 'method.type', 'log10.copies')%>%
		data_grid(taxa, method.type, eventID)%>%
		add_epred_draws(fit1)%>%
		ggplot(aes(x = .epred, y = interaction(taxa,method.type), fill = method.type))+
		stat_halfeye(alpha = 0.75)+
		scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
		geom_vline(xintercept=0,linetype=3)+
		labs(y = 'Interaction of Taxa & Method', x = 'Conditional Means')+
		guides(fill = 'none')+
		theme_bw()
	ggsave(postfit, file = 'EDNA/data_edna/brms/figures/postdraws_caterpillars_fixedeffects_model1_gaussian_log10data.png', device = 'png', unit = 'in', height = 5, dpi = 850)	


	## Diagnostics 

		## Overdispersion 
		n <- nrow(data)
			## fit1 
			res_fit1 <- resid(fit1)
			overdis_fit1 = sum(res_fit1^2)/(n-2)
			overdis_fit1 # 18.915, overdispersed (Gaussian on log10 copies)

		## Plots of residuals

			## fit 1
			res_fit1 <- as_tibble(resid(fit1)) # E1 in the book, page 58
			F1 <- as_tibble(predict(fit1))%>%rename(fitted = Estimate, fit.Q2.5 = Q2.5, fitQ97.5 = Q97.5, fitEst.Error = Est.Error)
			diag_fit1 <- bind_cols(F1, res_fit1)

			resVfit_fit1 <- ggplot(data = diag_fit1,aes(x= fitted, y = Estimate))+
				geom_point()+ 
				labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
				theme_bw()
			res_fit1_hist <- ggplot(data = res_fit1, aes(x = Estimate))+ 
				geom_histogram(binwidth = nrow(res_fit1)/1000,fill = 'black')+ 
				theme_bw()+
				labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
			res_fit1_qq <- ggplot(data=F1, aes(sample = fitted))+
				stat_qq(size=1,pch=21)+
				labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
				stat_qq_line(linetype=2, col='red')+
				theme_bw()

			diagnostics_fit1 <- resVfit_fit1+res_fit1_hist+res_fit1_qq
			
			ggsave(diagnostics_fit1, file = 'EDNA/data_edna/brms/figures/diagnostics_model1_gaussian_log10data.png', device = 'png', unit = 'in', height = 4, width = 8, dpi = 850)	





























