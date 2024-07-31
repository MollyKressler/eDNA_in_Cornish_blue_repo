## - Bayesian models for eDNA (2023 samples)


## created 5 July 2024
## by Molly M Kressler

########
## Load data  
########

pacman::p_load(sf,tidyverse,dplyr,ggplot2,flextable, modelsummary, readr, readxl, bayesplot, tidybayes, modelr, ggdist, nimble, MCMCvis)

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/')

data <- read.csv('EDNA/data_edna/qPCRresults/processedQPCRresults_cornwallspecies_june2024.csv')%>%
	filter(NTC.Amplified == 0)%>%
	dplyr::select(-NTC.Amplified)%>%
		as_tibble()%>%  # remove data from tests where th enegatie lab control amplified. 
		mutate(taxa = as.character(case_when(Target.Name == 'Engraulis' ~ 'Fish', Target.Name == 'Scomber' ~ 'Fish',Target.Name == 'Prionace' ~ 'Shark', Target.Name == 'Alopias' ~ 'Shark')))%>%
		mutate(log10.copies = log(copies.sampavg+0.1),
			round.copies = as.numeric(round(copies.sampavg,0)),
			re.level = as.numeric(as.numeric(factor(eventID))),
			methodID = as.numeric(as.numeric(factor(method.type))),
			taxaID = as.numeric(as.numeric(factor(taxa))))%>%
		dplyr::select('eventID', 're.level', 'Sample.Name', 'Target.Name','taxa', 'method.type', 'methodID', 'taxaID', 'copies','copies.techrepavg','log10.copies', 'round.copies')%>%
		group_by(Sample.Name, Target.Name)%>%
		slice(1) # nrow = 373, We are interested in the median technical replicate value, there are three tech reps per Sample per Target. The median value has already been determiend and put into a column that is repeated across tech reps, so we just need to slice the top of each group (Sample + Target), to get the metadata and the median value. 
	data # waterbottles = 2, metaprobes = 1,fish = 1, shark = 2. Adjust these in constants as Bayes/nimble indexes from 0. 

########
## Data exploration 
########

	hist<- ggplot(data = data, aes(x =log10.copies, fill = Target.Name))+
		geom_histogram(position = 'dodge', bins = 5)+
	  	scale_fill_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus'))+
		theme_bw()+
		labs(y = NULL, x = 'DNA yield (log10)')

	ggsave(hist, file = 'EDNA/data_edna/brms/figures/hist_bySpecies.png', device = 'png', units = 'in', height = 6, width = 6, dpi = 800)


########
## Nimble  models 
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

	