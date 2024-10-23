## - Bayesian models for eDNA (2023 samples)


## created 5 July 2024
## by Molly M Kressler

########
## Load data  
########

pacman::p_load(sf,tidyverse,dplyr,ggplot2,flextable, modelsummary, readr, readxl, bayesplot, tidybayes, modelr, ggdist, nimble, MCMCvis, patchwork)

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
		    muL[i] <- beta0 + beta1 * taxaID[i] + 
		                          beta2 * methodID[i] + 
		                          beta3 * taxaID[i] * methodID[i] + 
		                          e[re[i]]
		 }
		 # derived parameters #
             b0 <- exp(beta0) # method 0, taxa 0 - starting value against which others are compared: this is the population mean for metaprobe fish 
             b1 <- exp(beta0 + beta1) # method 1, taxa 0 - difference between metaprobe fish and wbt-fish, this is the effect of method
             b2 <- exp(beta0 + beta2) # method 0, taxa 1 - difference between metaprobe-fish and metaprobe-shark, this is the effect of taxa 
             b3 <- exp(beta0 + beta1 + beta2 + beta3) # method 1, taxa 1 - difference between meta-fish and wbt-shark, this is the effect of the interaction
             b4 <- exp(beta1 + beta3) # method 1, taxa 1 - difference between wbt-fish and wbt-shark, this is the effect of taxa for waterbottles
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
		  
		model_fitZ$calculate()  

		CmfZ<-compileNimble(model_fitZ) 
		confZ <- configureMCMC(model_fitZ, monitors = c('sigma.re', 'beta1','beta2', 'beta3', 'beta0', 'b1', 'b2', 'b3', 'b4', 'b0','pz', 'sdL'))#, onlySLICE = FALSE)
		MCMC_model_fitZ <- buildMCMC(confZ, na.rm = TRUE)
		ccMCMCZ <- compileNimble(MCMC_model_fitZ, project = model_fitZ)

		samplesZ <- runMCMC(ccMCMCZ,niter = 50000, nburnin = 12000,nchains = 3, samplesAsCodaMCMC = TRUE)

		MCMCsummary(samplesZ)		

		saveRDS(samplesZ, 'EDNA/data_edna/nimble_outputs_modelobjects/ZinflatedPoissonLogNormal_eDNA_50kiter_12kburn_3c_july2024.RDS')
	
########
## - Plots and Tables
#######

	# summary and caterpillars

		samplesListZ <- readRDS('EDNA/data_edna/nimble_outputs_modelobjects/ZinflatedPoissonLogNormal_eDNA_50kiter_12kburn_3c_july2024.RDS')

		draws <- samplesListZ %>%
			gather_draws(`b1`,`b2`, `b3`, `b4`, `b0`,`beta0`,`beta1`, `beta2`, `beta3`,)%>%
			median_hdi(.width = c(.95, .5)) # calculates the median and the 95% quantile interval for the parameter
		
		draws2 <- samplesListZ %>%
			gather_draws(`b1`,`b2`, `b3`, `b4`)%>%
			median_hdi(.width = c(.95, .5)) # calculates the median and the 95% quantile interval for the parameter

		caterpillars <- draws %>%
			filter(str_detect(.variable, 'eta'))%>%
			mutate(key = case_when(.variable == 'beta1' ~ 'Method', .variable == 'beta2' ~ 'Taxa', .variable == 'beta3' ~ 'Interaction', .variable == 'beta0' ~ 'Intercept'))%>%
			filter(key != 'Intercept')%>%
			mutate(key = fct_reorder(key,.value, .fun = 'median'))%>%
			ggplot(aes(y = key, x = .value, xmin = .lower, xmax = .upper))+
			geom_pointinterval()+
			geom_vline(xintercept=0,linetype=3)+
			theme_bw()+
			guides(col = 'none')+
			labs(x = 'Effect Size', y = NULL)
		caterpillars

		ggsave(caterpillars, file = 'EDNA/data_edna/figures_and_tables/nimble_outputs/caterpillars_modelZ_july2024.png', device = 'png', units = 'in', height = 4, width = 5, dpi = 1080)

		bayesplot_theme_set(theme_minimal(base_size = 8, base_family = "sans"))
		beta.traces <- mcmc_trace(samplesListZ, pars = c('beta0','beta1','beta2','beta3'))
		ggsave(beta.traces, file = 'EDNA/data_edna/figures_and_tables/nimble_outputs/betatraces_modelZ_july2024.png', device = 'png', units = 'in', height = 4, width = 5, dpi = 1080)

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

		maintext_fitZ <- MCMCsummary(samplesListZ,round=3,pg0=TRUE,prob=c(0.05,0.95)) %>%
		      tibble::rownames_to_column()%>%
		      rename_with(str_to_title)%>%
		      rename(Parameter = Rowname)%>%
			mutate(Parameter = case_when(Parameter == 'beta1' ~ 'Method',Parameter == 'beta0' ~ 'Intercept', Parameter == 'beta2' ~ 'Taxa', Parameter == 'beta3' ~ 'Interaction',Parameter == 'b0' ~ 'M-F', Parameter == 'b1' ~ 'W-F', Parameter == 'b2' ~ 'M-S', Parameter == 'b4' ~ 'W-S',  .default = NA), .before = 'Parameter')%>%
			filter(!is.na(Parameter))%>%
		      rename('pg0'='P>0')%>%
		      mutate(pg00 = case_when(Mean >= 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
		      rename('Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
		      mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
		      dplyr::select(-lower,-upper,-Sd, -pg0)%>%
		      arrange(-Estimate)%>%
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


	# plot mean estimates for effects 
		draws3 <- samplesListZ %>%
			gather_draws(`b0`,`b1`,`b2`, `b3`, `b4`)%>%
			median_hdi(.width = c(.5,.95))%>%
			mutate(key = case_when(.variable == 'b0' ~ 'M-F', .variable == 'b1' ~ 'W-F', .variable == 'b2' ~ 'M-S', .variable == 'b4' ~ 'W-S', .variable == 'b3' ~ 'Intx'))
		inter = as.numeric(draws3%>%filter(key == 'Intercept')%>%dplyr::select(.value))
		post <- draws3 %>%
			filter(key != 'Intx') # summary

		spread2 <- samplesListZ %>%
			spread_draws(`b0`,`b1`,`b2`, `b3`, `b4`)%>%
			pivot_longer(cols = c(`b0`,`b1`,`b2`, `b3`, `b4`), names_to = 'var', values_to = 'value')%>%
			mutate(key = case_when(var == 'b0' ~ 'M-F', var == 'b1' ~ 'W-F', var == 'b2' ~ 'M-S', var == 'b4' ~ 'W-S', var == 'b3' ~ 'Intx'),
				plot.order = case_when(var == 'b0' ~ 2, var == 'b1' ~ 1, var == 'b2' ~ 4, var == 'b4' ~ 3, var == 'b3' ~ 5
					))
		post2 <- spread2 %>%
			filter(key != 'Intx') # all draws 


		condition.means <-	ggplot()+
		 	geom_pointinterval(data = post, aes(y = .value, x = reorder(key, -.lower), ymax = .upper, ymin = .lower))+
		 	theme_bw()+
		 	labs(y = 'DNA yield (µg/µL)', x = 'Condition')+
		 	ylim(0,450)
		 condition.means


		lines <- post2 %>% 
			mutate(method = case_when(str_detect(key, 'M') ~ 'Metaprobe', str_detect(key,'W') ~ 'Water Bottle'))%>%
			group_by(method)%>%
			summarise(median = median(value))

		condition.means2 <- ggplot(data = post2, aes(y = value, x = reorder(key,plot.order)) )+
			stat_interval(geom = 'interval', .width = c(.5, .8, .95), point_interval = 'median_hdi', orientation = 'vertical', width = 0.9, stars = TRUE)+
		 	scale_color_grey(guide = NULL, start = 0.75,end = .15)+
		 	stat_summary(fun="median", geom="segment", mapping=aes(xend=..x.. - 0.1, yend=..y..), size = 1)+
		 	stat_summary(fun="median", geom="segment", mapping=aes(xend=..x.. + 0.1, yend=..y..), size = 1)+
		 	geom_hline(yintercept = lines$median[1], linetype = 2, col = '#DAA507')+ # 20.4
		 	geom_hline(yintercept = lines$median[2], linetype = 2, col = '#8EC7D2')+ # 40.5
		 	theme_bw()+
		 	labs(y = 'DNA yield', x = 'Condition', color = NULL)
	 	condition.means2 ## in the manuscript

		condition.means3 <-	ggplot()+
		 	stat_slabinterval(data = post2, aes(y = value, x = reorder(key, plot.order)), point_interval = 'median_hdi', position = 'dodge', fill = 'lightblue2', slab_type = 'pdf')+
		 	theme_bw()+
		 	labs(y = 'DNA yield (µg/µL)', x = 'Condition')+
		 	ylim(0,450)
		 condition.means3 # probability density function in blue

		ggsave(condition.means2, file = 'EDNA/data_edna/figures_and_tables/nimble_outputs/condition_means_CIsinBlue_medianLines_modelZ_july2024.png', device = 'png', units = 'in', height = 4, width = 5, dpi = 1080)


		pairwise.posterior.comparisons <- spread2 %>%
			compare_levels(mean, by = var)%>%
			ggplot(aes(x=mean,y = var))+
		 	stat_halfeye(.width = c(.5,.95), fill = 'lightblue2')+
		 	geom_vline(xintercept = 0, linetype = 2, lwd= 0.25)+
		 	theme_bw()+
		 	theme(panel.grid = element_blank())+
		 	labs(y = 'Pairs', x = 'Posteriors')

		ggsave(pairwise.posterior.comparisons, file = 'EDNA/data_edna/figures_and_tables/nimble_outputs/pairwise_posterior_comparisons_means_modelZ_july2024.png', device = 'png', units = 'in', height = 4, width = 5, dpi = 1080)
 
		both <- caterpillars + pairwise.posterior.comparisons + plot_layout(widths = c(.75,1)) + plot_annotation(tag_levels = 'a')
		ggsave(both, file = 'EDNA/data_edna/figures_and_tables/nimble_outputs/cater_AND_pairwise_posterior_comparisons_modelZ_july2024.png', device = 'png', units = 'in', height = 4, width = 8, dpi = 850)


	### Figure 7 a & b 
		fig7 <- caterpillars + condition.means2 + plot_layout(widths = c(.75,1)) + plot_annotation(tag_levels = 'a')

		ggsave(fig7, file = 'towardschapters_PhD/assembling the thesis/chp4_fig7_cater_AND_conditionmeans_modelZ_july2024.png', device = 'png', units = 'in', height = 4, width = 8, dpi = 1080)


	## Stack with cbind and remake conditional posteriors plot for Method. Plot the 50s+80s %s, leave out the 95%. 
		post2 ## spread long, so a cbind isn't what we need. We need to create a new ID based on 'var' which combines metaprobes and waterbottles groups.

		post3 <- post2 %>%
			mutate(group = case_when(
				var == 'b0' ~ 'Metaprobe',
				var == 'b1' ~ 'Water Bottle',
				var == 'b2' ~ 'Metaprobe',
				var == 'b4' ~ 'Water Bottle'
				), .after = 'var')

		stacked_groupmeans <- ggplot(data = post3, aes(y = value, x = group))+
			stat_interval(geom = 'interval', .width = c(.5, .8), point_interval = 'median_hdi', orientation = 'vertical', width = 0.9)+
		 	scale_color_grey(guide = NULL, start = 0.75,end = .15)+
		 	stat_summary(fun="median", geom="segment", mapping=aes(xend=..x.. - 0.08, yend=..y..), size = 1)+
		 	stat_summary(fun="median", geom="segment", mapping=aes(xend=..x.. + 0.08, yend=..y..), size = 1)+
		 	theme_bw()+
		 	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
		 	labs(y = 'Posterior DNA yield', x = NULL, color = NULL)
		
		ggsave(stacked_groupmeans, file = 'EDNA/data_edna/figures_and_tables/nimble_outputs/GlobalMeans_method_medianLines_modelZ_july2024.png', device = 'png', units = 'in', height = 4, width = 5, dpi = 1080)
		ggsave(stacked_groupmeans, file = 'EDNA/data_edna/figures_and_tables/nimble_outputs/GlobalMeans_method_medianLines_modelZ_july2024_forposter.png', device = 'png', units = 'in', height = 3.5, width = 3.5, dpi = 1080)















	