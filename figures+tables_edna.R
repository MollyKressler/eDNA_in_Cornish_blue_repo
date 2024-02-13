## - Figures and Tables for eDNA in Cornwall 
### non-research question specific, relating to qPCR results

## created 29 September 2023 
## by Molly M Kressler

# testign Git push
########
## Load data  
########
## cleaned in eDNA_datacleaning.R

pacman::p_load(sf,tidyverse,dplyr,ggplot2, patchwork, cowplot,lubridate,flextable)

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/')

## qPCR results, all species, not standards, test assay positive nor negative test controls 
sp1 <- st_as_sf(st_read('EDNA/data_edna/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.shp'),crs='WGS84')
sp2 <- read.csv('EDNA/data_edna/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv')
## qPCR standards, test assay positive and negative test controls 
c <- read.csv('EDNA/data_edna/resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv')# quick fix of label problem

########################################
## Sampling and Metadata Descriptive plots 
########################################

	method.palette <- c()

  ### Descriptive Stats
	## Samples per month - creates a bar chart of samples collected per month with a table showign the number of each type per month 
		samplespermonth <- ggplot(data=sp1,aes(x=month))+
			geom_histogram(stat='count',col='cadetblue4',fill='cadetblue4',alpha=0.8)+
			scale_x_discrete(labels=c('April','May','June','July','August','September'))+
			xlab('Month (2023)')+ ylab('Samples (count)')+
			guides(col='none',fill='none',alpha='none')+
			scale_y_continuous(breaks=c(3,6,9,12),limits=c(0,14))+
			theme_bw()
		methodtypepermonth_controlledformultiplerowsperday <- sp1%>%group_by(methodtype,month)%>%dplyr::select(methodtype,month)%>%tally()%>%ungroup()%>%as.tibble()%>%dplyr::select(-geometry)%>%pivot_wider(names_from=month,values_from=n)
		methodtypepermonth.table<-methodtypepermonth_controlledformultiplerowsperday%>%
			flextable()%>%
			set_header_labels('methodtype'='Method','4'='April','5'='May','6'='June','7'='July','8'='August','9'='September')%>%
			bold(bold=TRUE,part='header')%>%
			align(align = 'center', part = 'all')%>%
			font(fontname = 'Arial', part = 'all')%>%
			color(color='black',part='all')%>%
			fontsize(size = 10, part = 'all')%>%
			autofit()
		met <- flextable::gen_grob(methodtypepermonth.table)
		outreachANDvessels_permonth_apriltosept2023 <- cowplot::plot_grid(samplespermonth,met,nrow=2,ncol=1,rel_heights=c(6,1),scale=c(1,.75))
		ggsave(outreachANDvessels_permonth_apriltosept2023,file='EDNA/data_edna/figures_and_tables/samplesPERmonth_sampling_colourful_april2sept2023.png',device=png,units='in',height=5,width=5,dpi=600)
	
	## Outreach per month - creates a bar chart with number of passengers (unique) exposed to the project per month, and the number of participating (unique) vessel per month 

		passengerspermonth_controlledformultiplerowsperday <- sp1%>%group_by(eventID)%>%slice_head(n=1)%>%ungroup()%>%group_by(month)%>%tally(passengers)

		outreachpermonth<-ggplot(data=passengerspermonth_controlledformultiplerowsperday,aes(x=month,y=n))+
			geom_col(col='grey65',fill='grey65',alpha=0.8)+
			scale_x_discrete(labels=c('April','May','June','July','August','September'))+
			xlab('Month (2023)')+ ylab('Outreach (count)')+
			guides(col='none',fill='none',alpha='none')+
			theme_bw()
			# have to divide by two because the data contain 
		vesselspermonth_controlledformultiplerowsperday <- sp1%>%group_by(eventID)%>%slice_head(n=1)%>%ungroup()%>%group_by(vessel)%>%group_by(month)%>%tally()
		vesselspermonth.table<-vesselspermonth_controlledformultiplerowsperday%>%
			as.tibble()%>%
			dplyr::select(month,n)%>%
			pivot_wider(names_from=month,values_from=n)%>%
			flextable()%>%
			set_header_labels('4'='April','5'='May','6'='June','7'='July','8'='August','9'='September')%>%
			bold(bold=TRUE,part='header')%>%
			align(align = 'center', part = 'all')%>%
			font(fontname = 'Arial', part = 'all')%>%
			color(color='black',part='all')%>%
			fontsize(size = 10, part = 'all')%>%
			autofit()
		out <- flextable::gen_grob(vesselspermonth.table)
		outreachANDvessels_permonth_apriltosept2023 <- cowplot::plot_grid(outreachpermonth,out,nrow=2,ncol=1,rel_heights=c(6,1),scale=c(1,.75))
		ggsave(outreachANDvessels_permonth_apriltosept2023,file='EDNA/data_edna/figures_and_tables/outreach_and_vessels_permonth_april2september2023.png',device=png,units='in',height=5,width=5,dpi=600)

	## DNA concentration over time

		dnacont.overtime<-ggplot(data=d1,aes(x=sampling.date,y=dnaCont,group=methodtype,col=methodtype,fill=methodtype))+
				geom_smooth()+
				scale_fill_manual(values=c('cadetblue2','violetred4'))+
				scale_color_manual(values=c('cadetblue2','violetred4'))+
				scale_x_continuous(labels=c('April','May','June','July','August','September'))+
				scale_y_continuous(breaks=c(0,20,40,60,80,100))+
				xlab('Sampling Date (2023)')+ ylab('DNA Concentraton (double-stranded)')+
				theme_bw()
		ggsave(dnacont.overtime,file='dnaCont_oversamplingperiod_april2september2023.png',device=png,units='in',height=5,width=6.5,dpi=700)

  ### Maps 
	## LOCATIONS
		##  using sampling data only (more samples collected than processed at most times)
		locations_sampling <- ggplot()+
			geom_sf(data=kernios2,alpha=0.8,fill='grey72')+
			geom_sf(data=sp1, col='cadetblue4',size=2,pch=18,alpha=0.45)+
			theme_bw() 
			ggsave(locations_sampling,file='locations_sampling_april2sept2023.png',device=png,units='in',height=5,width=5,dpi=600)
		
		locations_sampling_CSorR <- ggplot()+
			geom_sf(data=kernios2,alpha=0.8,fill='grey72')+
			geom_sf(data=sp1%>%filter(type=='R'), col='cadetblue4',size=2,pch=18,alpha=0.45)+
			geom_sf(data=sp1%>%filter(type=='CS'), col='violetred4',size=2,pch=16)+
			theme_bw() 

		locations_sampling_bymonth <- ggplot()+
			geom_sf(data=kernios2,alpha=0.8,fill='grey72')+
			geom_sf(data=sp1, aes(colour=month),size=3,pch=18,alpha=0.75)+
			scale_colour_manual(values=sixmonths.palette)+
			theme_bw() 

		##  using extraction + metadata sf
		location_extractions <- ggplot()+
			geom_sf(data=kernios2,col='grey72')+
			geom_sf(data=sp1, col='cadetblue4',size=2,pch=18,alpha=0.45)+
			theme_bw() 

	# Cq.means of NTCs and samples
	c



########################################
## qPCR results figures & tables 
########################################

## examine NTCs and field controls for amplification 
	head(c)
	amplified <- c%>%filter(Cq!='NA') # includes Standards 
	nrow(amplified)

## import species results (copies per replicate sample)

	copies <- read_csv("data_downloads_notpublic/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv")
	ntc <- read_csv("data_downloads_notpublic/resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv", 
	                col_types = cols(...1 = col_skip(), X = col_skip()))%>%
	  mutate(Assay.Role = case_when(Assay.Role == 'Negative' ~ 'NTC',
	                                Assay.Role == 'Positive' ~ 'Standard'
	                                , .default = as.character(Assay.Role)))
	
	##########
	## - Separate the Field controls from the copies results
	##########
	## use 'copies', replicateID will include 'WBC' 
	
	field.controls <- copies %>% filter(grepl('WBC', replicateID)) %>%
	  mutate('Assay.Role' = 'Field.Control') # This is the field controls for the water bottles (and metaprobes from joint deployments)
	field.controls 
	
	samples <- copies %>% filter(!grepl('WBC', replicateID)) %>%
	  mutate('Assay.Role' = 'Field.Sample') # This is all field samples, water bottles and metaprobes, without the field controls 
	
	eng_sf <- st_as_sf(st_read('compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.shp'),crs='WGS84')%>% # spatial and meta 
		rename(sampleID = samplID,
			Sample.Name = Smpl_Nm,
			replicateID = rplctID,
			engraulis.cq.mean = engrl__,
			engraulis.copies = engrls_,
			sampDATE = smpDATE,
			methodtype = mthdtyp)
	
	##########
	## - Visualizing Amplification Results
	##########
	## stack datasets of field samples, field controls, and lab controls, & select down to only the columns we need - testID, Assay.Role, Sample.Name, species.cq.mean, species.copies
	## engraulis 
	eng.samples2 <- samples %>% dplyr::select(eventID,sampleID,testID, Assay.Role, Sample.Name, engraulis.cq.mean,engraulis.copies,methodtype)
	eng.field.controls2 <- field.controls %>% 
		dplyr::select(eventID,sampleID,testID, Assay.Role, Sample.Name, engraulis.cq.mean,engraulis.copies,methodtype)
	eng.ntc2 <- ntc %>% 
		dplyr::select(testID, Assay.Role, Sample.Name, Cq.Mean,engraulis.copies)%>%
		rename(engraulis.cq.mean=Cq.Mean)%>%
		mutate(sampleID = 'lab.controls',eventID = 'lab.controls',methodtype = 'lab.controls')
	
	eng.f2 <- bind_rows(eng.samples2, eng.field.controls2, eng.ntc2)
	
	## where Cq.Mean is zero matters but these wont show on a plot showing the amplification of samples and standards and NTCs - so we add in a false amplification one cycle LONGER than it ran for. 
	
	eng.ntcs.na <- eng.f2%>%
	  filter(is.na(engraulis.cq.mean) & Assay.Role == 'NTC')%>%
	  mutate(forplot = 43)
	eng.field.sample.nas <- eng.f2%>%
	  filter(is.na(engraulis.cq.mean) & Assay.Role == 'Field.Sample')%>%
	  mutate(forplot = 43)
	
	#### PLOT: Cq.value/mean by type (NTC, Standard, Field Sample, Field Control) - species specific
	## plot them as dots with colors based on Assay.Role
	cqmeans_allwells_engraulis <- ggplot()+
	       geom_point(data=eng.f2%>%filter(engraulis.cq.mean !='NA'), aes(x = as_factor(testID), y = engraulis.cq.mean, col = Assay.Role), size=3)+
	      scale_color_manual(values = c('#DAA507','#8EC7D2','#733A4F','#07475A'))+
	      geom_point(data=eng.ntcs.na, aes(x = as_factor(testID), y= forplot), size=5, pch=8, stroke = 1.5, col='#733A4F')+
	    geom_point(data=eng.field.sample.nas, aes(x = as_factor(testID), y= forplot), size=2, pch=8, col='#8EC7D2', position = 'jitter')+
	     ylim(0,50)+
	      ylab('CQ (mean)')+
	      xlab('Assay')+ 
	      theme_bw()
    
		cqmeans_allwells_engraulis
		ggsave(cqmeans_allwells_engraulis,file = 'EDNA/data_edna/figures_and_tables/cqmeans_allassays_engraulis_plot.png', device=png,units='in',height=5,width=6.5,dpi=700)



  ## Species standard curve test results
		  ## engraulis
	


  ## Glmer model of copy number between methods type, random effect of event ID
    lm1 <- readRDS('data_edna/modelRDS/glmer_RandomEventID_allspecies_copies_byMethodType.RDS') ## all species
    lm2 <- readRDS('data_edna/modelRDS/glmer_RandomEventID_fish_copies_byMethodType.RDS') ## fish 
    lm3 <- readRDS('data_edna/modelRDS/glmer_RandomEventID_sharks_copies_byMethodType.RDS') ## sharks   

    models <- c(lm1, lm2, lm3)

    summary <- modelsummary(models, coef_rename = c('methodtypewaterbottle' = 'Water Bottle', 'methodmetaprobe' = 'Metaprobe'),fmt=3,estimate='estimate', statistic='conf.int',stars=FALSE,conf_level=0.95,output='flextable')

    summary_flex <- summary %>%
        set_header_labels('Model 1' = 'All Species', 'Model 2' = 'Fish Species', 'Model 3' = 'Shark Species')%>%
        theme_zebra()%>%
        align(align = 'center', part = 'all')%>%
        font(fontname = 'Arial', part = 'all')%>%
        fontsize(size = 10, part = 'all')%>%
        autofit()

    save_as_image(summary_flex, 'summary_glmers_eventIDrandom_copies_byMethod.png', webshot = 'webshot2')

     

  ## Mega Table, grouped by Event ID: replicate copies and CQ results, corresponding experiment TestID, NTCs and Standards CQ values. 









########
## - Visualising Species Detections and copy number
########

## DENSITY OF DETECTIONS, method specific then species specfic (so 12 plots in total) - creates a map plot with points of sampling locations (grouped by sampleID), where the size of the point ~ the relative abundance of a species  
	## data = extraction data frame.



## SPECIES PROPORTION OF DETECTIONS, method specific (2 plots in total) - creates a map plot with points of sampling locations (grouped by sampleID), where the point is a pie chart, and each species is a different colour, and the slice width is ~ the relative abundance of a species at THAT location. 
















