## - Figures and Tables for eDNA in Cornwall 

## created 29 September 2023 
## by Molly M Kressler


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

	engraul <- read.csv('EDNA/data_edna/copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv') # replicates and copies results 
	allspp <- read.csv('copies_perReplicate_notStandards_allspecies_2023eDNACornishBlue.csv')
	summary(engraul)
	eng_sf <- st_as_sf(st_read('compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.shp'),crs='WGS84')%>% # spatial and meta 
		rename(sampleID = samplID,
			Sample.Name = Smpl_Nm,
			replicateID = rplctID,
			engraulis.cq.mean = engrl__,
			engraulis.copies = engrls_,
			sampDATE = smpDATE,
			methodtype = mthdtyp)

## range and medians Cq.means for samples, (1) all and (2) 'real' detections

	# stack the data sets - select down for testID, Sample.Name, Cq.Mean, replicateID, engraulils.copies
	samples <- engraul %>%
		dplyr::select(testID, Assay.Role, Sample.Name, Cq.Mean,engraulis.copies)
	c2 <- c %>%
		dplyr::select(testID, Assay.Role, Sample.Name, Cq.Mean, engraulis.copies)
	full <- bind_rows(samples,c2)
	# NTC and Negataive are the same but are called different sometimes. And a standard is a positive. 
	f2 <- full %>%
		mutate(Assay.Role = case_when(Assay.Role == 'Negative' ~ 'NTC',
			Assay.Role == 'Positive' ~ 'Standard'
			, .default = as.character(Assay.Role)))
	f2
	# where Cq.Mean is zero matters but these wont show ona plot shoing the amplifcation of samples and standards and NTCs 

	ntcs_na <- f2%>%filter(is.na(Cq.Mean) & Assay.Role == 'NTC')%>%
		mutate(forplot = 42)

	## some samples have Cq.Mean == NA..do the same for them, and add a jigger effect
	
	samples_na <- f2%>%filter(is.na(Cq.Mean) & Assay.Role == 'Unknown')%>%
		mutate(forplot = 42)

		scale_color_manual(values = c('#733A4F','#DAA507','#8EC7D2','#0D6A87','#07475A'))+

	## plot them as dots with colors based on Assay.Role
	cqmeans_allwells_engraulis <- ggplot()+
		geom_point(data=f2%>%filter(Cq.Mean !='NA'), aes(x = as_factor(testID), y = Cq.Mean, col = Assay.Role), size=3)+
		scale_color_manual(values = c('#733A4F','#DAA507','#8EC7D2'))+
		geom_point(data=ntcs_na, aes(x = as_factor(testID), y= forplot), size=5, pch=8, stroke = 1.5, col='#733A4F')+
		geom_point(data=samples_na, aes(x = as_factor(testID), y= forplot), size=2, pch=8, col='#8EC7D2', position = 'jitter')+
		ylim(0,50)+
		ylab('CQ (mean)')+
		xlab('Assay')+ 
		theme_bw()
	cqmeans_allwells_engraulis

		ggsave(cqmeans_allwells_engraulis,file = 'EDNA/data_edna/figures_and_tables/cqmeans_allassays_engraulis_plot.png', device=png,units='in',height=5,width=6.5,dpi=700)



## tables of species standard curve test results
	## table per species of standards results (6 tables)

	## table per species of replicates Cq and copies estimates (6 tables)


## copies per replicate - species reuslts stacked (1 plot)
	## 







## DENSITY OF DETECTIONS, method specific then species specfic (so 12 plots in total) - creates a map plot with points of sampling locations (grouped by sampleID), where the size of the point ~ the relative abundance of a species  
	## data = extraction data frame.



## SPECIES PROPORTION OF DETECTIONS, method specific (2 plots in total) - creates a map plot with points of sampling locations (grouped by sampleID), where the point is a pie chart, and each species is a different colour, and the slice width is ~ the relative abundance of a species at THAT location. 





##################
## Investigating values of field controls - n.4s 
	a <- sp2 %>% 
		filter(methodtype == 'waterbottle')%>%
		filter(grepl('.4',Sample.Name)) # grabs the rows  that have the string '.4' in the Sample.Name column
	ggplot(data=c%>%filter(Assay.Role == c('Negative')), aes(x=Experiment, y=Cq.Mean))+
	geom_bar(stat='identity')+ 
	theme_bw()





















