## - Analysis and results for eDNA in Cornwall 


## created 29 September 2023 
## by Molly M Kressler

########
## Load data  
########
## cleaned in pipeline documents

pacman::p_load(sf,tidyverse,dplyr,ggplot2, patchwork, cowplot,lubridate,flextable, rnaturalearth)

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/')

## qPCR results, all species, not standards, test assay positive nor negative test controls 
fieldsamples.sf <- st_as_sf(st_read('EDNA/data_edna/qPCRresults/processedQPCRresults_cornwallspecies_june2024.shp'),crs='WGS84')%>%
	rename(Assay.Role = Assy_Rl, Target.Name = Trgt_Nm, sampleID = samplID, ratio260.280 = r260_28, ratio260.230 = r260_23, sampledate = sampldt, method.subsample = mthd_sb, sampling.date = smplng_, recorded.by = rcrdd_b, passengers = pssngrs, time.timeIN = tm_tmIN, method.type = mthd_ty, Reporter = Reportr, Quencher = Quenchr, Sample.Name = Smpl_Nm, intercept = intrcpt, r.squared = r_squrd, efficiency = effcncy, loq_check = lq_chck, reliable = reliabl, copies.techrepavg = cps_tch, copies.sampavg = cps_smp, Well.Position = Wll_Pst)
fieldsamples <- read.csv('EDNA/data_edna/qPCRresults/processedQPCRresults_cornwallspecies_june2024.csv')

sp3 <- sp1 %>% distinct(Sample.Name, .keep_all = TRUE) %>%
		dplyr::select(Sample.Name, eventID:geometry)

## qPCR standards, test assay positive and negative test controls 
controls <- read.csv('EDNA/data_edna/qPCRresults/processedqPCRresults_2024_labANDfield_controls_standards_cornwallspecies.csv')

## land shapes
coastalwater <- st_as_sf(st_read('EDNA/data_edna/kernow_ios_coastalwater.shp'), crs='WGS84')
kernow <- st_as_sf(st_read('EDNA/data_edna/kernow.shp'), crs='WGS84')
ios <- st_as_sf(st_read('EDNA/data_edna/islesofscilly.shp'), crs='WGS84')
kernios <- st_as_sf(st_union(kernow, ios), crs='WGS84')
kernios <- st_transform(kernios, crs='WGS84')
## inset globe 
uk <- ne_countries(scale='medium', returnclass = 'sf') %>% dplyr::select(sovereignt,region_un, region_wb, subregion,name_sort,abbrev,admin, geometry) %>% filter(admin == 'United Kingdom')
uk
	# ggplot()+geom_sf(data = uk, fill = 'grey72', col='black')


#########
## - Locations 
#########

	location_extractions <- ggplot()+
		geom_sf(data=coastalwater,fill='cadetblue4', lwd=0.5, alpha = 0.2)+
		geom_sf(data=kernios,fill='grey82', lwd=0.5)+
		geom_sf(data=st_jitter(st_geometry(sp3), factor=0.01),size=2,pch=20, alpha = 0.5)+
		theme_bw() 
	#location_extractions # this is at the sample replicate level

	uk.map <- ggplot()+
		geom_sf(data = uk)+
		geom_rect(aes(xmin=-6.5,xmax=-4,ymin=49.5,ymax=51),col='deeppink4',fill='deeppink2',alpha=0.4,lwd=.8)+
		theme_void()+
		theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(), panel.border = element_rect(fill=NA, colour = 'black'), panel.background = element_rect(fill = 'white', color = NA))
	uk.map

	sample_location_extractions <- ggplot()+
		geom_sf(data=coastalwater,fill='cadetblue4', col=NA, lwd=0.5, alpha = 0.5)+
		geom_sf(data=kernios, lwd=0.5)+
		geom_sf(data=st_jitter(st_geometry(sp3), factor=0.0075), col='#000A52',size=2,pch=20, alpha = 0.5)+
		theme_bw()
	#sample_location_extractions # this is at the sample replicate level

	map <- ggdraw()+
		draw_plot(sample_location_extractions)+
		draw_plot(uk.map, height=0.4, x = -0.29, y = 0.52)
	
	ggsave(map,file = 'EDNA/data_edna/figures_and_tables/samplinglocations_2023_withUKinset.png',device='png',units='in',dpi=450,height=5.5,width=5.5)

	## Distance from shore 

	sp11 <- sp1 %>% mutate(dist2shore.km = st_distance(sp1, kernios)/1000)
	min(sp11$dist2shore.km) # 0.03 km
	median(sp11$dist2shore.km) # 6.149 km
	max(sp11$dist2shore.km) # 50.585 km


#########
## - qPCR results - standard curve test results
#########

	neb <- read_excel('qPCRresults/standard_curve_equations_qPCR_2023species.xlsx', sheet = 'Sheet1', range = 'A1:I8',.name_repair = 'universal', col_types = c('text','text','numeric','numeric','numeric','numeric','date','numeric','numeric'))%>%
		mutate(dateofassay=as_date(dateofassay))%>%
		filter(year(dateofassay) == 2024) 
	neb

	neb.table <- neb %>% 
		dplyr::select(-sp)%>%
		flextable()%>%
		autofit()%>%
		theme_zebra()

		 save_as_docx(neb.table, path = 'qPCRresults/cornwallednaspecies_standardcurvetestresults_2024.docx') 
		 save_as_image(neb.table, 'qPCRresults/cornwallednaspecies_standardcurvetestresults_2024.png', webshot = 'webshot2')


#########
## - qPCR results - amplification plots, by species
#########
## done in a separate R doc because it takes up a lot of lines. 


#########
## - qPCR results - amplification reliability, table 
#########

	head(fieldsamples)
	head(controls)

	# select down to only what we need

	fs <- fieldsamples %>% dplyr::select('Assay.Role','Target.Name','Sample.Name', 'eventID', 'method.type', starts_with('copies.'), 'testID')%>%as_tibble()
	c <- controls %>% dplyr::select('Assay.Role','Sample.Name', 'Target.Name','eventID', 'method.type', starts_with('copies.'), 'testID')%>%as_tibble()%>%
		filter(Assay.Role == 'UNKNOWN')

	st <- controls %>% dplyr::select('Assay.Role','Sample.Name', 'Target.Name','eventID', 'method.type', starts_with('copies.'), 'testID', Quantity)%>%as_tibble()%>%
		filter(Assay.Role == 'STANDARD')

	# calculate the number of tech reps that amplified

	t.fs <- fs %>% group_by(Sample.Name, Target.Name) %>% 
        tally(copies.techrepavg != 'NA') #404 rows
    
    fs2 <- left_join(fs%>%group_by(Sample.Name, Target.Name), t.fs, relationship = 'many-to-one', by=c('Sample.Name', 'Target.Name'))%>%
    	ungroup()%>%
      mutate(Amplified = paste0(n,'/3'))%>%
      mutate_if(is.numeric, round, digits = 3)%>%
      arrange('EventID')%>%
      group_by(Sample.Name, Target.Name)%>%
      slice(1)%>%
      dplyr::select(-n)
      fs2


	t.c <- c %>% group_by(Sample.Name, Target.Name) %>% 
        tally(copies.techrepavg != 'NA')
    
    c2 <- left_join(c%>%group_by(Sample.Name, Target.Name), t.c, relationship = 'many-to-one', by=c('Sample.Name', 'Target.Name'))%>%
    	ungroup()%>%
      mutate(Amplified = paste0(n,'/3'))%>%
      mutate_if(is.numeric, round, digits = 3)%>%
      group_by(Sample.Name, Target.Name)%>%
      slice(1)%>% 
      arrange('EventID')%>%
      dplyr::select(-n)
      c2

	t.st <- st %>% group_by(testID, Target.Name, Quantity) %>% 
        tally(copies.techrepavg != 'NA')
    
    st2 <- left_join(st%>%group_by(testID, Target.Name, Quantity), t.st, relationship = 'one-to-one', by=c('testID', 'Target.Name', 'Quantity'))%>%
    	ungroup()%>%
      mutate(Amplified = paste0(n,'/1'))%>%
      mutate_if(is.numeric, round, digits = 3)%>%
      arrange('EventID')%>%
      dplyr::select(-n)
      st2

    # make table 

      allampRatio <- bind_rows(fs2, c2, st2) %>%
      	dplyr::select(Sample.Name, Target.Name, Assay.Role, eventID, method.type, testID, copies.sampavg, Amplified)%>%
      	flextable()%>%
      	autofit()%>%
      	theme_zebra()
      allampRatio # for supplmentary


        save_as_image(allampRatio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_all.png', webshot = 'webshot2')     	
       	save_as_docx(allampRatio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_all.docx', webshot = 'webshot2')   

      onlyamp_ampratio <- bind_rows(fs2, c2, st2) %>%
      	filter(copies.sampavg != 'NA')%>%
      	dplyr::select(Sample.Name, Target.Name, Assay.Role, eventID, method.type, testID, copies.sampavg, Amplified)%>%
      	arrange(Sample.Name)%>%      	
       	flextable()%>%
       	set_header_labels('Sample.Name' = 'Replicate','Target.Name' = 'Target', 'method.type' = 'Method', Assay.Role = 'Assay Role', eventID = 'Sampling \n\ Event ID', testID = 'qPCR \n\ Assay ID', copies.sampavg = 'Replicate \n\ Avg. Copy No.')%>%
      	autofit()%>%
      	theme_zebra()%>% 
      	align(align = 'center', part ='all')

       	save_as_image(onlyamp_ampratio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_onlyAmplifiedones.png', webshot = 'webshot2')     	
       	save_as_docx(onlyamp_ampratio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_onlyAmplifiedones.docx')


#########
## - Comparing methods - reliable amplifications
#########

      ## How many metaprobes/waterbottles had 3/3 of amplifications

		a <- fs2 %>%
      	filter(copies.sampavg != 'NA')
      	aa <- a %>% group_by(Sample.Name, method.type)%>%
      	count()%>%
      	group_by(method.type)%>%
      	count()
      	aa # 40 metaprobes field replicates amplified; 48 waterbottle field replicates amplified, irrespective of species
      	
      	aaa <- a %>% group_by(Target.Name,method.type)%>%
      	count()
    	aaa.table<-aaa %>%      	
       	flextable()%>%
       	set_header_labels('Target.Name' = 'Target', 'method.type' = 'Method', n = 'No. of\n\ 3/3 Amp.' )%>%
      	autofit()%>%
      	theme_zebra()%>% 
      	align(align = 'center', part ='all')

      	save_as_image(aaa.table, path ='EDNA/data_edna/figures_and_tables/no_of_replicates_per_method_amplified3of3.png', webshot = 'webshot2')


#########
## - Comparing methods - concentration over time
#########








#########
## - Comparing methods - Copies over time: all, fish grouped, sharks grouped
#########







#########
## -  Glmer (stat. tests): effect of method on copies over time for all species and tazon groups
#########


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
    autofit() # this table goes with the geom_smooth graphs of copies over time by taxon

save_as_image(summary_flex, 'EDNA/data_edna/qPCRresults/figures_and_tables/summary_glmers_eventIDrandom_copies_byMethod.png', webshot = 'webshot2')

  








#########
## -  Effect of Soak Time 
#########

  soaks <- samples %>%
    filter(methodtype=='metaprobe')%>%
    mutate(timeOUT = case_when(timeOUT == 'TBD' ~ as.character(timeIN), is.na(timeOUT) ~ as.character(timeIN), .default = as.character(paste0(timeOUT,':00'))))%>%
  mutate(timeOUT2 = hms(timeOUT), timeIN2 = hms(timeIN), soaktime = timeOUT2-timeIN2)%>%
    dplyr::select(Sample.Name, replicateID, eventID, soaktime)%>%
    mutate(soaktime = seconds(soaktime))
 soaks.flex <-  soaks %>% flextable() %>%
    set_header_labels(Sample.Name = 'Sample ID', replicateID = 'replicateID', 'eventID' = 'Event ID', soaktime = 'Soak Time (sec)')%>%
    theme_zebra()%>%
    align(align = 'center', part = 'all')%>%
    font(fontname = 'Arial', part = 'all')%>%
    fontsize(size = 10, part = 'all')%>%
    autofit()
    
  save_as_image(soaks.flex, 'soaktimes_metaprobes_seconds_2023.png', webshot = 'webshot2')














