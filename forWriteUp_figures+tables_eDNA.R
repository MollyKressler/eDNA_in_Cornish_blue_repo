## - Analysis and results for eDNA in Cornwall 


## created 29 September 2023 
## by Molly M Kressler

########
## Load data  
########
## cleaned in pipeline documents

pacman::p_load(sf,tidyverse,dplyr,ggplot2, patchwork, cowplot,lubridate,flextable, rnaturalearth, lme4, modelsummary, readr, readxl)

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/')

## qPCR results, all species, not standards, test assay positive nor negative test controls 
fieldsamples.sf <- st_as_sf(st_read('EDNA/data_edna/qPCRresults/processedQPCRresults_cornwallspecies_june2024.shp'),crs='WGS84')%>%
		rename(Assay.Role = Assy_Rl, Target.Name = Trgt_Nm, sampleID = samplID, ratio260.280 = r260_28, ratio260.230 = r260_23, sampledate = sampldt, method.subsample = mthd_sb, sampling.date = smplng_, recorded.by = rcrdd_b, passengers = pssngrs, time.timeIN = tm_tmIN, method.type = mthd_ty, Reporter = Reportr, Quencher = Quenchr, Sample.Name = Smpl_Nm, intercept = intrcpt, r.squared = r_squrd, efficiency = effcncy, loq_check = lq_chck, reliable = reliabl, copies.techrepavg = cps_tch, copies.sampavg = cps_smp, Well.Position = Wll_Pst, NTC.Amplified = NTC_Amp)%>%
		filter(NTC.Amplified == 0)%>%
		dplyr::select(-NTC.Amplified)

fieldsamples <- read.csv('EDNA/data_edna/qPCRresults/processedQPCRresults_cornwallspecies_june2024.csv')%>%
	filter(NTC.Amplified == 0)%>%
	dplyr::select(-NTC.Amplified) # st_write converted the TRUE/FALSE to 0s and 1s. 0 = TRUE and 1 = FALSE. This removes any plate where the NTC amplified before the LOQ. 

sp3 <- fieldsamples.sf %>% distinct(Sample.Name, .keep_all = TRUE) %>%
		dplyr::select(Sample.Name, eventID:geometry)

## qPCR standards, test assay positive and negative test controls 
controls <- read.csv('EDNA/data_edna/qPCRresults/processedqPCRresults_2024_lab_controls_standards_formattedForTables_cornwallspecies.csv')

## glmer RDS

lm1 <- readRDS('EDNA/data_edna/glmer1_methodtype_copiesTechRep_intxTaxa_REeventid_Poisson.RDS')
lm2 <- readRDS('EDNA/data_edna/glmer2_soaktime_copiesTechRep_REeventid.RDS')

## land shapes
coastalwater <- st_as_sf(st_read('EDNA/data_edna/kernow_ios_coastalwater.shp'), crs='WGS84') %>% st_transform(crs = 'WGS84')
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

	sp11 <- fieldsamples.sf %>% mutate(dist2shore.km = st_distance(fieldsamples.sf, kernios)/1000)
	min(sp11$dist2shore.km) # 0.03 km
	median(sp11$dist2shore.km) # 6.149 km
	max(sp11$dist2shore.km) # 50.585 km


#########
## - qPCR results - standard curve test results
#########

	neb <- read_excel('EDNA/data_edna/qPCRresults/standard_curve_equations_qPCR_2023species.xlsx', sheet = 'Sheet1', range = 'A1:I8',.name_repair = 'universal', col_types = c('text','text','numeric','numeric','numeric','numeric','date','numeric','numeric'))%>%
		mutate(dateofassay=as_date(dateofassay))%>%
		filter(year(dateofassay) == 2024) 
	neb

	neb.table <- neb %>% 
		dplyr::select(-sp)%>%
		flextable()%>%
		autofit()%>%
		theme_zebra()

		 save_as_docx(neb.table, path = 'EDNA/data_edna/qPCRresults/cornwallednaspecies_standardcurvetestresults_2024.docx') 
		 save_as_image(neb.table, 'qPCRresults/cornwallednaspecies_standardcurvetestresults_2024.png', webshot = 'webshot2')


#########
## - qPCR results - amplification plots, by species
#########
  
  ## load data for results and amplifications, join parts to provide Sample.Name and Task
	eng1 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/05062024-KRESSLER-ENGRAULIS-ASSAY1.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp1 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/05062024-KRESSLER-ENGRAULIS-ASSAY1.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., eng1%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')


	eng2 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/05062024-KRESSLER-ENGRAULIS-ASSAY2.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp2 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/05062024-KRESSLER-ENGRAULIS-ASSAY2.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., eng2%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
	eng3 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-ENGRAULIS-ASSAY3.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp3 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-ENGRAULIS-ASSAY3.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., eng3%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	

	sco1 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-SCOMBER-ASSAY4.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp4 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-SCOMBER-ASSAY4.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., sco1%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')


	sco2 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-SCOMBER-ASSAY5.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp5 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-SCOMBER-ASSAY5.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., sco2%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
	sco3 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/11062024-KRESSLER-SCOMBER-ASSAY6.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp6 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/11062024-KRESSLER-SCOMBER-ASSAY6.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., sco3%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
	
	prio1 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-PRIONACE-ASSAY7.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp7 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-PRIONACE-ASSAY7.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., prio1%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')


	prio2 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY8.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp8 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY8.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., prio2%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
	prio3 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY9.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp9 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY9.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., prio3%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
	
	alo1 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-ALOPIAS-ASSAY10.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp10 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-ALOPIAS-ASSAY10.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., alo1%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')

	alo2 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY11.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp11 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY11.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., alo2%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
	alo3 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY12.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp12 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY12.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., alo3%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
	
	combo1 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY17.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp13 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY17.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., combo1%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
	
	combo2 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY18.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp14 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY18.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
	left_join(., combo2%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
	
	combo3 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY19.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types=c('numeric','text','text', 'text','text','text','text','text','text','numeric','numeric','numeric'))

	amp15 <- read_excel('EDNA/data_edna/qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY19.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal', col_types=c('numeric', 'numeric', 'text', 'numeric', 'numeric'))%>%
		left_join(., combo3%>%dplyr::select(Well, Sample.Name, Task), by = 'Well')
	
  ## bind_rows of amplification data, filter for only 'Unknown'/Field samples and controls
	
	amps <- bind_rows(amp1,amp2,amp3,amp4,amp5,amp6,amp7,amp8,amp9,amp10,amp11,amp12,amp13,amp14,amp15, .id = 'assay')%>% 
		mutate(assay = as.numeric(assay))%>%
		filter(Task == 'UNKNOWN')
  
  ## plot & save

	amps.samples.byassay <- ggplot(amps, aes(x = Cycle, y = Delta.Rn, col = Sample.Name))+
	  geom_smooth(span = 0.3, method = 'gam', se = FALSE, lwd = 0.4) +
	  theme_bw() +
	  theme(legend.position = 'none')+
	  scale_colour_viridis_d()+
	  facet_wrap(~assay, scales = 'free', dir = 'h', labeller = labeller(assay = c(
	  	'1' = 'Assay 1, all E. encrasicolus', 
	  	'2' = 'Assay 2, all E. encrasicolus',
	  	'3' = 'Assay 3, all E. encrasicolus',
	  	'4' = 'Assay 4, all S. scombrus',
	  	'5' = 'Assay 5, all S. scombrus',
	  	'6' = 'Assay 6, all S. scombrus',
	  	'7' = 'Assay 7, all P. glauca',
	  	'8' = 'Assay 8, all P. glauca',
	  	'9' = 'Assay 9, all P. glauca',
	  	'10' = 'Assay 10, all A. vulpinas',
	  	'11' = 'Assay 11, all A. vulpinas',
	  	'12' = 'Assay 12, all A. vulpinas',
	  	'13' = 'Assay 13, combination\n\ of species',
	  	'14' = 'Assay 14, combination\n\ of species',
	  	'15' = 'Assay 15, combination\n\ of species'
	  	)))+
      theme(strip.text = element_text(face = "italic"),strip.background = element_rect(fill = NA,colour = NA),plot.title = element_text(size=10, face='italic'))
	ggsave(amps.samples.byassay, file = 'EDNA/data_edna/figures_and_tables/fieldsample_andcontrols_Amplifications_by_assay.png', device = 'png', units = 'in', height = 8, width = 11, dpi = 850)






#########
## - qPCR results - amplification reliability, table 
#########

	head(fieldsamples)
	head(controls)

	# select down to only what we need

	fs <- fieldsamples %>% dplyr::select('Assay.Role','Target.Name','Sample.Name', 'eventID', 'method.type', starts_with('copies.'), 'testID', 'Cq','reliable')%>%as_tibble()
	c <- fieldsamples %>% dplyr::select('Assay.Role','Sample.Name', 'Target.Name','eventID', 'method.type', starts_with('copies.'), 'testID', 'Cq','reliable')%>%as_tibble()%>%
		filter(grepl("\\.4$", Sample.Name))

	st <- controls %>% dplyr::select('Assay.Role','Sample.Name', 'Target.Name','eventID', 'method.type', starts_with('copies.'), 'testID', 'Cq', Quantity,'reliable')%>%as_tibble()%>%
		filter(Assay.Role != 'UNKNOWN')

	# calculate the number of tech reps that amplified

	t.fs <- fs %>% group_by(Sample.Name, Target.Name) %>% 
        tally(Cq != 'NA' & reliable =='0') 
    
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
        tally(Cq != 'NA' & reliable =='0') 
    
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
        tally(Cq != 'NA' & reliable =='TRUE') 
    
    st2 <- left_join(st%>%group_by(testID, Target.Name, Quantity), t.st, relationship = 'many-to-one', by=c('testID', 'Target.Name', 'Quantity'))%>%
    	ungroup()%>%
      mutate(Amplified = paste0(n,'/1'))%>%
      mutate_if(is.numeric, round, digits = 3)%>%
      arrange('EventID')%>%
      dplyr::select(-n)
      st2

    # make table 

      allampRatio <- bind_rows(fs2, c2, st2) %>%
      	dplyr::select(Sample.Name, Target.Name, Assay.Role, eventID, method.type, testID, copies.sampavg, Amplified)%>%
      	arrange(testID)%>%
      	flextable()%>%
      	autofit()%>%
      	theme_zebra()
      allampRatio # for supplmentary


        save_as_image(allampRatio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_all.png', webshot = 'webshot2')     	
       	save_as_docx(allampRatio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_all.docx', webshot = 'webshot2')   

      onlyamp_ampratio <- bind_rows(fs2, c2) %>%
      	filter(Cq != 'NA')%>%
      	dplyr::select(Sample.Name, Target.Name, Assay.Role, eventID, method.type, testID, copies.sampavg, Amplified)%>%
      	arrange(testID)%>%      	
       	flextable()%>%
       	set_header_labels('Sample.Name' = 'Replicate','Target.Name' = 'Target', 'method.type' = 'Method', Assay.Role = 'Assay Role', eventID = 'Sampling \n\ Event ID', testID = 'qPCR \n\ Assay ID', copies.sampavg = 'Replicate \n\ Avg. Copy No.')%>%
      	autofit()%>%
      	theme_zebra()%>% 
      	align(align = 'center', part ='all')

       	save_as_image(onlyamp_ampratio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_onlyAmplifiedones.png', webshot = 'webshot2')     	
       	save_as_docx(onlyamp_ampratio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_onlyAmplifiedones.docx')


       	onlylabcontrols_check.table <- st2  %>%
      	dplyr::select(Sample.Name, Target.Name, Assay.Role, eventID, method.type, testID, copies.sampavg, Amplified)%>%
      	arrange(testID)%>%
      	flextable()%>%
      	autofit()%>%
      	theme_zebra()
      	onlylabcontrols_check.table


 		save_as_image(allampRatio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_all.png', webshot = 'webshot2')     	
   	save_as_docx(allampRatio, path ='EDNA/data_edna/figures_and_tables/qPCR_amplificationresults_bySamplingReplicate_all.docx', webshot = 'webshot2')  


#########
## - Comparing methods - reliable amplifications
#########

      ## How many metaprobes/waterbottles had 3/3 of amplifications for field samples

		a <- fs2 %>%
		filter(grepl("\\.1$|\\.2$|\\.3$", Sample.Name)) %>%
        filter(reliable =='1') 
      	aa <- a %>% group_by(Sample.Name, method.type)%>%
      	count()%>%
      	group_by(method.type)%>%
      	count()
      	aa # 10 metaprobes field replicates amplified; 31 waterbottle field replicates amplified, irrespective of species, only field sampes (not field controls) - tech reps, not samples
      	
      	aaa <- a %>% group_by(Target.Name,method.type)%>%
      	count()
    	#aaa.table<-aaa %>%      	
       	flextable()%>%
       	set_header_labels('Target.Name' = 'Target', 'method.type' = 'Method', n = 'No. of\n\ 3/3 Amp.' )%>%
      	autofit()%>%
      	theme_zebra()%>% 
      	align(align = 'center', part ='all')

      	save_as_image(aaa.table, path ='EDNA/data_edna/figures_and_tables/no_of_replicates_per_method_amplified3of3.png', webshot = 'webshot2')

     ## would be better as a hist/bar plot
      	a.plot <- ggplot(data = aaa, aes(y = n, x = Target.Name, fill = method.type))+
      		geom_col()+
      		theme_bw()+    
      		scale_fill_manual(labels = c('Metaprobe', 'Water Bottle'), values=c('#DAA507','#8EC7D2'))+    
      		labs(fill="Method", x = 'Species', y = 'Detections')+
      		theme(legend.position = 'bottom', legend.direction = 'horizontal') 
  	  	
  	  	ggsave(a.plot, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/no_of_replicates_per_method_amplified3of3_plot.png', device = 'png', units = 'in', height = 4.5, width = 4.5, dpi = 550)


#########
## - Comparing methods - concentration over time
#########
  
  d <- fieldsamples %>% distinct(Sample.Name, .keep_all = TRUE)%>%
  		mutate(date = ymd(sampling.date)) # THere is only one concentration value per Sampling replicate so we don't need every row, just one per. (n= 101)
  	head(d)

  dnacont_overtime<- ggplot(d, aes(date,dnaCont, group = method.type, pch = method.type, col = method.type, fill = method.type))+
   	geom_point()+
    geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
 	scale_color_manual(labels = c('Metaprobe', 'Water Bottle'), values=c('#DAA507','#8EC7D2'))+
    scale_fill_manual(labels = c('Metaprobe', 'Water Bottle'), values=c('#DAA507','#8EC7D2'))+
    scale_shape_manual(labels = c('Metaprobe', 'Water Bottle'), values = c(19,17))+
    xlab('Sampling Date')+ 
    ylab('Nanodrop Concentration')+
    theme_bw()+
    theme(plot.title = element_text(size=10, face='italic'), legend.position = 'bottom')+
    labs(fill="Method", col= "Method", pch = "Method")+
    scale_x_date(date_breaks = 'months', date_labels = '%b-%y')

    dnacont_overtime

    ggsave(dnacont_overtime, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/dnaCont_overtime_bymethod_ALlspecies.png', device = 'png', units = 'in', height = 5, width = 4.5, dpi = 550)

 

#########
## - Comparing methods - Copies by sampling event (stat_sum of the replicates) 
#########

    taxa1 <- fieldsamples %>%
    	mutate(taxa = as.character(case_when(Target.Name == 'Engraulis' ~ 'Fish', Target.Name == 'Scomber' ~ 'Fish',Target.Name == 'Prionace' ~ 'Shark', Target.Name == 'Alopias' ~ 'Shark')))%>%
  		mutate(date = ymd(sampling.date))%>%
  		as_tibble()%>%
  		dplyr::select(Sample.Name, Target.Name, Cq, copies,copies.techrepavg, copies.sampavg,eventID,taxa, date, method.type)
  	taxa1

 
   stat_sum_event_speciesbycolor <- ggplot(data = taxa1)+
      stat_summary(aes(x = eventID, y = log(copies.techrepavg+1),col = Target.Name, pch = method.type), fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = mean, position = position_jitter(width = 0.3, height = .3))+
	  scale_color_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus'), guide = 'none')+
      ylab('Copies (log)')+
      xlab('Sampling Event')+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 40, hjust=1), plot.title = element_text(size=10, face='italic'), legend.position = 'bottom',legend.box="vertical")+
      labs(color="Species", pch = 'Method')#+
      #scale_x_date(date_breaks = 'weeks', date_labels = '%d-%b-%y')
	  ggsave(stat_sum_event_speciesbycolor, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/stat_sum_copies_bydate_medianTECHRepcopies.png', device = 'png', units = 'in', height = 6, width = 6, dpi = 800)


      taxa2 <- taxa1 %>% 
      	mutate(row = row_number()) %>%
		arrange(Target.Name, method.type) %>% 
      	pivot_wider(names_from = method.type, names_sort = TRUE, values_from = copies.sampavg)%>%
      	arrange(row) %>%
  		select(-row) %>%
  		mutate_at(9:10, ~replace_na(., 0)) %>%
  		ungroup()%>%
		group_by(eventID, Target.Name, Sample.Name)%>%
  		slice(1)%>%
      	ungroup()%>%
      	dplyr::select(eventID, Target.Name, Sample.Name, metaprobe, waterbottle)
  		taxa2
  		##### save the taxa2 df because that's a handy one
  		write.csv(taxa2, 'EDNA/data_edna/copies_permethod_perEventID_perSpecies_foreachTechRep.csv')
  		d <- read.csv('EDNA/data_edna/copies_permethod_perEventID_perSpecies_foreachTechRep.csv')

	  segments <-  ggplot(d, aes(x = 'metaprobe', xend = 'waterbottle', y = log(metaprobe+1), yend = log(waterbottle+1), col = Target.Name, group = eventID))+
	      geom_segment(position = position_jitter(width = 0, 1), lwd = 0.25)+
	      scale_color_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus'))+
	      scale_x_discrete(labels = c('Metaprobe', 'Water Bottle'))+
	      theme_bw()+
	      theme(legend.position = 'bottom',legend.text = element_text(face='italic'))+
	      ylab(NULL)+
	      xlab(NULL)+
	      labs(col= "Species")+
	      guides(color = guide_legend(override.aes=list(size = 1.5)))
	    ggsave(segments, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/withinsamplingEvent_MEANcopies_byMethod_bySpecies.png', device = 'png', units = 'in', height = 5, width = 6, dpi = 800)


	both <- stat_sum_event_speciesbycolor + segments + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
	   ggsave(both, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/stat_sum_and_segments_speciesBYColor.png', device = 'png', units = 'in', height = 5, width = 11.5, dpi = 800)


#########
## - Comparing methods - Copies over time: all, fish grouped, sharks grouped
#########

    taxa <- fieldsamples %>%
    	mutate(taxa = as.character(case_when(Target.Name == 'Engraulis' ~ 'Fish', Target.Name == 'Scomber' ~ 'Fish',Target.Name == 'Prionace' ~ 'Shark', Target.Name == 'Alopias' ~ 'Shark')))%>%
  		mutate(date = ymd(sampling.date))%>%
  		as_tibble()
  	taxa

	## all species together 
     copies_time_allspp <- ggplot(taxa, aes(date,log(copies.techrepavg+.1), group = method.type, pch = method.type, col = method.type, fill = method.type))+
      geom_point()+
      geom_smooth(span = 2, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(labels = c('Metaprobe', 'Water Bottle'),values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(labels = c('Metaprobe', 'Water Bottle'),values=c('#DAA507','#8EC7D2'))+
      scale_shape_manual(labels = c('Metaprobe', 'Water Bottle'), values = c(19,17))+
      xlab('Sampling Date')+ 
      ylab('Copies (log)')+
      ggtitle('All Species')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'),legend.position = 'bottom')+
      labs(fill="Method", col= "Method", pch = 'Method')+
    scale_x_date(date_breaks = 'months', date_labels = '%b-%y')
	
	  ggsave(copies_time_allspp, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/copiesovertime_bymethod_ALlspecies_logTechRepCopies.png', device = 'png', units = 'in', height = 4, width = 4.5, dpi = 450)

	## fish and shark separately 

	  copies_time_fish_shark <- ggplot(taxa, aes(date,log(copies.techrepavg+0.1), pch = method.type,group = method.type, col = method.type, fill = method.type))+ geom_smooth(span = 2, se = TRUE, level = 0.95, alpha=0.3)+
	      scale_color_manual(labels = c('Metaprobe', 'Water Bottle'),values=c('#DAA507','#8EC7D2'))+
	      scale_fill_manual(labels = c('Metaprobe', 'Water Bottle'),values=c('#DAA507','#8EC7D2'))+
      		scale_shape_manual(labels = c('Metaprobe', 'Water Bottle'), values = c(19,17))+
	      geom_point()+
		  facet_wrap(~taxa, nrow = 1,scales = 'free')+
	      theme_bw()+
	      theme(strip.background = element_rect(fill = NA,colour = NA),plot.title = element_text(size=10, face='italic'),legend.position = 'bottom')+
		      xlab('Sampling Date')+ 
	      ylab('Copies (log)')+
	      labs(fill="Method", col= "Method", pch = 'Method')+
	    scale_x_date(date_breaks = 'months', date_labels = '%b-%y')

	## all with faceted taxa plots combined

		   tt <- taxa %>%mutate(taxa = "All")
		   ttt <- bind_rows(taxa,tt)

 		 copies_by_taxa <- ggplot(ttt, aes(date,log(copies.techrepavg+0.1), group = method.type, pch = method.type,col = method.type, fill = method.type))+ geom_smooth(span = 2, se = TRUE, level = 0.95, alpha=0.3)+
		      scale_color_manual(labels = c('Metaprobe', 'Water Bottle'),values=c('#DAA507','#8EC7D2'))+
		      scale_fill_manual(labels = c('Metaprobe', 'Water Bottle'),values=c('#DAA507','#8EC7D2'))+
      		  scale_shape_manual(labels = c('Metaprobe', 'Water Bottle'), values = c(19,17))+
		      geom_point()+
		      facet_wrap(~taxa, nrow = 1,scales = 'free')+
		      theme_bw()+
		      theme(strip.background = element_rect(fill = NA,colour = NA),plot.title = element_text(size=10, face='italic'),legend.position = 'bottom')+
 		      xlab('Sampling Date')+ 
		      ylab('Copies (log)')+
		      labs(fill="Method", col= "Method", pch = 'Method')+
		    scale_x_date(date_breaks = 'months', date_labels = '%b-%y')
		   #copies_by_taxa


	## stack species ones 
		all_spp_solo <- ggplot(taxa, aes(date,log(copies.techrepavg+0.1), group = method.type, pch = method.type,col = method.type, fill = method.type))+ geom_smooth(span = 2, se = TRUE, level = 0.95, alpha=0.3)+
		      scale_color_manual(labels = c('Metaprobe', 'Water Bottle'),values=c('#DAA507','#8EC7D2'))+
		      scale_fill_manual(labels = c('Metaprobe', 'Water Bottle'),values=c('#DAA507','#8EC7D2'))+
      		  scale_shape_manual(labels = c('Metaprobe', 'Water Bottle'), values = c(19,17))+
		      geom_point()+
		      facet_wrap(~Target.Name, nrow = 1,scales = 'free',labeller = labeller(Target.Name = c(Engraulis = 'Engraulis encrasicolus', 'Scomber' = 'Scomber scombrus', 'Prionace' = 'Prionace glauca', 'Alopias' = 'Alopias vulpinas')))+
		      theme_bw()+
		      theme(strip.text = element_text(face = "italic"),strip.background = element_rect(fill = NA,colour = NA),plot.title = element_text(size=10, face='italic'),legend.position = 'bottom')+
 		      xlab('Sampling Date')+ 
		      ylab('Copies (log)')+
		      labs(fill="Method", col= "Method", pch = 'Method')+
		    scale_x_date(date_breaks = 'months', date_labels = '%b-%y')
		   # all_spp_solo


   #######
   ## - save
	  # for presentations/talks, individual with axes and guides
	  ggsave(copies_time_allspp, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/copiesovertime_bymethod_ALlspecies_logTechRepCopies.png', device = 'png', units = 'in', height = 4, width = 4.5, dpi = 450)
	  ggsave(copies_time_fish_shark, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/copiesovertime_bymethod_FISH_andSHARK_species_logTechRepCopies.png', device = 'png', units = 'in', height = 4, width = 4.5, dpi = 450)

	  # for pub, faceted and with axes and guides collected 
	  ggsave(copies_by_taxa, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/copiesovertime_bymethod_byTaxa_logTechRepCopies.png', device = 'png', units = 'in', height = 4, width = 8.5, dpi = 450)
	  ggsave(all_spp_solo, file = 'EDNA/data_edna/figures_and_tables/comparingmethods/copiesovertime_bymethod_SpeciesIndividually_logTechRepCopies.png', device = 'png', units = 'in', height = 4, width = 8.5, dpi = 450)



#########
## -  Glmer (stat. tests): effect of method on copies over time for all species and tazon groups
#########

# normal with log10 data 
# presence/absence species with binomial fam. 
# hurdle models, binomial yes/no detection, then if yes how much (gamma binomial or beta binomial)

# change NAs to zeros - might do back at pipe. do models in brms. 
	  # when all 3 tech reps are na should be zero for sampling replicate 
	  # think about the loq-lod threshold and reassess/re-work it. overly conservative 
# occupancy model 

	# select down the df 
		
	data <- fieldsamples %>%
		dplyr::select(Sample.Name, Target.Name, eventID, copies, copies.techrepavg,method.type)%>%
		as_tibble()%>%
    	mutate(taxa = as.character(case_when(Target.Name == 'Engraulis' ~ 'Fish', Target.Name == 'Scomber' ~ 'Fish',Target.Name == 'Prionace' ~ 'Shark', Target.Name == 'Alopias' ~ 'Shark')))
	data
	stopifnot(nrow(data) == 1119) # check 

	# random effect on eventID

	# irregardless of target species
	lm1 <- glmer(copies.techrepavg ~ method.type * taxa + (1|eventID), data = data, family = poisson(link = 'log'))
	lm1
	saveRDS(lm1, 'EDNA/data_edna/glmer1_methodtype_copiesTechRep_intxTaxa_REeventid_Poisson.RDS')

	summary <- modelsummary(lm1, coef_rename = c('method.typewaterbottle' = 'Water Bottle','taxaShark' = 'Sharks (detectability)','method.typewaterbottletaxaShark' = 'Water Bottle*Sharks', 'method.metaprobe' = 'Metaprobe'),fmt=3,estimate='estimate', statistic='conf.int',stars=TRUE,conf_level=0.95,output='flextable', gof_omit = 'BIC|ICC|RMSE')%>%
		theme_zebra()%>%
		set_header_labels('Model 1' = 'Method & Target\n\ Taxa GLMM')%>%
	    align(align = 'center', part = 'all')%>%
	    font(fontname = 'Arial', part = 'all')%>%
	    fontsize(size = 10, part = 'all')%>%
	    autofit()
	summary # this table pairs with the geom_smooth plots 

	save_as_image(summary, 'EDNA/data_edna/figures_and_tables/glmer/summary_glmer_eventIDrandom_copies_byMethod_byTaxa.png', webshot = 'webshot2')
	save_as_docx(summary, path = 'EDNA/data_edna/figures_and_tables/glmer/summary_glmer_eventIDrandom_copies_byMethod_byTaxa.docx')

  
	## t-test between techrep copies and sample rep copies, grouped by method type. 
		## Testing for a significant variation between sampling replicates (the tech rep copies umber is the average of tech reps per smapling replicate, the sample rep copies is the sampling replicate average).
		## insignificant difference  = high similarity of replicates to the mean = agreement between replicates
		## significant difference = low similarity of replicates to the mean = higher variability between replicates of a sampling replicate

		data <- read.csv('EDNA/data_edna/qPCRresults/processedQPCRresults_cornwallspecies_june2024.csv')%>%
  			mutate(Sample.Name = as.factor(Sample.Name))%>%
  			as_tibble()
		data

		# grouping level = method.type 
		cor.test.methods<- data %>% 
		  nest(data = -method.type) %>%
		  mutate(cor=map(data,~cor.test(.$copies, .$copies.techrepavg, method = "kendall"))) %>%
		  mutate(tidied = map(cor, tidy)) %>% 
		  unnest(tidied) %>% 
		  select(-data, -cor)

		 cor.test.methods.table <- cor.test.methods %>%
		 	dplyr::select(-method, -alternative)%>%
		 	mutate(p.value = case_when(p.value<0.001 ~ 'p<0.001', p.value>=0.001 ~ as.character(p.value)))%>%
		 	flextable() %>%
		    set_header_labels(method.type = 'Method', 'estimate' = 'Estimate', statistic = 'Statistic',p.value = 'p-value')%>%
		    theme_zebra()%>%
		    align(align = 'center', part = 'all')%>%
		    font(fontname = 'Arial', part = 'all')%>%
		    fontsize(size = 10, part = 'all')%>%
		    autofit()

		    save_as_image(cor.test.methods.table, 'EDNA/data_edna/figures_and_tables/cortest_methods_replicates_vs_sampRepMeans.png', webshot = 'webshot2')
  			save_as_docx(cor.test.methods.table, path = 'EDNA/data_edna/figures_and_tables/cortest_methods_replicates_vs_sampRepMeans.docx')

#########
## -  Effect of Soak Time 
#########

  soaks <- fieldsamples %>%
  	as_tibble%>%
    filter(method.type=='metaprobe')%>%
		dplyr::select(Sample.Name, Target.Name, eventID, copies, copies.techrepavg,method.type, time.timeIN, timeOUT)%>%
		mutate(timeOUT = as.character(paste0(timeOUT,':00')), time.timeIN = as.character(paste0(time.timeIN,':00')))%>%
  	mutate(timeOUT2 = hms::as_hms(strptime(timeOUT, "%H:%M:%S")), timeIN2 = hms::as_hms(strptime(time.timeIN, "%H:%M:%S")), soaktime = as.duration(timeOUT2-timeIN2))%>%
    mutate(soaktime.min = as.numeric(soaktime)/60)%>%
    mutate(soaktime.hr = as.numeric(soaktime.min)/60)

 soaks.flex <-  soaks %>% 
 	dplyr::select(Sample.Name,eventID, Target.Name,copies.techrepavg,soaktime.min)%>%
 	group_by(Sample.Name, Target.Name)%>%
 	distinct(Sample.Name, Target.Name, .keep_all = TRUE)%>%
 	flextable() %>%
    set_header_labels(Sample.Name = 'Sample ID', 'eventID' = 'Event ID', Target.Name = 'Target Species',soaktime.min = 'Soak Time\n\ (min)', copies.techrepavg = 'Copies per Sampling\n\ Replicate (mean)')%>%
    theme_zebra()%>%
    align(align = 'center', part = 'all')%>%
    font(fontname = 'Arial', part = 'all')%>%
    fontsize(size = 10, part = 'all')%>%
    autofit()
    
  save_as_image(soaks.flex, 'EDNA/data_edna/qPCRresults/figures_and_tables/soaktimes_metaprobes_seconds_2023.png', webshot = 'webshot2')
  save_as_docx(soaks.flex, path = 'EDNA/data_edna/qPCRresults/figures_and_tables/soaktimes_metaprobes_seconds_2023.docx')


  ## glmer 

	lm2 <- glmer(copies.techrepavg ~ soaktime.hr + (1|eventID), data = soaks, family = poisson(link = 'log'))
	lm2
	saveRDS(lm2, 'EDNA/data_edna/glmer2_soaktime_copiesTechRep_REeventid.RDS')
	
	summary2 <- modelsummary(lm2, coef_rename = c('soak.time.hr' = 'Soak Time (hr)'),fmt=3,estimate='estimate', statistic='conf.int',stars=TRUE,conf_level=0.95,output='flextable', gof_omit = 'BIC|ICC|RMSE')%>%
		theme_zebra()%>%
		set_header_labels('Model 1' = 'Soak time GLMM')%>%
	    align(align = 'center', part = 'all')%>%
	    font(fontname = 'Arial', part = 'all')%>%
	    fontsize(size = 10, part = 'all')%>%
	    autofit()
	summary2 # this table pairs with the geom_smooth plots 

	save_as_image(summary2, 'EDNA/data_edna/figures_and_tables/glmer/summary_glmer_metaprobe_soaktime_eventIDrandomeffect.png', webshot = 'webshot2')
	save_as_docx(summary2, path = 'EDNA/data_edna/figures_and_tables/glmer/summary_glmer_metaprobe_soaktime_eventIDrandomeffect.docx')

	# plot - geom_smooth 

 	soaktime.. <- ggplot(soaks, aes(soaktime.hr,log(copies.techrepavg+1)))+
 	  geom_point(pch = 19, col = 'goldenrod4')+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.15, fill = '#DAA507', col = '#DAA507')+
      xlab('Soak time (min)')+ 
      ylab('Copies (log + 1)')+
      theme_bw()+
      theme(legend.position = 'bottom')
    # soaktime..

	 ggsave(soaktime.., file = 'EDNA/data_edna/figures_and_tables/glmer/soaktime_copies_techrepavg_geomsmooth.png', device = 'png', units = 'in', height = 4, width = 4, dpi = 850)


#########
## -  Model Validation - glmers for method*taxa and soaktime
#########

	ypred1 <- as_tibble(predict(lm1))%>%rename(fitted = value)
	res1 <- as_tibble(residuals(lm1, type = 'pearson')) %>% rename(resid = value)
	diag1 <- bind_cols(ypred1, res1)
	n1 <- nrow(data)
	p <- length(fixef(lm1))+1
	r1 <-resid(lm1, type = 'pearson')
	Overdispersion <- sum(r1^2)/(n1-p)
	

	lm1.resid.fitted <- ggplot(data = diag1, aes(x = fitted, y = resid))+
		geom_point()+ 
		labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
		theme_bw()
	lm1.hist.resid <- ggplot(data = res1, aes(x = resid))+ 
		geom_histogram(binwidth = 25000, fill = 'black')+ 
		theme_bw()+
		labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
	qq1 <- ggplot(data=ypred1, aes(sample = fitted))+
	stat_qq(size=1,pch=21)+
	labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
	stat_qq_line(linetype=2, col='red')+
	theme_bw()

	diagnostics1 <- lm1.resid.fitted+lm1.hist.resid+qq1


	ypred2 <- as_tibble(predict(lm2))%>%rename(fitted = value)
	res2 <- as_tibble(residuals(lm2, type = 'pearson')) %>% rename(resid = value)
	diag2 <- bind_cols(ypred2, res2)
	
	lm2.resid.fitted <- ggplot(data = diag2, aes(x = fitted, y = resid))+
		geom_point()+ 
		labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
		theme_bw()
	lm2.hist.resid <- ggplot(data = res2, aes(x = resid))+ 
		geom_histogram(binwidth = 15,fill = 'black')+ 
		theme_bw()+
		labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
	qq2 <- ggplot(data=ypred2, aes(sample = fitted))+
	stat_qq(size=1,pch=21)+
	labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
	stat_qq_line(linetype=2, col='red')+
	theme_bw()

	diagnostics2 <- lm2.resid.fitted+lm2.hist.resid+qq2

	######
	## - save

	ggsave(diagnostics1, file = 'EDNA/data_edna/figures_and_tables/glmer/diagnostics_glmer1_methodtype_copiestechrepavg.png', device = 'png', units = 'in', height = 3.5, width = 8.5, dpi = 450)
	ggsave(diagnostics2, file = 'EDNA/data_edna/figures_and_tables/glmer/diagnostics_glmer2_soaktime_copiestechrepavg.png', device = 'png', units = 'in', height = 3.5, width = 8.5, dpi = 450)
 

#########
## -  Detections of species at locations
#########

# start simple, one point per amplification (sampling replicates)

	p <- fieldsamples.sf %>% 
			dplyr::select(eventID, Sample.Name, Target.Name,copies.techrepavg, method.type, geometry)%>%
			filter(!is.na(copies.techrepavg))
	p

		spp.locations <-ggplot()+
				geom_sf(data=coastalwater,fill='cadetblue4', col=NA, lwd=0.5, alpha = 0.5)+
				geom_sf(data=kernios, lwd=0.5)+
				geom_sf(data=st_jitter(p, factor = 0.01), aes(col = Target.Name, fill = Target.Name, size = log(copies.techrepavg+1)), alpha =0.25)+
	      scale_color_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus'))+
	      scale_fill_manual(values = alpha(c('#477939', '#799ecb', '#85cb7c', '#003a78'),0.25), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus'))+
	      scale_size(name = 'Copies (log)', limits=c(1,25), breaks = c(5,10,15,25))+
	      theme_bw()+
	      theme(legend.position = 'right', legend.direction = 'vertical', legend.title = element_text(size=8), legend.text=element_text(size=8))+
	      labs(fill="Species", col= "Species")
		ggsave(spp.locations, file = 'EDNA/data_edna/figures_and_tables/species_detections_locations_copies2size_jittered.png', device = 'png', units = 'in', height = 5.5, width = 6, dpi = 800)
  

		north <- st_crop(p, y=c(ymax =50.68334, ymin = 50.5, xmin = -4.75, xmax = -5.2))
		water.north <- st_crop(coastalwater, y=c(ymax =50.68334, ymin = 50.5, xmin = -4.75, xmax = -5.2))
		land.north <- st_crop(kernios, y=c(ymax =50.68334, ymin = 50.5, xmin = -4.75, xmax = -5.2))
		south <- st_crop(p, y=c(ymax =50.25, ymin = 49.46, xmin = -6.15, xmax = -4.85))
		water.south <- st_crop(coastalwater, y=c(ymax =50.25, ymin = 49.46, xmin = -6.15, xmax = -4.85))
		land.south <- st_crop(kernios, y=c(ymax =50.25, ymin = 49.46, xmin = -6.15, xmax = -4.85))
	
		magnified.south <-ggplot()+
				geom_sf(data=water.south,fill='cadetblue4', col=NA, lwd=0.5, alpha = 0.5)+
				geom_sf(data=land.south, lwd=0.5)+
				geom_sf(data=st_jitter(south, factor = 0.01), aes(col = Target.Name, fill = Target.Name, size = log(copies.techrepavg)), alpha=0.25)+
	      scale_color_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus', guide = 'none'))+
	      scale_fill_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus', guide = 'none'))+
				guides(fill = 'none', color = 'none', size = 'none')+
	      theme_bw()+
	      scale_x_continuous(breaks = seq(-6, -4.8, by = 0.4))
		ggsave(magnified.south, file = 'EDNA/data_edna/figures_and_tables/species_detections_locations_copies2size_jittered_magnified_south.png', device = 'png', units = 'in', height = 5.5, width = 6, dpi = 800)
  
		magnified.north <-ggplot()+
				geom_sf(data=water.north,fill='cadetblue4', col=NA, lwd=0.5, alpha = 0.5)+
				geom_sf(data=land.north, lwd=0.5)+
				geom_sf(data=st_jitter(north, factor = 0.01), aes(col = Target.Name, fill = Target.Name, size = log(copies.techrepavg)), alpha=0.25)+
	      scale_color_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus', guide = 'none'))+
	      scale_fill_manual(values = c('#477939', '#799ecb', '#85cb7c', '#003a78'), labels = c('A. vulpinas', 'E. encrasicolus','P. glauca','S.scombrus', guide = 'none'))+
				guides(fill = 'none', color = 'none', size = 'none')+
	      theme_bw()+
	      scale_x_continuous(breaks = seq(-5.1, -4.75, by = 0.15))
		ggsave(magnified.north, file = 'EDNA/data_edna/figures_and_tables/species_detections_locations_copies2size_jittered_magnified_north.png', device = 'png', units = 'in', height = 5.5, width = 6, dpi = 800)
  

		all <- spp.locations + magnified.north + magnified.south + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.justification = 'center')
		ggsave(all, file = 'EDNA/data_edna/figures_and_tables/species_detections_locations_copies2size_jittered_setof3.png', device = 'png', units = 'in', height = 4, width = 8.5, dpi = 850)


    taxa1 <- fieldsamples %>%
    	mutate(taxa = as.character(case_when(Target.Name == 'Engraulis' ~ 'Fish', Target.Name == 'Scomber' ~ 'Fish',Target.Name == 'Prionace' ~ 'Shark', Target.Name == 'Alopias' ~ 'Shark')))%>%
  		mutate(date = ymd(sampling.date))%>%
  		as_tibble()%>%
  		filter(copies.sampavg != 'NA')%>%
  		dplyr::select(Sample.Name, Target.Name, Cq, copies,copies.techrepavg, copies.sampavg,eventID,taxa, date, method.type)
  	taxa1


