## - Data cleaning, eDNA in Cornwall 

## created 29 September 2023. updated april 2024 for results tables from the Quant Flex Studio 7 in the ESI. 
## by Molly M Kressler 

#################################
## - eDNA sampling locations 
#################################
##  REPEAT WITH NEW DATA DOWNLOADS
#################################
## created by Molly Kressler :: September 2023
#################################

setwd('/Users/mollykressler/Documents/EDNA/data_edna/')
pacman::p_load(sf,dplyr,lubridate,readr,readxl,lubridate,ggplot2,patchwork,cowplot)

## Load at Start
	 ext <- read.csv('data_downloads_notpublic/extractiondata_edna_endofyear2023.csv')%>%
	 	dplyr::select(sampleID,eppendorfID,dnaCont,ratio260.280,ratio260.230,date)
	 samp <- read.csv('data_downloads_notpublic/samplingdata_edna_endofyear2023.csv')%>%rename(Latitude=longitude,Longitude=latitude)%>%
	 	mutate(Latitude = as.character(Latitude),Longitude = as.character(Longitude))
	 head(samp)

## Calculate number of passengers reached - outreach metrics
	 outreach <- samp %>% 
	 	group_by(date) %>%
	 	tally(passengers) %>%
	 	tally(n)
	 outreach # number of people reached

	 vesselcount <- samp %>% 
	 	group_by(vessel) %>%
	 	tally()
	 vesselcount # number of people reached
	
	 sampligntripscounted <- samp %>% 
	 	group_by(date) %>%
	 	slice_head(n=1) %>%
	 	tally() %>%
	 	tally(n) 
	 sampligntripscounted # number of people reached


## Clean and add location and eventID to each eppendorf row in extraction.csv
	 e2 <- ext %>%
	 	rename(extraction.date=date)%>%
	 	mutate(sampledat = str_split(sampleID,'-',simplify=TRUE)[,1],
	 		type = str_split(sampleID,'-',simplify=TRUE)[,2],
	 		AorB = str_split(sampleID,'-',simplify=TRUE)[,3], 
	 		method.subsample = str_split(sampleID,'-',simplify=TRUE)[,4])%>%
	 	mutate(eventID = paste0(sampledat,"-",type,"-",AorB), .before=sampleID)%>%
	 	mutate(sampling.date=dmy(sampledat),month=month(sampling.date))%>%
	 	mutate(across(contains('.quantity'),~replace(.,is.na(.),0)))%>%
	 	mutate(across(contains('.presence'),~replace(.,is.na(.),0)))%>%
	 	filter(!is.na(dnaCont)) # empty last row
	 	head(e2)
	 	sapply(e2,class)

## 	Left_join locations from samp to e2 using the event ID

	d1 <- left_join(e2, samp%>%
		dplyr::select(eventID,Latitude,Longitude,recorded.by,vessel,passengers,time.timeIN,timeOUT), by = 'eventID',relationship='many-to-one',multiple='first')%>%
		dplyr::select(-blueshark.presence,-blueshark.Cq,-blueshark.copies,-engraulis.presence, -engraulis.Cq, -engraulis.copies)%>%
		mutate(timeOUT = ifelse(is.na(timeOUT),time.timeIN,timeOUT))%>%
	 	mutate(method.type = case_when(startsWith(as.character(method.subsample),'W') ~ 'waterbottle', TRUE ~ 'metaprobe'))	 	## mutate for timeOUT fills in the timeout column for water bottles to match the time.timeIN value. 
	head(d1)
	summary(d1)



## Save as CSV and SHP
	
	write.csv(d1,'eDNA_data_meta_and_qPCR_2023.csv') # extraction + meta
	
	sf1 <- st_as_sf(d1,coords=c('Longitude','Latitude'),crs='WGS84')%>%
		dplyr::select(-blueshark.presence,-blueshark.Cq,-blueshark.copies,-engraulis.presence, -engraulis.Cq, -engraulis.copies) #extraction + meta
		unique(sf1$geometry) # will be FEWER than in samp/sampsf1, but equal to how many locations you've done extractions and qPCR for 

	sampsf1 <- 	st_as_sf(samp,coords=c('Longitude','Latitude'),crs='WGS84')%>%
	 	mutate(type = str_split(eventID,'-',simplify=TRUE)[,2],
	 		sampledat = str_split(eventID,'-',simplify=TRUE)[,1],.before='geometry')%>%
	 	mutate(sampling.date=dmy(sampledat), month = as.character(month(sampling.date)),.before='geometry')%>%
	 	rename()%>%
	 	mutate(method.type = case_when(startsWith(as.character(sampleID),'W') ~ 'waterbottle', TRUE ~ 'metaprobe'),.before='geometry')		# sampling data
	 	sampsf1
		unique(sampsf1$geometry)

	st_write(sampsf1, 'eDNA_data_meta_and_qPCR_2023.shp',driver='ESRI Shapefile')

##############
## Dataset for Undergaduate project 2023/34
	## sample ID, metadata, concentration and engraulis copies 

	data <- read.csv('compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv')

	write.csv(data, 'teachingdataframes/eDNACornishBlue_2023undergraduateproject_concentrations_metadata_Engraulis_encrasicolus_qPCRresults.csv') # same exact data, but need to save twice so it's distinct. Later, if re-making this, you'll need to remove columsn for other species 

## Sub-sampled dataset for MSc project 2023/24
	## for them to come to grips with the data format 
	## sample ID, concentrations, metadata, and engraulis copies for 20 random samples	

	subsample1 <- data %>% slice_sample(n=20)

	write.csv(subsample1, 'teachingdataframes/eDNACornishBlue_2023MScproject_SUBSAMPLE_Engraulis_encrasicolus_qPCRresults.csv') # 


#################################
## - Tidy qPCR Results
#################################
## files are in weird formats, so need a workflow for dissecting out the data frame

	 ext2 <- read.csv('data_downloads_notpublic/extractiondata_edna_endofyear2023.csv')%>%
	 	dplyr::select(sampleID,eppendorfID,dnaCont,ratio260.280,ratio260.230,date)%>%
	 	rename(Sample.Name = eppendorfID) # use this to match smapleID to eppendorf IDs

######### SAMPLES & STANDARDS, repeated for all species 

	# NOVEMBER 2023 Engraulis encrasicolus 
		ent1 <- read_excel('qPCRresults/2023Engraulisencrasicolus/qPCRmax_output_excelfiles/En.encras_TEST1_NOV23_2023.xlsx', sheet = 'Results',range = 'A11:R59', .name_repair = 'universal', col_types='text')%>%
			mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))
		
		ent2 <- read_excel('qPCRresults/2023Engraulisencrasicolus/qPCRmax_output_excelfiles/En.encras_TEST2_NOV23_2023.xlsx', sheet = 'Results',range = 'A11:R59', .name_repair = 'universal', col_types='text')%>% 
			mutate(Sample.Name = case_when(Assay.Role != 'Unknown' ~ paste0(Assay.Role,'_',Quantity), Assay.Role == 'Unknown' ~ Sample.Name))%>%
			mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')
		
		ent3 <- read_excel('qPCRresults/2023Engraulisencrasicolus/qPCRmax_output_excelfiles/En.encras_TEST3_NOV23_2023.xlsx', sheet = 'Results',range = 'A11:R59', .name_repair = 'universal', col_types='text')
		
		ent4 <-read_excel('qPCRresults/2023Engraulisencrasicolus/qPCRmax_output_excelfiles/En.encras_TEST4_NOV23_2023.xlsx', sheet = 'Results',range = 'A11:R59', .name_repair = 'universal', col_types='text')%>%
			mutate(Sample.Name = case_when(Assay.Role != 'Unknown' ~ paste0(Assay.Role,'_',Quantity), Assay.Role == 'Unknown' ~ Sample.Name))%>%
			mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')
		
		ent5 <- read_excel('qPCRresults/2023Engraulisencrasicolus/qPCRmax_output_excelfiles/En.encras_TEST5_24112023.xlsx', sheet = 'Results',range = 'A11:R59', .name_repair = 'universal', col_types='text')%>%
			mutate(Sample.Name = case_when(Assay.Role != 'Unknown' ~ paste0(Assay.Role,'_',Quantity), Assay.Role == 'Unknown' ~ Sample.Name))
		
		ent6 <- read_excel('qPCRresults/2023Engraulisencrasicolus/qPCRmax_output_excelfiles/En.encras_TEST6_NOV24_2023.xlsx', sheet = 'Results',range = 'A11:R59', .name_repair = 'universal', col_types='text')%>%
			mutate(Sample.Name = case_when(Assay.Role != 'Unknown' ~ paste0(Assay.Role,'_',Quantity), Assay.Role == 'Unknown' ~ Sample.Name))

		ent7 <- read_excel('qPCRresults/2023Engraulisencrasicolus/qPCRmax_output_excelfiles/En.encras_TEST7_NOV24_2023.xlsx', sheet = 'Results',range = 'A11:R59', .name_repair = 'universal', col_types='text')%>%
			mutate(Sample.Name = case_when(Assay.Role != 'Unknown' ~ paste0(Assay.Role,'_',Quantity), Assay.Role == 'Unknown' ~ Sample.Name))%>%
			mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		all_engraulis <- bind_rows(ent1, ent2, ent3, ent4, ent5, ent6, ent7)%>%
			dplyr::select(-Exclude, -Plate.Control)%>%
			mutate(sp = 'engraulis')
		all_engraulis
		all <- left_join(all_engraulis, ext2, by='Sample.Name')%>%
			 mutate(sampleID = case_when(Assay.Role != 'Unknown' ~ paste0(Assay.Role,'_',Quantity), Assay.Role == 'Unknown' ~ sampleID))%>%
			 rename(replicateID=sampleID)

	# APRIL 2024 
		# Standard Curve tests 
		standards1 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		standards2 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		standards3 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')


		# Engraulis encrasicolus
		assay1 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		assay2 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		assay3 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		# Mixed Species 
		assay4 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		# Sprattus sprattus
		assay5 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		assay6 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		assay7 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')
		
		# Scomber scombrus
		assay8 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		assay9 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')

		assay10 <- read_excel('qPCRresults/2024assays/.xlsx', sheet = 'Results',range = 'A41:X138', .name_repair = 'universal', col_types='text')%>%
		rename(Sample.Name = 'Sample Name')%>%
		mutate(Sample.Name = case_when(Task != 'Unknown' ~ paste0(Task,'_',Quantity), Task == 'Unknown' ~ Sample.Name))%>%
		mutate_if(is.character,str_replace_all,pattern='M',replacement='MP')



	### save tidyed dfs 
		write.csv(ent1,'qPCRresults/2023Engraulisencrasicolus/individual_assays/tidyed_results_En.encras_TEST1_NOV23_2023.csv')
		write.csv(ent2,'qPCRresults/2023Engraulisencrasicolus/individual_assays/tidyed_results_En.encras_TEST2_NOV23_2023.csv')
		write.csv(ent3,'qPCRresults/2023Engraulisencrasicolus/individual_assays/tidyed_results_En.encras_TEST3_NOV23_2023.csv')
		write.csv(ent4,'qPCRresults/2023Engraulisencrasicolus/individual_assays/tidyed_results_En.encras_TEST4_NOV23_2023.csv')
		write.csv(ent5,'qPCRresults/2023Engraulisencrasicolus/individual_assays/tidyed_results_En.encras_TEST5_24112023.csv')
		write.csv(ent6,'qPCRresults/2023Engraulisencrasicolus/individual_assays/tidyed_results_En.encras_TEST6_NOV24_2023.csv')
		write.csv(ent7,'qPCRresults/2023Engraulisencrasicolus/individual_assays/tidyed_results_En.encras_TEST7_NOV24_2023.csv')

		write.csv(all,'qPCRresults/2023Engraulisencrasicolus/tidyed_results_En.encras_TESTS_1to7_NOV2023.csv')


######### STANDARDS 

	ensc <- read_excel('qPCRresults/2023Engraulisencrasicolus/qPCRmax_output_excelfiles/En.encras_STANDARDCURVETEST_NOV24_2023.xlsx', sheet = 'Results',range = 'A11:R59', .name_repair = 'universal', col_types='text')%>%
		dplyr::select(-Exclude, -Plate.Control)%>%
		filter(Assay.Name != 'NA')%>% # filter out empty wells 
		mutate(Sample.Name = paste0(Assay.Role,'_',Quantity)) # adjust sample.name to include standard and the concentration
		
		stopifnot(nrow(ensc)==25) # 24 wells of standards, 1 well of NTC water, 1 row per well 

	prio <- read_excel('data_edna/qPCRresults/2023Prinoaceglauca/Pglauca_STANDARDCURVETEST_JUL312023.xlsx', sheet = 'Results',range = 'A11:R59', .name_repair = 'universal', col_types='text')%>%
		dplyr::select(-Exclude, -Plate.Control)%>%
		filter(Assay.Name == 'P.glauca')%>% # filter out empty wells 
		mutate(Sample.Name = paste0(Assay.Role,'_',Quantity)) # adjust sample.name to include standard and the concentration
		
		stopifnot(nrow(prio)==19) 

	## save standards csv
		write.csv(prio,'data_edna/qPCRresults/2023Prinoaceglauca/Pglauca_STANDARDCURVETEST_JUL312023.csv')


































