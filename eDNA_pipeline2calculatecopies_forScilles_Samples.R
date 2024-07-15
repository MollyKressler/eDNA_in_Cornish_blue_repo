## Pipeline to Calculate copy number for technical and sample replicates from seawater samples analysesed qith qPCR

## Molly M Kressler :: 28 May 2024

#########
## - Background
#########

## specifically written for xlsx files as produced by Applied Biosciences Quant 7 Studio Flex qPCR machine, but could be adapted with modifications to initial data importing. 

## Caculates copies for technical replicates with conditions regarding the limit of detection (LOD) and the limit of quantification (LOQ). The LOD and LOQ are identified in the Standard Curve Tests performed for each species and are stored with the standard curve equation data in a separate data file. 

## Definitions
	## LOD: limit of detection, lowest concentration (or, practically, the corresponding Cq) of standard template DNA with at least one amplified technical replicate.  Species specific, identified by standard curve tests with qPCR assay
	## LOQ: limit of quantification, lowest concentration (or, practcally, correspodning Cq) of standard template DNA with 3 positive technical replicates. Species specific, identified by standard curve tests with qPCR assay.
	## technical replicate, or tech rep: the replicates within a qPCR assay for each sampling replicate, n=3.
	## sample replicate, or samp rep: the replicates within a single sampling event, n=3 for metaprobe and water bottle sampling


#########
## - Broad steps
#########

## (1) Import and clean datasets from qPCR assay output files, filtering out Laboratory NTCs and Standards/Positives. Done for each species, saving individual species files, before collating species and saving as one large file. 
## (2) Calculate copies per technical replicate, and per sampling replicate, for field samples and field controls. Then, calculate sampling event average copy number excluding the field controls, but preserving these

#########
## - Step 0: Define and Load workspace
#########

pacman::p_load(sf,dplyr,tidyverse,lubridate,readr,readxl,stringr, flextable)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/EDNA/data_edna/')

#########
## - Step 1: tidying exports from qPCR machine
#########
	## Standard Curves, all species
	neb <- read_excel('qPCRresults/standard_curve_equations_qPCR_2023species.xlsx', sheet = 'Sheet1', range = 'A1:I10',.name_repair = 'universal', col_types = c('text','text','numeric','numeric','numeric','numeric','date','numeric','numeric'))%>%
		mutate(dateofassay=as_date(dateofassay))%>%
		filter(dateofassay >= '2024-06-28') 
	neb

	## Assays
	assay1 <- read_excel('scillies_data_edna/qPCR_results_scillies/08072024-KRESSLER-SCILLES-ASSAY1.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'assay1')%>%
		left_join(., neb, by = c('Target.Name'='sp'))
	assay2 <- read_excel('scillies_data_edna/qPCR_results_scillies/08072024-KRESSLER-SCILLES-ASSAY2.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'assay2')%>%
		left_join(., neb, by = c('Target.Name'='sp'))
	assay3 <- read_excel('scillies_data_edna/qPCR_results_scillies/10072024-KRESSLER-SCILLES-ASSAY3.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'assay3')%>%
		left_join(., neb, by = c('Target.Name'='sp'))

	## combine into one df, all species, all assays

	data <- rbind(assay1, assay2, assay3) 
	data

	#########
	## save
	
		write.csv(data, 'scillies_data_edna/qPCR_results_scillies/collatedqPCRresults_2024_IOSspecies.csv')
		data <- read.csv('scillies_data_edna/qPCR_results_scillies/collatedqPCRresults_2024_IOSspecies.csv')%>%
		dplyr::select(-X)

#########
## - Step 2: calculate copies
#########
## write function to calculate copies
 	 calculate_copies <- function(intercept, slope, Cq.adj) {
    	copies <- 10^((Cq.adj - intercept)/slope)
    	return(copies)
    	}

    ## apply function in 2 parts in accordance with LOQ and LOD criteria. (A) First to calculate tech rep average copies, and then (B) the Field Sample and Field Control average (average of tech reps), then (C) join dataframes.
    ## updated on 05.07.2024: NAs at the end need to be zeros (0).  

	#(A) Identify where the NTC amplified before the LOQ

		o <- data  %>%
		  group_by(testID, Target.Name) %>%
		  mutate(NTC.Cq = ifelse(any(Assay.Role == "NTC", na.rm = TRUE), Cq[Assay.Role == "NTC"], NA)) %>%
 		  mutate(NTC.Amplified = if_else(Assay.Role == "NTC" & !is.na(NTC.Cq) & NTC.Cq < LOQ, TRUE, FALSE)) %>%
 		  ungroup() %>%
 		  group_by(testID, Target.Name) %>%
		  mutate(NTC.Amplified = if_else(Assay.Role != "NTC" & any(Assay.Role == "NTC", na.rm = TRUE), 
		      first(NTC.Amplified[Assay.Role == "NTC"], default = FALSE),NTC.Amplified)) %>%
		  select(-NTC.Cq)%>%
		  filter(Assay.Role != 'NTC') # filter NTCs back out, all we want are the Samples for this file.
  		o

	# (B) Where technical replicates are reliable, this will adjust NA Cq values to the LOD which is arbitrarily low. Where reliable but non-NA the Cq stays the same. 
		output <- group_by(o, Target.Name, Sample.Name) %>%
		    mutate(., loq_check = ifelse(Cq <= LOQ & !is.na(Cq), 1, NA)) %>%
		    mutate(., loq_check = ifelse(any(loq_check == 1), 1, NA)) %>%
		    mutate(Cq.adj = as.numeric(ifelse(loq_check == 1 & is.na(Cq), LOD, Cq)))%>%
		    mutate(., reliable = ifelse(loq_check == is.numeric(loq_check), TRUE))%>%
		    mutate(reliable = replace_na(reliable,FALSE))%>%
		    mutate(Cq.adj = case_when(reliable == 'TRUE' ~ as.numeric(Cq), reliable == 'FALSE' ~ LOD))%>% 
		    mutate(copies = ifelse(reliable, calculate_copies(intercept, slope, Cq.adj), 0))%>%
		    mutate(copies = ifelse(is.na(Cq), 0, copies))%>% 
		    mutate(copies.techrepavg = if_else(reliable, median(copies), 0))
		    summary(output)
	# (C) 
		fieldsamp_averagecopies <- output %>%
		  filter(grepl("\\.1$|\\.2$|\\.3$", Sample.Name)) %>%
		  group_by(Target.Name, Sample.Name) %>%
		  summarise(fieldsamples.copies = median(copies.techrepavg, na.rm = TRUE))
		fieldsamp_averagecopies
		summary(fieldsamp_averagecopies)
	# (D) 
		output2 <- output %>%
		  left_join(fieldsamp_averagecopies, by = c("Target.Name", "Sample.Name")) %>%
		  mutate(copies.sampavg = if_else(grepl("\\.4$", Sample.Name), copies.techrepavg, fieldsamples.copies))%>%
		  dplyr::select(-fieldsamples.copies)
		output2%>%print(n=20)
		summary(output2)
	

#########
## - Step 3: join to metadata and spatial information
#########

	## Only run for first open. match eppendorf/Sample.Name to metadata, convert to sf object.
		meta <- read_excel('scillies_data_edna/Field_METADATA.xlsx', sheet = 'Sheet1',range = 'A1:N17', .name_repair = 'universal', col_types='text')%>%
			dplyr::select(eventID, timeIN, timeOUT, methodtype, Latitude, Longitude, Transect)
		key <- read_excel('scillies_data_edna/extractions_scilles.xlsx', sheet = 'extractions',range = 'A1:D53', .name_repair = 'universal', col_types='text')%>%
			mutate(extr.date = as.numeric(extr.date))
		mm <- left_join(key, meta,by = 'eventID', relationship = 'many-to-many') %>% 
				mutate(extr.date = as_date(extr.date), extr.date = update(extr.date, year = 2024))%>%
				filter(!str_detect(eppendorfID, 'extb'))
		meta.sf <- st_as_sf(mm,coords = c('Longitude', 'Latitude'), crs = 'WGS84') 

		######
		## - save
		######

			st_write(meta.sf,'scillies_data_edna/IOS_metaprobe_deployments_metadata.shp', driver = 'ESRI Shapefile')
			st_write(meta.sf,'scillies_data_edna/IOS_metaprobe_deployments_metadata.csv', driver = 'CSV')

	## Import metadata (run every open)

	meta.sf <- st_as_sf(st_read(('scillies_data_edna/IOS_metaprobe_deployments_metadata.shp'),crs='WGS84'))%>%
		rename(sampleID = smplnID, eppendorfID = eppndID, extr.date = extr_dt, transect = Transct, method.type = mthdtyp)
	## join

	joined <- left_join(output2, meta.sf%>%mutate(eppendorfID = str_remove(eppendorfID, '[P]')), by = c('Sample.Name' = 'eppendorfID'))%>%
		st_as_sf%>%
		mutate(sampleID = case_when(str_detect(Sample.Name, 'extbl') == TRUE ~ Sample.Name, str_detect(Sample.Name, 'extbl') == FALSE ~ sampleID))

	joined
	joined%>%filter(is.na(sampleID)) #should be 0 rows

#########
## - Step 4: save data 
#########

	st_write(joined, 'scillies_data_edna/qPCR_results_scillies/processedQPCRresults_IOS_june2024.csv', driver = 'CSV', delete_dsn = TRUE, delete_layer = TRUE)
	st_write(joined, 'scillies_data_edna/qPCR_results_scillies/processedQPCRresults_IOS_june2024.shp', driver = 'ESRI Shapefile', delete_dsn = TRUE, delete_layer = TRUE)

	## check it saved appropriately 
  	s <- st_as_sf(st_read('scillies_data_edna/qPCR_results_scillies/processedQPCRresults_IOS_june2024.shp'))%>%
		rename(Assay.Role = Assy_Rl, Target.Name = Trgt_Nm, sampleID = samplID,  method.type = mthd_ty, Reporter = Reportr, Quencher = Quenchr, Sample.Name = Smpl_Nm, intercept = intrcpt, r.squared = r_squrd, efficiency = effcncy, loq_check = lq_chck, reliable = reliabl, copies.techrepavg = cps_tch, copies.sampavg = cps_smp, Well.Position = Wll_Pst)












