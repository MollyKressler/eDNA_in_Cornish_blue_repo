## Pipeline to Calculate copy number for technical and sample replicates from seawater samples analysesed qith qPCR

## Molly M Kressler :: 28 May 2024
## method for calculating copies updated after meetign with Rich 05.07.2024

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

	## Engraulis
	eng1 <- read_excel('qPCRresults/ESI_assays_spring2024/05062024-KRESSLER-ENGRAULIS-ASSAY1.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'eng1')%>%
		left_join(., neb, by = c('Target.Name'='sp'))

	eng2 <- read_excel('qPCRresults/ESI_assays_spring2024/05062024-KRESSLER-ENGRAULIS-ASSAY2.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'eng2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	eng3 <- read_excel('qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-ENGRAULIS-ASSAY3.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'eng3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Scomber
	sco1 <- read_excel('qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-SCOMBER-ASSAY4.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'sco1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	sco2 <- read_excel('qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-SCOMBER-ASSAY5.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'sco2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	sco3 <- read_excel('qPCRresults/ESI_assays_spring2024/11062024-KRESSLER-SCOMBER-ASSAY6.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'sco3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Prionace
	prio1 <- read_excel('qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-PRIONACE-ASSAY7.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'prio1', Target.Name = 'Prionace')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	prio2 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY8.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'prio2', Target.Name = 'Prionace')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	prio3 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY9.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'prio3', Target.Name = 'Prionace')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Alopias
	alo1 <- read_excel('qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-ALOPIAS-ASSAY10.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'alo1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	alo2 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY11.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'alo2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	alo3 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY12.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'alo3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 


	## one combo plate, all four species - for chapter 4 remove the BRUV ones
	combo1 <-  read_excel('qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY17.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'combo1', Target.Name = case_when(Target.Name == 'Prionance' ~ 'Prionace', Target.Name != 'Prionance' ~ Target.Name), Assay.Role = case_when(Assay.Role == 'NTC_NA' ~ 'NTC', Assay.Role != 'NTC_NA' ~ Assay.Role))%>%
		left_join(., neb, by = c('Target.Name'='sp')) %>%
		filter(!str_detect(Sample.Name, 'B'))
	
	combo2 <-  read_excel('qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY18.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'combo2', Target.Name = case_when(Target.Name == 'Prionance' ~ 'Prionace', Target.Name != 'Prionance' ~ Target.Name))%>%
		left_join(., neb, by = c('Target.Name'='sp')) %>%
		filter(!str_detect(Sample.Name, 'B'))


	combo3 <-  read_excel('qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY19.xlsx', sheet = 'Results',range = 'A41:L107', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role !='STANDARD')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'combo3', Target.Name = case_when(Target.Name == 'Prionance' ~ 'Prionace', Target.Name != 'Prionance' ~ Target.Name))%>%
		left_join(., neb, by = c('Target.Name'='sp')) %>%
		filter(!str_detect(Sample.Name, 'B'))

	## combine into one df, all species, all assays

	data <- rbind(eng1,eng2,eng3,sco1,sco2,sco3,prio1, prio2, prio3,alo1,alo2,alo3,combo1,combo2,combo3)
	data

	table <- data %>% 
		mutate(Cq = round(as.numeric(Cq),3),Ct.Mean = round(as.numeric(Ct.Mean),3),Ct.SD = round(as.numeric(Ct.SD),3),Quantity = round(as.numeric(Quantity),3))%>%
		dplyr::select(-Omit, -Well, -Assay.Role, -species, -dateofassay, -slope, -intercept, -efficiency, -LOD, -LOQ, -Quantity, -r.squared)%>%
		flextable()%>%
		autofit()%>%
		theme_zebra()

		 save_as_docx(table, path = 'qPCRresults/fieldsamples_cornwalledna_amprawdata.docx') 
		 # takes a while so dont run everytime
	#########
	## save
		data2 <- data %>% filter(Assay.Role !='NTC')
		#write.csv(data2, 'qPCRresults/collatedqPCRresults_2024_cornwallspecies.csv')
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
		    mutate(copies = ifelse(is.na(Cq), 0, copies))%>% # this says if the original Cq is NA, set the copies to 0, but it will not directly affect the calculations in the next step which calculate the tech rep averages based on loq_check reliability. The results for copies mean is a slight reduction in the mean overall, but no differrence to the min or max. 
		    mutate(copies.techrepavg = if_else(reliable, median(copies), 0))
		    summary(output)
	# (C) - July 2024: I dont think this is logical anymore, calculatign the sampling event mean. bc we treat each sampling relicate as individual samples. it doesn't changes results because we only use the techrepavg in the analysis dcuments. But I'm leaving this as a note here for future reference. 
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

	## Import metadata 

	meta.sf <- st_as_sf(st_read(('eDNA_data_meta_and_qPCR_2023.shp'),crs='WGS84'))%>%
		rename(sampleID = samplID, eppendorfID = eppndID, ratio260.280 = r260_28, ratio260.230 = r260_23, sampledate = sampldt, method.subsample = mthd_sb, sampling.date = smplng_, recorded.by = rcrdd_b, passengers = pssngrs, time.timeIN = tm_tmIN, method.type = mthd_ty)

		## Don't re-run. Here for documentation. meta.sf is missing a data entry for eppendorf 13.4 which is the field control for 15082023-R-B. Construct it and row bind it to meta.sf. 

			m <- meta.sf %>% 
				filter(eppendorfID == '13.3') %>%
				mutate(sampleID = '15082023-R-B-WBT13.4', dnaCont = 0, ratio260.280 = 1.250, ratio260.230 = -0.068, method.subsample = 'WBT1.4', eppendorfID = '13.4')
			mmm <- rbind(m, meta.sf)

		## Don't re-run. Here for documentation. Meta.sf timeIN/timeOUT are incorrect based on original data on OneDrive. Updatign here, then will re-join to copies dataset. 

			raw <- read.csv('raw_sampling_metadata_2023.csv')%>%
				dplyr::select(eventID, sampleID, time.timeIN, timeOUT)%>%
					as_tibble()%>%
				mutate(timeOUT = case_when(is.na(timeOUT) ~ time.timeIN, !is.na(timeOUT) ~ timeOUT))%>%
				filter(str_detect(sampleID,'M'))
			raw

			meta2 <- meta.sf%>%
				left_join(raw, by='eventID', suffix = c('','_join'))%>%
				mutate(timeOUT = if_else(method.type == 'metaprobe', timeOUT_join, timeOUT))%>%
				select(-timeOUT_join)

			meta.sf <- meta2 %>%
				dplyr::select(-time.timeIN_join, -sampleID_join)		
			## overwrite the old file
			
			#st_write(meta.sf, 'eDNA_data_meta_and_qPCR_2023.shp', driver = 'ESRI Shapefile', append = FALSE)
			#st_write(meta.sf, 'eDNA_data_meta_and_qPCR_2023.csv', driver = 'CSV', append = FALSE)


	## join

	joined <- left_join(output2, meta.sf%>%mutate(eppendorfID = str_remove(eppendorfID, '[P]')), by = c('Sample.Name' = 'eppendorfID'))%>%
		st_as_sf 

	joined
	joined%>%filter(is.na(sampleID)) #should be 0 rows

#########
## - Step 4: save data 
#########

	st_write(joined, 'qPCRresults/processedQPCRresults_cornwallspecies_june2024.csv', driver = 'CSV', delete_dsn = TRUE, delete_layer = TRUE)
	st_write(joined, 'qPCRresults/processedQPCRresults_cornwallspecies_june2024.shp', driver = 'ESRI Shapefile', delete_dsn = TRUE, delete_layer = TRUE)

  s <- st_as_sf(st_read('qPCRresults/processedQPCRresults_cornwallspecies_june2024.shp'))%>%
		rename(Assay.Role = Assy_Rl, Target.Name = Trgt_Nm, sampleID = samplID, ratio260.280 = r260_28, ratio260.230 = r260_23, sampledate = sampldt, method.subsample = mthd_sb, sampling.date = smplng_, recorded.by = rcrdd_b, passengers = pssngrs, time.timeIN = tm_tmIN, method.type = mthd_ty, Reporter = Reportr, Quencher = Quenchr, Sample.Name = Smpl_Nm, intercept = intrcpt, r.squared = r_squrd, efficiency = effcncy, loq_check = lq_chck, reliable = reliabl, copies.techrepavg = cps_tch, copies.sampavg = cps_smp, Well.Position = Wll_Pst)












