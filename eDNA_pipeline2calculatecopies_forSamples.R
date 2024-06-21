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
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'eng1')%>%
		left_join(., neb, by = c('Target.Name'='sp'))

	eng2 <- read_excel('qPCRresults/ESI_assays_spring2024/05062024-KRESSLER-ENGRAULIS-ASSAY2.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'eng2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	eng3 <- read_excel('qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-ENGRAULIS-ASSAY3.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'eng3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Scomber
	sco1 <- read_excel('qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-SCOMBER-ASSAY4.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'sco1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	sco2 <- read_excel('qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-SCOMBER-ASSAY5.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'sco2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	sco3 <- read_excel('qPCRresults/ESI_assays_spring2024/11062024-KRESSLER-SCOMBER-ASSAY6.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'sco3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Prionace
	prio1 <- read_excel('qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-PRIONACE-ASSAY7.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'prio1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	prio2 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY8.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'prio2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	prio3 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY9.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'prio3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Alopias
	alo1 <- read_excel('qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-ALOPIAS-ASSAY10.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'alo1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	alo2 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY11.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'alo2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	alo3 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY12.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'alo3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Lamna - appears to be contaminated (the PPM, could remake but not enough sample volume. may re-run with aliquot later)
	#lam1 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'Unknown')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA, TRUE ~ Cq))%>%
		mutate(testID = 'lam1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	#lam2 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'Unknown')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA, TRUE ~ Cq))%>%
		mutate(testID = 'lam2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	#lam3 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'Unknown')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA, TRUE ~ Cq))%>%
		mutate(testID = 'lam3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 

	## one combo plate, all four species - for chapter 4 remove the BRUV ones
	combo1 <-  read_excel('qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY17.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'combo1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) %>%
		filter(str_detect(Sample.Name, 'M'))
	
	combo2 <-  read_excel('qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY18.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'combo')%>%
		left_join(., neb, by = c('Target.Name'='sp')) %>%
		filter(str_detect(Sample.Name, 'M'))


	combo3 <-  read_excel('qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY19.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role == 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'combo')%>%
		left_join(., neb, by = c('Target.Name'='sp')) %>%
		filter(str_detect(Sample.Name, 'M'))

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
	
		write.csv(data, 'qPCRresults/collatedqPCRresults_2024_cornwallspecies.csv')

#########
## - Step 2: calculate copies
#########

	## write function to calculate copies
 	 calculate_copies <- function(intercept, slope, Cq.adj) {
    	copies <- 10^((Cq.adj - intercept)/slope)
    	return(copies)
    	}

    ## apply function in 2 parts in accordance with LOQ and LOD criteria. (A) First to calculate tech rep average copies, and then (B) the Field Sample and Field Control average (average of tech reps), then (C) join dataframes. 
	# (A)
		output <- group_by(data, Target.Name, Sample.Name) %>%
		    mutate(., loq_check = ifelse(Cq <= LOQ & !is.na(Cq), 1, NA)) %>%
		    mutate(., loq_check = ifelse(any(loq_check == 1), 1, NA)) %>%
		    mutate(Cq = as.numeric(ifelse(loq_check == 1 & is.na(Cq), LOD, Cq)))%>%
		    mutate(., reliable = ifelse(loq_check == is.numeric(loq_check), TRUE))%>%
		    mutate(reliable = replace_na(reliable,FALSE))%>%
		    mutate(Cq.adj = case_when(reliable == 'TRUE' ~ as.numeric(Cq), reliable == 'FALSE' ~ LOD))%>%
		    mutate(copies = if_else(reliable, calculate_copies(intercept, slope, Cq.adj), NA_real_))%>%
		    mutate(copies.techrepavg = if_else(reliable, mean(copies), NA_real_))
		   summary(output)
	# (B) 
		fieldsamp_averagecopies <- output %>%
		  filter(grepl("\\.1$|\\.2$|\\.3$", Sample.Name)) %>%
		  group_by(Target.Name, Sample.Name) %>%
		  summarise(fieldsamples.copies = mean(copies.techrepavg, na.rm = TRUE))
		fieldsamp_averagecopies
		summary(fieldsamp_averagecopies)
	# (C) 
		output2 <- output %>%
		  left_join(fieldsamp_averagecopies, by = c("Target.Name", "Sample.Name")) %>%
		  mutate(copies.sampavg = if_else(grepl("\\.4$", Sample.Name), copies.techrepavg, fieldsamples.copies))%>%
		  dplyr::select(-fieldsamples.copies)
		output2%>%print(n=20)


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
			
			#st_write(mmm, 'eDNA_data_meta_and_qPCR_2023.shp', driver = 'ESRI Shapefile', append = FALSE)
			#st_write(mmm, 'eDNA_data_meta_and_qPCR_2023.csv', driver = 'CSV', append = FALSE)

	## join

	joined <- left_join(output2, meta.sf%>%mutate(eppendorfID = str_remove(eppendorfID, '[P]')), by = c('Sample.Name' = 'eppendorfID')) # this might have to be done with st_join. tbd. 

	joined
	joined%>%filter(is.na(sampleID))

#########
## - Step 4: save data 
#########

	saveRDS(joined, 'qPCRresults/processedQPCRresults_cornwallspecies_june2024.Rdata')

















