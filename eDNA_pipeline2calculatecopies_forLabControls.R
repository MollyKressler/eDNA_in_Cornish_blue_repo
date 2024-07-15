## Pipeline to Calculate copy number for technical and sample replicates from assays of seawater samples analysed with qPCR but only the Negative and Positive Laboratory controls from these assays

## Molly M Kressler :: 28 May 2024

#########
## - Background
#########

## specifically written for .xlsx files as produced by Applied Biosciences Quant 7 Studio Flex qPCR machine, but could be adapted with modifications to initial data importing. 

## Calculates copies for laboratory controls (NTCs and Positives/Standards) with conditions regarding the limit of detection (LOD) and the limit of quantification (LOQ). The LOD and LOQ are identified in the Standard Curve Tests performed for each species and are stored with the standard curve equation data in a separate data file. 

## Definitions
	## LOD: limit of detection, lowest concentration (or, practically, the corresponding Cq) of standard template DNA with at least one amplified technical replicate.  Species specific, identified by standard curve tests with qPCR assay
	## LOQ: limit of quantification, lowest concentration (or, practcally, correspodning Cq) of standard template DNA with 3 positive technical replicates. Species specific, identified by standard curve tests with qPCR assay.
	## technical replicate, or tech rep: the replicates within a qPCR assay for each sampling replicate, n=3.
	## sample replicate, or samp rep: the replicates within a single sampling event, n=3 for metaprobe and water bottle sampling



#########
## - Broad steps
#########

## (1) Import and tidy datasets from qPCR assay output files, specifically segregating down to only Lab Controls (NTCs and Standards/Positive Controls). Done for each species, saving individual species files, before collating species and saving as one large file. 
## (2) Calculate copies per technical replicate, and per sampling replicate. 
## (3) Join to metadata and spatial data


#########
## - Step 0: Define and Load workspace
#########

pacman::p_load(sf,dplyr,lubridate,readr,readxl, tidyverse, stringr, flextable)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/EDNA/data_edna/')

#########
## - Step 1: tidying exports from qPCR machine
#########
	## Standard Curves
	neb <- read_excel('qPCRresults/standard_curve_equations_qPCR_2023species.xlsx', sheet = 'Sheet1', range = 'A1:I8',.name_repair = 'universal', col_types = c('text','text','numeric','numeric','numeric','numeric','date','numeric','numeric'))%>%
		mutate(dateofassay=as_date(dateofassay))%>%
		filter(year(dateofassay) == 2024) 

	## Engraulis
	eng1 <- read_excel('qPCRresults/ESI_assays_spring2024/05062024-KRESSLER-ENGRAULIS-ASSAY1.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'eng1')%>%
		left_join(., neb, by = c('Target.Name'='sp'))

	eng2 <- read_excel('qPCRresults/ESI_assays_spring2024/05062024-KRESSLER-ENGRAULIS-ASSAY2.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'eng2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	eng3 <- read_excel('qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-ENGRAULIS-ASSAY3.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'eng3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Scomber
	sco1 <- read_excel('qPCRresults/ESI_assays_spring2024/07062024-KRESSLER-SCOMBER-ASSAY4.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'sco1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	sco2 <- read_excel('qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-SCOMBER-ASSAY5.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'sco2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	sco3 <- read_excel('qPCRresults/ESI_assays_spring2024/11062024-KRESSLER-SCOMBER-ASSAY6.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'sco3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Prionace
	prio1 <- read_excel('qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-PRIONACE-ASSAY7.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'prio1', Target.Name = 'Prionace')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	prio2 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY8.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'prio2', Target.Name = 'Prionace')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	prio3 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-PRIONACE-ASSAY9.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'prio3', Target.Name = 'Prionace')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Alopias
	alo1 <- read_excel('qPCRresults/ESI_assays_spring2024/10062024-KRESSLER-ALOPIAS-ASSAY10.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'alo1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	alo2 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY11.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'alo2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	alo3 <- read_excel('qPCRresults/ESI_assays_spring2024/12062024-KRESSLER-ALOPIAS-ASSAY12.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'alo3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Lamna - appears to be contaminated (the PPM, could remake but not enough sample volume. may re-run with aliquot later)
	#lam1 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA, TRUE ~ Cq))%>%
		mutate(testID = 'lam1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	#lam2 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA, TRUE ~ Cq))%>%
		mutate(testID = 'lam2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	#lam3 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA, TRUE ~ Cq))%>%
		mutate(testID = 'lam3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 

	## one combo plate, all four species - for chapter 4 remove the BRUV ones
	combo1 <-  read_excel('qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY17.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
			rename(Cq = CT, Assay.Role = Task)%>%
			mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
			filter(Assay.Role != 'UNKNOWN')%>%
			mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'combo1', Target.Name = case_when(Target.Name == 'Prionance' ~ 'Prionace', Target.Name != 'Prionance' ~ Target.Name))%>%
			left_join(., neb, by = c('Target.Name'='sp')) 
	
	combo2 <-  read_excel('qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY18.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'combo2', Target.Name = case_when(Target.Name == 'Prionance' ~ 'Prionace', Target.Name != 'Prionance' ~ Target.Name))%>%
		left_join(., neb, by = c('Target.Name'='sp')) 

	combo3 <-  read_excel('qPCRresults/ESI_assays_spring2024/14062024-KRESSLER-MIXED-ASSAY19.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'combo3', Target.Name = case_when(Target.Name == 'Prionance' ~ 'Prionace', Target.Name != 'Prionance' ~ Target.Name))%>%
		left_join(., neb, by = c('Target.Name'='sp'))
	
	## combine into one df, all species, all assays

	data <- rbind(eng1,eng2,eng3,sco1,sco2,sco3,prio1, prio2, prio3,alo1,alo2,alo3,combo1,combo2,combo3)
	summary(data)

	table <- data %>% 
		mutate(Cq = round(as.numeric(Cq),3),Ct.Mean = round(as.numeric(Ct.Mean),3),Ct.SD = round(as.numeric(Ct.SD),3),Quantity = round(as.numeric(Quantity),3))%>%
		dplyr::select(-Omit, -Well, -Assay.Role, -species, -dateofassay, -slope, -intercept, -efficiency, -LOD, -LOQ, -Quantity, -r.squared)%>%
		flextable()%>%
		autofit()%>%
		theme_zebra()

		 save_as_docx(table, path = 'qPCRresults/labcontrols_standards_cornwalledna_amprawdata.docx') 

	#########
	## save
	
		write.csv(data, 'qPCRresults/collatedqPCRresults_2024_labcontrols_standards_cornwallspecies.csv')


#########
## - Step 2: calculate copies
#########

		data <- read.csv('qPCRresults/collatedqPCRresults_2024_labcontrols_standards_cornwallspecies.csv')%>%
			dplyr::select(-X)

  calculate_copies <- function(intercept, slope, Cq.adj) {
    copies <- 10^((Cq.adj - intercept)/slope)
    return(copies)
    }

    ## Don't need to run.
			c <- read.csv('testset_of_controls_for_writing_controlsPipe_may2024.csv')%>%
	      dplyr::select(-X)%>%
	      mutate(Assay.Role = case_when(Assay.Role == 'Negative' ~ 'NTC', Assay.Role != 'Negative' ~ Assay.Role))

	output3 <- group_by(data, Target.Name, testID) %>%
  	mutate(Cq.adj = as.numeric(ifelse(Assay.Role == 'NTC' & Cq >= LOQ, NA_real_, Cq)))%>% 
  	mutate(copies = calculate_copies(intercept, slope, Cq.adj))

  summary(output3)

  #########
	## save
	
		write.csv(output3, 'qPCRresults/processedqPCRresults_2024_labcontrols_standards_cornwallspecies.csv')


#########
## - Step 3: format for tables later
#########
  output3 <- read.csv('qPCRresults/processedqPCRresults_2024_labcontrols_standards_cornwallspecies.csv')%>%
  	dplyr::select(-X)%>%
  	mutate(eventID = NA_character_,method.type = NA_character_,sampleID = NA_character_,dnaCont = case_when(Assay.Role == 'STANDARD' ~ Quantity, Assay.Role == 'NTC' ~ 0, FALSE ~ NA_real_),ratio260.280 = NA_real_,ratio260.230 = NA_real_, loq_check = 1, reliable = TRUE,copies.sampavg = copies,copies.techrepavg= copies)
 
  #########
	## save
	
		write.csv(output3, 'qPCRresults/processedqPCRresults_2024_lab_controls_standards_formattedForTables_cornwallspecies.csv')






