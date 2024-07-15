## Pipeline to Calculate copy number for Negative and Positive Laboratory controls from qPCR assays

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
	neb <- read_excel('qPCRresults/standard_curve_equations_qPCR_2023species.xlsx', sheet = 'Sheet1', range = 'A1:I10',.name_repair = 'universal', col_types = c('text','text','numeric','numeric','numeric','numeric','date','numeric','numeric'))%>%
		mutate(dateofassay=as_date(dateofassay))%>%
		filter(dateofassay == '2024-06-28') 

	## Assays
	assay1 <- read_excel('scillies_data_edna/qPCR_results_scillies/08072024-KRESSLER-SCILLES-ASSAY1.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'assay1')%>%
		left_join(., neb, by = c('Target.Name'='sp'))
	assay2 <- read_excel('scillies_data_edna/qPCR_results_scillies/08072024-KRESSLER-SCILLES-ASSAY2.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'assay2')%>%
		left_join(., neb, by = c('Target.Name'='sp'))
	assay3 <- read_excel('scillies_data_edna/qPCR_results_scillies/10072024-KRESSLER-SCILLES-ASSAY3.xlsx', sheet = 'Results',range = 'A41:L137', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'UNKNOWN')%>%
		mutate(Cq = case_when(Cq == 'Undetermined' ~ NA_character_, TRUE ~ Cq))%>%
		mutate(testID = 'assay3')%>%
		left_join(., neb, by = c('Target.Name'='sp'))


	## combine into one df, all species, all assays

	data <- rbind(assay1, assay2, assay3)

	#########
	## save
	
		write.csv(data, 'scillies_data_edna/qPCR_results_scillies/collatedqPCRresults_2024_labcontrols_standards_cornwallspecies.csv')


#########
## - Step 2: calculate copies
#########

		data <- read.csv('scillies_data_edna/qPCR_results_scillies/collatedqPCRresults_2024_labcontrols_standards_cornwallspecies.csv')%>%
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
	
		write.csv(output3, 'scillies_data_edna/qPCR_results_scillies/processedqPCRresults_2024_labcontrols_standards_cornwallspecies.csv')


#########
## - Step 3: format for tables later
#########
  output3 <- read.csv('scillies_data_edna/qPCR_results_scillies/processedqPCRresults_2024_labcontrols_standards_cornwallspecies.csv')%>%
  	dplyr::select(-X)%>%
  	mutate(eventID = NA_character_,method.type = NA_character_,sampleID = NA_character_,dnaCont = case_when(Assay.Role == 'STANDARD' ~ Quantity, Assay.Role == 'NTC' ~ 0, FALSE ~ NA_real_),ratio260.280 = NA_real_,ratio260.230 = NA_real_, loq_check = 1, reliable = TRUE,copies.sampavg = copies,copies.techrepavg= copies)
 
  #########
	## save
	
		write.csv(output3, 'scillies_data_edna/qPCR_results_scillies/processedqPCRresults_2024_lab_controls_standards_formattedForTables_cornwallspecies.csv')






