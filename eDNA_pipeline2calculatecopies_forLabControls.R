## Pipeline to Calculate copy number for technical and sample replicates from assays of seawater samples analysed with qPCR but only the Negative and Positive Laboratory controls from these assays

## Molly M Kressler :: 28 May 2024

#########
## - Background
#########

## specifically written for xlsx files as produced by Applied Biosciences Quant 7 Studio Flex qPCR machine, but could be adapted with modifications to initial data importing. 

## Caculates copies for laboratory controls (NTCs and Positives/Standards) with conditions regarding the limit of detection (LOD) and the limit of quantification (LOQ). The LOD and LOQ are identified in the Standard Curve Tests performed for each species and are stored with the standard curve equation data in a separate data file. 

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

pacman::p_load(sf,dplyr,lubridate,readr,readxl)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/EDNA/data_edna/')

#########
## - Step 1: tidying exports from qPCR machine
#########
	## Standard Curves
	neb <- read_excel('qPCRresults/standard_curve_equations_qPCR_2023species.xlsx',.name_repair = 'universal', col_types = c('text','text','numeric','numeric','numeric','numeric','date'))%>%
		mutate(dateofassay=as_date(dateofassay), LOQ = as.numeric(LOQ),LOD = as.numeric(LOD))
	neb

	## Engraulis
	eng1 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'eng1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	eng2 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'eng2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	eng3 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'eng3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Scomber
	sco1 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'sco1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	sco2 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'sco2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	sco3 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'sco3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Alopias
	alo1 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'alo1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	alo2 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'alo2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	alo3 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'alo3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## Lamna 
	lam1 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'lam1')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	lam2 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'lam2')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	lam3 <- read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'lam3')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	## one combo plate, all four species
	combo <-  read_excel('qPCRresults/ESI_assays_spring2024/.xlsx', sheet = 'Results',range = 'A41:S104', .name_repair = 'universal', col_types='text')%>%
		rename(Cq = CT, Assay.Role = Task)%>%
		mutate(Sample.Name = case_when(!is.na(Sample.Name) ~ Sample.Name, is.na(Sample.Name) ~ paste0(Assay.Role,'_',Quantity)))%>%
		filter(Assay.Role != 'Unknown')%>%
		mutate(testID = 'combo')%>%
		left_join(., neb, by = c('Target.Name'='sp')) 
	
	## combine into one df, all species, all assays

	data <- rbind(eng1,eng2,eng3,sco1,sco2,sco3,alo1,alo2,alo3,lam1,lam2,lam3,combo)
	stopifnot(ncol(data) = 28)


#########
## - Step 2: calculate copies
#########

  calculate_copies <- function(intercept, slope, Cq.adj) {
    copies <- 10^((Cq.adj - intercept)/slope)
    return(copies)
    }


c <- read.csv('testset_of_controls_for_writing_controlsPipe_may2024.csv')%>%
      dplyr::select(-X)%>%
      mutate(Assay.Role = case_when(Assay.Role == 'Negative' ~ 'NTC', Assay.Role != 'Negative' ~ Assay.Role))
c

output <- group_by(c, Target.Name, testID) %>%
  mutate(Cq.adj = ifelse(Assay.Role == 'NTC' & Cq >= LOQ, NA_real_, Cq))%>% # have to add this line in for the cases when NTC's (or anything) self-prime past the LOQ
  mutate(copies = calculate_copies(intercept, slope, Cq.adj))


#########
## - Step 3: join to metadata and spatial information
#########

	## Import metadata 

	meta.sf <- read.csv('eDNA_data_meta_and_qPCR_2023.csv')%>%
		dplyr::select(-X)%>%
		st_as_sf(., coords=c('Longitude', 'Latitude'), crs = 'WGS84')

	## join

	joined <- left_join(output, meta.sf, by = c(Sample.Name = eppendorfID)) # this might have to be done with st_join. tbd. 
	stopifnot(nrow(joined = nrow(ouput)))











