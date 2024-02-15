## - Calculate copies of eDNA, per replicate sample, per species 
## Project: eDNA in Cornwall 

## created 07 December 2023
## by Molly M Kressler 

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/EDNA/data_edna/')
pacman::p_load(tidyverse,lubridate,readr,readxl,lubridate,stringr)

#################################
## - Tidy NEBio Calculator results from Standard Curve Tests  
#################################

## import excel sheet with numbers 

	neb <- read_excel('qPCRresults/standard_curve_equations_qPCR_2023species.xlsx',.name_repair = 'universal', col_types = c('text','text','numeric','numeric','numeric','numeric','date'))%>%
		mutate(dateofassay=as_date(dateofassay))
	neb

#################################
## - Calculate Copies   
#################################

## import species tidyed qPCR results 
	all_engr <- read.csv('qPCRresults/2023Engraulisencrasicolus/tidyed_results_En.encras_TESTS_1to7_NOV2023.csv')%>%
		as_tibble()%>%
		mutate(Cq.Mean = as.numeric(Cq.Mean))%>%
		dplyr::select(-X)
	all_engr

## define Test Assay quick-names and add it to all_engr (samples & standards/controls dataframe)
	## engraulis
	assays <- as_tibble(all_engr%>%distinct(Experiment)) %>%
	mutate(testID = str_extract(string = Experiment,pattern = 'TEST[0-9]'))%>%
	mutate(testID = replace_na(testID, 'TEST1')) # test1 experiment name is missing a T so this fills in the blank new variable Test1
	
	all_engr <- left_join(all_engr, assays, by='Experiment',relationship = 'many-to-one')


## Separate Standards/NTCs from Field Samples and Controls
	## engraulis
	controls_engraulis <- all_engr%>%filter(Assay.Role != 'Unknown')%>%
		mutate(species = 'engraulis')
	all_eng <- all_engr %>% filter(Assay.Role == 'Unknown')%>%
		mutate(species = 'engraulis')


## row by row, write tidy function using NEBio standard curve equation (species specific), into new column Species.Copies, e.g. Engraulis.copies. ## Updating 15 Feb 2024 based on readings. New method: (1) calculate copies for each technical rep, save this df separately; (2) take the mean copy number for all the technical reps per replicate sample; (3) Calculate copy number for standards and NTCs 

	## create a function to calculate copies with species-specific standard curve data 
	calculate_copies <- function(intercept, slope, Cq.Mean) {
		copies <- 10^((Cq.Mean - intercept)/slope)
		return(copies)
		}

	## (1) Calculate copies for each technical replicate of replicate samples, use all_eng df. As of 15/2/24, where Cq = NA, set this equal to 0. 
		e.techreps <- all_eng %>%
		group_by(Sample.Name) %>%
		mutate(intercept = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(intercept)),
			slope = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(slope)), 
			r.squared = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(r.squared)), 
			efficiency = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(efficiency)))%>% # bring in the standard curve information and statistics
		mutate(engraulis.copies = calculate_copies(intercept, slope, Cq.Mean))%>% # calculate copies using linear formula
		replace_na(list(engraulis.copies = 0)) # where Cq.Mean was NA, and therefore copies is NA, replace with a 0. 
		
		summary(e.techreps) 

		## save df with copy number for each techincal replicate 
			write.csv(e.techreps, 'copies_perTechnicalReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')

	## (2) Calculate the mean copy number across technical replicates for samples. 
		e.repSamps <- e.techreps %>%
			group_by(Sample.Name)%>%
			summarise(mean.engraulis.copies = mean(engraulis.copies))%>%
			left_join(., all_eng, by='Sample.Name', multiple = 'first', relationship = 'many-to-one') # multiple argument says only keep the first join-match - this removes the technical replicates.
		e.repSamps

		## save df with mean copy number for each replicateID
			write.csv(e.repSamps, 'copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')
	
	## (3) Calculate the copy number for Standards and NTCs 
		controls_engrauli <- controls_engraulis %>%
			mutate(intercept = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(intercept)),
				slope = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(slope)), 
				r.squared = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(r.squared)), 
				efficiency = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(efficiency)))%>% # bring in the standard curve information and statistics
			mutate(engraulis.copies = calculate_copies(intercept, slope, Cq.Mean))%>% # calculate copies using linear formula
			replace_na(list(engraulis.copies = 0)) # where Cq.Mean was NA, and therefore copies is NA, replace with a 0. 
			
		## save
			write.csv(controls_engrauli, 'resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv') 



#################################
## - Combine multiple species copies datasets with spatial data and meta data
#################################

## import spatial metadata 
	sf <- st_as_sf(st_read('eDNA_data_meta_and_qPCR_2023.shp'),crs='WGS84')%>%
		rename(timeIN=tm_tmIN,
			typeID=samplID,
			recBY=rcrdd_b,
			passengers=pssngrs,
			sampDAT=sampldt,
			sampDATE=smplng_,
			methodtype=mthd_ty)%>%	
		mutate(sampleID = paste0(eventID,'-',typeID),.before= swt_ID1)


## import species copies_perReplicates & copies per technical replicates dfs - select down to only the Sample.Name, Cq.Mean, species.copies. And link these replicates back to their eventID for spatial joining. 
	# engraulis 
	e <- read.csv('copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')%>%
		dplyr::select(Sample.Name, sampleID, Cq.Mean, mean.engraulis.copies, testID, species)%>%
		rename(replicateID = sampleID,engraulis.cq.mean = Cq.Mean) %>%
		mutate(sampleID = str_remove(replicateID,'[.]\\w+|:'),.before=Sample.Name)%>%
		mutate(sampleID = str_replace_all(sampleID,'WBT1','WBT1.WBC1'))%>% # corrects WBT1.4s to WBT1.WBC1 in sampleID
		mutate(sampleID = str_replace_all(sampleID,'-WBC1','-WBT1.WBC1'))%>% # add WBT1.WBC1 if WBT 
		mutate(sampleID = str_replace_all(sampleID,' ','')) # removes erroneous spaces
	head(e)

	e.techs <- read.csv('copies_perTechnicalReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')%>%
		dplyr::select(Sample.Name, sampleID, Cq, Cq.Mean, engraulis.copies, testID, species)%>%
		rename(replicateID = sampleID,engraulis.cq.mean = Cq.Mean, engraulis.Cq = Cq) %>%
		mutate(sampleID = str_remove(replicateID,'[.]\\w+|:'),.before=Sample.Name)%>%
		mutate(sampleID = str_replace_all(sampleID,'WBT1','WBT1.WBC1'))%>% # corrects WBT1.4s to WBT1.WBC1 in sampleID
		mutate(sampleID = str_replace_all(sampleID,'-WBC1','-WBT1.WBC1'))%>% # add WBT1.WBC1 if WBT 
		mutate(sampleID = str_replace_all(sampleID,' ','')) # removes erroneous spaces
		head(e.techs)

	# sprattus 
	t <- read.csv('')%>%
		dplyr::select(Sample.Name, replicateID, Cq.Mean, sprattus.copies)%>%
		rename(engraulis.cq.mean = Cq.Mean) %>%
		mutate(sampleID = str_remove(replicateID,'[.]\\w+|:'),.before=Sample.Name)%>%
		mutate(sampleID = str_replace_all(sampleID,'WBT1','WBT1.WBC1'))%>% # corrects WBT1.4s to WBT1.WBC1 in sampleID
		mutate(sampleID = str_replace_all(sampleID,'-WBC1','-WBT1.WBC1'))%>% # add WBT1.WBC1 if WBT 
		mutate(sampleID = str_replace_all(sampleID,' ',''))
	# scombrus 
	m  <- read.csv('')%>%
		dplyr::select(Sample.Name, replicateID, Cq.Mean, scombrus.copies)%>%
		rename(engraulis.cq.mean = Cq.Mean) %>%
		mutate(sampleID = str_remove(replicateID,'[.]\\w+|:'),.before=Sample.Name)%>%
		mutate(sampleID = str_replace_all(sampleID,'WBT1','WBT1.WBC1'))%>% # corrects WBT1.4s to WBT1.WBC1 in sampleID
		mutate(sampleID = str_replace_all(sampleID,'-WBC1','-WBT1.WBC1'))%>% # add WBT1.WBC1 if WBT 
		mutate(sampleID = str_replace_all(sampleID,' ',''))
	# lamna
	l  <- read.csv('')%>%
		dplyr::select(Sample.Name, replicateID, Cq.Mean, lamna.copies)%>%
		rename(engraulis.cq.mean = Cq.Mean) %>%
		mutate(sampleID = str_remove(replicateID,'[.]\\w+|:'),.before=Sample.Name)%>%
		mutate(sampleID = str_replace_all(sampleID,'WBT1','WBT1.WBC1'))%>% # corrects WBT1.4s to WBT1.WBC1 in sampleID
		mutate(sampleID = str_replace_all(sampleID,'-WBC1','-WBT1.WBC1'))%>% # add WBT1.WBC1 if WBT 
		mutate(sampleID = str_replace_all(sampleID,' ',''))
	# alopias
	a  <- read.csv('')%>%
		dplyr::select(Sample.Name, replicateID, Cq.Mean, alopias.copies)%>%
		rename(engraulis.cq.mean = Cq.Mean) %>%
		mutate(sampleID = str_remove(replicateID,'[.]\\w+|:'),.before=Sample.Name)%>%
		mutate(sampleID = str_replace_all(sampleID,'WBT1','WBT1.WBC1'))%>% # corrects WBT1.4s to WBT1.WBC1 in sampleID
		mutate(sampleID = str_replace_all(sampleID,'-WBC1','-WBT1.WBC1'))%>% # add WBT1.WBC1 if WBT 
		mutate(sampleID = str_replace_all(sampleID,' ',''))
	# prionace
	p  <- read.csv('')%>%
		dplyr::select(Sample.Name, replicateID, Cq.Mean, prionace.copies)%>%
		rename(engraulis.cq.mean = Cq.Mean) %>%
		mutate(sampleID = str_remove(replicateID,'[.]\\w+|:'),.before=Sample.Name)%>%
		mutate(sampleID = str_replace_all(sampleID,'WBT1','WBT1.WBC1'))%>% # corrects WBT1.4s to WBT1.WBC1 in sampleID
		mutate(sampleID = str_replace_all(sampleID,'-WBC1','-WBT1.WBC1'))%>% # add WBT1.WBC1 if WBT 
		mutate(sampleID = str_replace_all(sampleID,' ',''))

## UPDATE THIS WITH MORE SPECIES - join dfs - want one row per replicate - 'allsppcop'
	allsppcop <-left_join(e,t,m,l,a,p ,by ='sampleID')

	allsppcop.techs <-left_join(e.techs,t.techs,m.techs,l.techs,a.techs,p.techs ,by ='sampleID')


## UPDATE THIS WITH MORE SPECIES: join allsppcop to spatial sf. keep all metadata so its all in one place 
	allsf <- st_as_sf(left_join(e,sf, by='sampleID', relationship = 'many-to-one'))
	allsf 

	allsf.techs <- st_as_sf(left_join(e.techs,sf, by='sampleID', relationship = 'many-to-one'))
	allsf.techs 

## UPDATE WITH MORE SPECIES, save as one shp, and csv
	write_sf(allsf, 'compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.shp', driver='ESRI Shapefile')
	write_sf(allsf.techs, 'compiledspeciesqPCR_copiesTECHNICALreps_withmetadata_2023eDNACornishBlue.shp', driver='ESRI Shapefile')
		# rename code: 
			#rename(sampleID = samplID,
				#Sample.Name = Smpl_Nm,
				#replicateID = rplctID,
				#engraulis.cq.mean = engrl__,
				#engraulis.copies = engrls_,
				#sampDATE = smpDATE,
				#methodtype = mthdtyp)
	st_write(allsf, 'compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv', driver = 'CSV')
	st_write(allsf.techs, 'compiledspeciesqPCR_copiesTECHNICALreps_withmetadata_2023eDNACornishBlue.csv', driver = 'CSV')
	write.csv(allsppcop, 'copies_perReplicate_notStandards_allspecies_2023eDNACornishBlue.csv')
## Repeat above 3 steps for standards, test assay positive controls and test assay negative controls 




#################################
## - Calculate sampleID copy numbers (mean of replicate samples values)
#################################

	copies <- read_csv("data_edna/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv")
	ntc <- read_csv("data_edna/resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv", 
	                col_types = cols(...1 = col_skip(), X = col_skip()))%>%
	  mutate(Assay.Role = case_when(Assay.Role == 'Negative' ~ 'NTC',
	                                Assay.Role == 'Positive' ~ 'Standard'
	                                , .default = as.character(Assay.Role)))

	## objective: copies per method per sampling event - UPDATE WITH MORE SPECIES

	c2 <- copies %>%
		group_by(sampleID)%>% #group the sampling events together, sampleID also has the method information
		summarise(sample.mean.engraulis.copies = mean(mean.engraulis.copies, na.rm = TRUE), sample.mean.engraulis.Cq = mean(engraulis.cq.mean, na.rm = TRUE))%>% # calculate the average of the sampleIDs mean.engraulis.copies. Add species to this as you get the data as successive X = mean(xx). 
		left_join(.,copies, by = 'sampleID', relationship = 'one-to-many', multiple = 'first')%>%
		dplyr::select(-mean.engraulis.copies, -engraulis.cq.mean, -replicateID, -Sample.Name)  # remove columns for data that is specific to singlular replicates of a sampling event. 
	c2 # this is now a dataframe with the sampleID and the sample.mean copy number for the species, with the metadata left_joined 

	
	## join with spatial info. UPDATE WITH MORE SPECIES 
	sf <- st_as_sf(st_read('eDNA_data_meta_and_qPCR_2023.shp'),crs='WGS84')%>%
		rename(timeIN=tm_tmIN,
			typeID=samplID,
			recBY=rcrdd_b,
			passengers=pssngrs,
			sampDAT=sampldt,
			sampDATE=smplng_,
			methodtype=mthd_ty)%>%	
		mutate(sampleID = paste0(eventID,'-',typeID),.before= swt_ID1)

	compiledreps <- st_as_sf(left_join(c2,sf%>%dplyr::select(sampleID,geometry), by='sampleID', relationship = 'many-to-one'))
	compiledreps

	## write and save

	write_sf(c2, 'compiledspeciesqPCR_SamplingEventCompiledCopiesandCq_withmetadata_2023eDNACornishBlue.shp', driver='ESRI Shapefile')

	st_write(compiledreps, 'compiledspeciesqPCR_SamplingEventCompiledCopiesandCq_withmetadata_2023eDNACornishBlue.csv', driver = 'CSV')












