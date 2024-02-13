## - Calculate copies of eDNA, per replicate sample, per species 
## Project: eDNA in Cornwall 

## created 07 December 2023
## by Molly M Kressler 

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/EDNA/data_edna/')
pacman::p_load(dplyr,lubridate,readr,readxl,lubridate,stringr)

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
	## need to separate out the standards and negative controls bc they get lost in the grouping and slicing because of duplicate names
	controls_engraulis <- all_engr%>%filter(Assay.Role != 'Unknown')
	all_eng <- all_engr %>% filter(Assay.Role == 'Unknown')


## row by row, write tidy function using NEBio standard curve equation (species specific), into new column Species.Copies, e.g. Engraulis.copies 

	# write a function to calculate copies with species-specific standard curve data 
	calculate_copies <- function(intercept, slope, Cq.Mean) {
		copies <- 10^((Cq.Mean - intercept)/slope)
		return(copies)
		}

	e <- all_eng %>%
		group_by(Sample.Name) %>%
		slice_head()%>% ## first: need to filter down to only one row for each replicate and use the Cq mean
		mutate(intercept = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(intercept)),
			slope = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(slope)), 
			r.squared = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(r.squared)), 
			efficiency = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(efficiency)))%>% # bring in the standard curve information and statistics
		mutate(engraulis.copies = calculate_copies(intercept, slope, Cq.Mean))%>% # calculate copies using linear formula
		replace_na(list(engraulis.copies = 0)) # where Cq.Mean was NA, and therefore copies is NA, replace with a 0. 
		
		stopifnot(nrow(e) == 100)
	
		summary(e) # NAs are fine in Cq.Mean and Cq because they are from samples that never amplified (i.e. they have no detectable traces of that species). In the pipe above, the final line adjusts the copies estimate to 0 for these cases.

	controls_engrauli <- controls_engraulis %>%
		mutate(intercept = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(intercept)),
			slope = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(slope)), 
			r.squared = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(r.squared)), 
			efficiency = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(efficiency)))%>% # bring in the standard curve information and statistics
		mutate(engraulis.copies = calculate_copies(intercept, slope, Cq.Mean))%>% # calculate copies using linear formula
		replace_na(list(engraulis.copies = 0)) # where Cq.Mean was NA, and therefore copies is NA, replace with a 0. 
		

## save as a new results file, because of the slicing  
	write.csv(e, 'copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')
	write.csv(controls_engrauli, 'resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv') ### NEED TO COMPILE FOR ALL SPP LIKE I DO FOR SAMPLES BELOW

 ## Add Test Assay name variables, e.g. 'TEST3'
	engraul <- read.csv('copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')
	
	assays <- as_tibble(engraul%>%distinct(Experiment)) %>%
		mutate(testID = str_extract(string = Experiment,pattern = 'TEST[0-9]'))%>%
		mutate(testID = replace_na(testID, 'TEST1')) # test1 experiment name is missing a T so this fills in the blank new variable Test1

	engraul2 <- left_join(engraul, assays, by='Experiment',relationship = 'many-to-one') # yep it worked.

	write.csv(engraul2, 'copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')

	## repeat for controls 
	controls <- read.csv('resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv')
	assays2 <- as_tibble(controls%>%distinct(Experiment)) %>%
		mutate(testID = str_extract(string = Experiment,pattern = 'TEST[0-9]'))%>%
		mutate(testID = replace_na(testID, 'TEST1')) # test1 experiment name is missing a T so this fills in the blank new variable Test1

	controls2 <- left_join(controls, assays2, by='Experiment',relationship = 'many-to-one') 

	write.csv(controls2, 'resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv') ### NEED TO COMPILE FOR ALL SPP LIKE I DO FOR SAMPLES BELOW


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


## import species copies_perReplicates dfs - select down to only the Sample.Name, Cq.Mean, species.copies. And link these replicates back to their eventID for spatial joining. 
	# engraulis 
	e <- read.csv('copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')%>%
		dplyr::select(Sample.Name, replicateID, Cq.Mean, engraulis.copies, testID)%>%
		rename(engraulis.cq.mean = Cq.Mean) %>%
		mutate(sampleID = str_remove(replicateID,'[.]\\w+|:'),.before=Sample.Name)%>%
		mutate(sampleID = str_replace_all(sampleID,'WBT1','WBT1.WBC1'))%>% # corrects WBT1.4s to WBT1.WBC1 in sampleID
		mutate(sampleID = str_replace_all(sampleID,'-WBC1','-WBT1.WBC1'))%>% # add WBT1.WBC1 if WBT 
		mutate(sampleID = str_replace_all(sampleID,' ','')) # removes erroneous spaces
	head(e)

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


## UPDATE THIS WITH MORE SPECIES: join allsppcop to spatial sf. keep all metadata so its all in one place 
	allsf <- st_as_sf(left_join(e,sf, by='sampleID', relationship = 'many-to-one')) # this join is not workinng
	allsf 

## UPDATE WITH MORE SPECIES, save as one shp, and csv
	write_sf(allsf, 'compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.shp', driver='ESRI Shapefile')
		# rename code: 
			#rename(sampleID = samplID,
				#Sample.Name = Smpl_Nm,
				#replicateID = rplctID,
				#engraulis.cq.mean = engrl__,
				#engraulis.copies = engrls_,
				#sampDATE = smpDATE,
				#methodtype = mthdtyp)
	st_write(allsf, 'compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv', driver = 'CSV')
	write.csv(allsppcop, 'copies_perReplicate_notStandards_allspecies_2023eDNACornishBlue.csv')
## Repeat above 3 steps for standards, test assay positive controls and test assay negative controls 



















## - Calculate copies of eDNA, per replicate sample, per species 
## Project: eDNA in Cornwall 

## created 07 December 2023
## by Molly M Kressler 

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/EDNA/data_edna/')
pacman::p_load(dplyr,lubridate,readr,readxl,lubridate,stringr)

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
	## need to separate out the standards and negative controls bc they get lost in the grouping and slicing because of duplicate names
	controls_engraulis <- all_engr%>%filter(Assay.Role != 'Unknown')
	all_eng <- all_engr %>% filter(Assay.Role == 'Unknown')


## row by row, write tidy function using NEBio standard curve equation (species specific), into new column Species.Copies, e.g. Engraulis.copies 

	# write a function to calculate copies with species-specific standard curve data 
	calculate_copies <- function(intercept, slope, Cq.Mean) {
		copies <- 10^((Cq.Mean - intercept)/slope)
		return(copies)
		}

	e <- all_eng %>%
		group_by(Sample.Name) %>%
		slice_head()%>% ## first: need to filter down to only one row for each replicate and use the Cq mean
		mutate(intercept = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(intercept)),
			slope = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(slope)), 
			r.squared = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(r.squared)), 
			efficiency = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(efficiency)))%>% # bring in the standard curve information and statistics
		mutate(engraulis.copies = calculate_copies(intercept, slope, Cq.Mean))%>% # calculate copies using linear formula
		replace_na(list(engraulis.copies = 0)) # where Cq.Mean was NA, and therefore copies is NA, replace with a 0. 
		
		stopifnot(nrow(e) == 100)
	
		summary(e) # NAs are fine in Cq.Mean and Cq because they are from samples that never amplified (i.e. they have no detectable traces of that species). In the pipe above, the final line adjusts the copies estimate to 0 for these cases.

	controls_engrauli <- controls_engraulis %>%
		mutate(intercept = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(intercept)),
			slope = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(slope)), 
			r.squared = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(r.squared)), 
			efficiency = as.numeric(neb %>% filter(sp == 'engraulis') %>% select(efficiency)))%>% # bring in the standard curve information and statistics
		mutate(engraulis.copies = calculate_copies(intercept, slope, Cq.Mean))%>% # calculate copies using linear formula
		replace_na(list(engraulis.copies = 0)) # where Cq.Mean was NA, and therefore copies is NA, replace with a 0. 
		

## save as a new results file, because of the slicing  
	write.csv(e, 'copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')
	write.csv(controls_engraulis, 'resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv') ### NEED TO COMPILE FOR ALL SPP LIKE I DO FOR SAMPLES BELOW

 ## Add Test Assay name variables, e.g. 'TEST3'
	engraul <- read.csv('copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')
	
	assays <- as_tibble(engraul%>%distinct(Experiment)) %>%
		mutate(testID = str_extract(string = Experiment,pattern = 'TEST[0-9]'))%>%
		mutate(testID = replace_na(testID, 'TEST1')) # test1 experiment name is missing a T so this fills in the blank new variable Test1

	engraul2 <- left_join(engraul, assays, by='Experiment',relationship = 'many-to-one') # yep it worked.

	write.csv(engraul2, 'copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')

	## repeat for controls 
	controls <- read.csv('resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv')
	assays2 <- as_tibble(controls%>%distinct(Experiment)) %>%
		mutate(testID = str_extract(string = Experiment,pattern = 'TEST[0-9]'))%>%
		mutate(testID = replace_na(testID, 'TEST1')) # test1 experiment name is missing a T so this fills in the blank new variable Test1

	controls2 <- left_join(controls2, assays2, by='Experiment',relationship = 'many-to-one') 

	write.csv(controls2, 'resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv') ### NEED TO COMPILE FOR ALL SPP LIKE I DO FOR SAMPLES BELOW


#################################
## - Combine multiple species copies datasets with spatial data and meta data
#################################

## import spatial metadata 
	sf <- st_as_sf(st_read('data_edna/eDNA_data_meta_and_qPCR_2023.shp'),crs='WGS84')%>%
		rename(timeIN=tm_tmIN,
			typeID=samplID,
			recBY=rcrdd_b,
			passengers=pssngrs,
			sampDAT=sampldt,
			sampDATE=smplng_,
			methodtype=mthd_ty)%>%	
		mutate(sampleID = paste0(eventID,'-',typeID),.before= swt_ID1)


## import species copies_perReplicates dfs - select down to only the Sample.Name, Cq.Mean, species.copies. And link these replicates back to their eventID for spatial joining. 
	# engraulis 
	e <- read.csv('data_edna/copies_perReplicate_notStandards_En.encras_TESTS_1to7_NOV2023.csv')%>%
		dplyr::select(Sample.Name, replicateID, Cq.Mean, engraulis.copies, testID, dnaCont)%>%
		rename(engraulis.cq.mean = Cq.Mean,engraulis.dnaCont = dnaCont) %>%
		mutate(sampleID = str_remove(replicateID,'[.]\\w+|:'),.before=Sample.Name)%>%
		mutate(sampleID = str_replace_all(sampleID,'WBT1','WBT1.WBC1'))%>% # corrects WBT1.4s to WBT1.WBC1 in sampleID
		mutate(sampleID = str_replace_all(sampleID,'-WBC1','-WBT1.WBC1'))%>% # add WBT1.WBC1 if WBT 
		mutate(sampleID = str_replace_all(sampleID,' ','')) # removes erroneous spaces
	head(e)

	# sprattus 
	t <- read.csv('')%>%
		dplyr::select(Sample.Name, replicateID, Cq.Mean, sprattus.copies)%>%
		rename(engraulis.cq.mean = Cq.Mean, engraulis.dnaCont = dnaCont) %>%
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


## UPDATE THIS WITH MORE SPECIES: join allsppcop to spatial sf. keep all metadata so its all in one place 
	allsf <- st_as_sf(left_join(e,sf, by='sampleID', relationship = 'many-to-one'))
	allsf 

## UPDATE WITH MORE SPECIES, save as one shp, and csv
	write_sf(allsf, 'data_edna/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.shp', driver='ESRI Shapefile')
		# rename code: 
			#rename(sampleID = samplID,
				#Sample.Name = Smpl_Nm,
				#replicateID = rplctID,
				#engraulis.cq.mean = engrl__,
				#engraulis.copies = engrls_,
				#sampDATE = smpDATE,
				#methodtype = mthdtyp)
	st_write(allsf, 'data_edna/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv', driver = 'CSV')
	write.csv(allsppcop, 'copies_perReplicate_notStandards_allspecies_2023eDNACornishBlue.csv')
## Repeat above 3 steps for standards, test assay positive controls and test assay negative controls 


















