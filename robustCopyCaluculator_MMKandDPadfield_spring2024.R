## Hello Daniel 

library(tidyverse)

## Little dataset
rep <- c(10.2, 10.2,10.2,10.2,10.3,10.3,10.3,10.3,9.1,9.1,9.1,9.1)
well <- c('A1','A2','A3','A4','B1','B2','B3','B4','C1','C2','C3','C4')
Cq <- c(NA,NA,32,29,NA,15,16,14,NA,NA,NA,34)

data <- cbind(rep, well, Cq)
	## each sample (e.g. 10.2) has four technical replicates, identifiable by well-ID.
	## the Cq is the cycle at which the well amplified.

## Quantiative Standards
LOD = 35
LOQ = 29

## Objective: when Cq == NA set the Cq equal to the LOD IF the Cqs of any other technical replicates are !=NA and ≤LOQ (less than because lower Cq is earlier amplification). Even one that meets the non-NA & ≤LOQ criteria is enough, like with tech reps for 10.2. 
## Where technical replicates Cqs != NA and ≥LOQ, i.e. here any value above 29, ignore/don't mutate the Cq value.


## thanks for the help.  

# 1. convert data to dataframe
data <- data.frame(data)
head(data)
str(data)

# 2. When Cq == NA - set the Cq equal to the LOD if the CQs of any other technical replicates are !=NA and ≤LOQ

# make a column to check whether any of the values are <LOQ
# Give them a 1 and if not replace with NA
data <- group_by(data, rep) %>%
  mutate(., loq_check = ifelse(Cq <= LOQ & !is.na(Cq), 1, NA)) %>%
  ungroup()

# where a 1 exists replace all NAs with 1
data <- group_by(data, rep) %>%
  mutate(., loq_check = ifelse(any(loq_check == 1), 1, NA)) %>%
  ungroup()

# where loq check == 1 and Cq == NA replace with LOD
data <- mutate(data, Cq = ifelse(loq_check == 1 & is.na(Cq), LOD, Cq))

## Molly - 2-3 may 2024

# 3. Where the sampling replicate is unreliable, i.e. it has NAs still after this (2) step, it needs to not be used in the averaging. 

data2 <- group_by(data, rep) %>%
  mutate(., reliable = ifelse(loq_check == is.numeric(loq_check), TRUE))%>%
  mutate(reliable = replace_na(reliable,FALSE))%>%
  mutate(Cq.adj = case_when(reliable == 'TRUE' ~ as.numeric(Cq), reliable == 'FALSE' ~ 0))%>%
  ungroup()
data2

# 4. Next step would be calculating the copy numbers for each tech rep. 

  calculate_copies <- function(intercept, slope, Cq.adj) {
    copies <- 10^((Cq.adj - intercept)/slope)
    return(copies)
    }

slope = -5
intercept = 0.05

data3 <- group_by(data2, rep) %>%
 mutate(copies = if_else(reliable, calculate_copies(intercept, slope, Cq.adj), NA_real_))
data3

# 5. Calculate sampling replicate average copy number 

data4 <- group_by(data3, rep) %>%
  mutate(copies.avg = if_else(reliable, mean(copies), NA_real_))
data4





##############
## - Make it pipe friendly
##############
## stes 1-5 above achieve the target: calculate copies for tech reps based on sequential conditions and their outcomes.
  ## (a) that non-amplified tech reps are set to the LOD when the Cqs of any other tech reps (within he same sampling rep) have amplified and are less than or equalt to the LOQ. This is called the LOQ-check.
  ## (b) determine if a tech rep is reliable based on the LOQ-check. Assign the Cq based on reliability: where TRUE, the Cq is not changed; where FALSE, the Cq is set to 0. 
  ## (c) calculate copies based on the Cq.adj using the formula. 

## define the calculate copies function before running pipe.
  calculate_copies <- function(intercept, slope, Cq.adj) {
    copies <- 10^((Cq.adj - intercept)/slope)
    return(copies)
    }
    ## it would pull values from the NEB df but her we'll feed it values. 
    slope = -5
    intercept = 0.05

rep <- c(10.1,10.1,10.1,10.1,10.2, 10.2,10.2,10.2,10.3,10.3,10.3,10.3,10.4,10.4,10.4,10.4)
well <- c('A1','A2','A3','A4','B1','B2','B3','B4', 'C1','C2', 'C3', 'C4', 'D1','D2', 'D3', 'D4')
Cq <- c(NA,NA,32,29,NA,15,16,14,NA,26,28,25,NA,14,18,15)
Target.Name = rep('sp1',16)
sp1 <- cbind(rep, well, Cq, Target.Name)
Target.Name = rep('sp2',16)
sp2 <- cbind(rep, well, Cq*.9, Target.Name)
Target.Name = rep('sp3',16)
sp3 <- cbind(rep, well, Cq*.5, Target.Name)
data <- rbind(sp1, sp2, sp3)

input <- as_tibble(data.frame(data))

output <- group_by(input, Target.Name, rep) %>%
    mutate(samp = sub("\\..*", "", rep)) %>%
    mutate(., loq_check = ifelse(Cq <= LOQ & !is.na(Cq), 1, NA)) %>%
    mutate(., loq_check = ifelse(any(loq_check == 1), 1, NA)) %>%
    mutate(Cq = ifelse(loq_check == 1 & is.na(Cq), LOD, Cq))%>%
    mutate(., reliable = ifelse(loq_check == is.numeric(loq_check), TRUE))%>%
    mutate(reliable = replace_na(reliable,FALSE))%>%
    mutate(Cq.adj = case_when(reliable == 'TRUE' ~ as.numeric(Cq), reliable == 'FALSE' ~ 0))%>%
    mutate(copies = if_else(reliable, calculate_copies(intercept, slope, Cq.adj), NA_real_))%>%
    group_by(Target.Name,rep)%>%
    mutate(copies.techrepavg = if_else(reliable, mean(copies), NA_real_))# calculates the average for each sampling replicate (mean(tech rep 1, tech rep 2, tech rep 3))
output

# Calculate the average value for Field samples (reps ending with .1, .2., (and .3 in real data but not here))
fieldsamp_averagecopies <- output %>%
  filter(grepl("\\.1$|\\.2$|\\.3$", rep)) %>%
  group_by(Target.Name, samp) %>%
  summarise(fieldsamples.copies = mean(copies.techrepavg, na.rm = TRUE))
fieldsamp_averagecopies

# Join the average values back to the processed data and carry over the Field Control average from the average of the technical replicates.
output2 <- output %>%
  left_join(fieldsamp_averagecopies, by = c("Target.Name", "samp")) %>%
  mutate(copies.sampavg = if_else(grepl("\\.4$", rep), copies.techrepavg, fieldsamples.copies))%>%
  select(-fieldsamples.copies)
output2%>%print(n=20)

## the next step would be to calculate the sample ID average copy number. 

# data pre-requisites
  ## data be in a data frame
  ## has the following columns: sample replicate ID number, Cq. 
  ## tech reps are in rows, and don't need an identifying ID.



##### 
## - controls need their own version because their structure is different
#####
## Still need the LOQ checks but you don't need the reliability testing. Could leave it in for checks. But in this case, if the NTC comes back as 0, then it is 0 because we KNOW it should be 0. Need to modify. 
## modify for NTCs: where Assay.Role == NTC & Cq != NA & Cq >= LOQ, set to 0.

c <- read.csv('testset_of_controls_for_writing_controlsPipe_may2024.csv')%>%
      dplyr::select(-X)%>%
      mutate(Assay.Role = case_when(Assay.Role == 'Negative' ~ 'NTC', Assay.Role != 'Negative' ~ Assay.Role))
c

output <- group_by(c, Target.Name, testID) %>%
  mutate(Cq.adj = ifelse(Assay.Role == 'NTC' & Cq >= LOQ, NA_real_, Cq))%>% # have to add this line in for the cases when NTC's (or anything) self-prime past the LOQ
  mutate(copies = calculate_copies(intercept, slope, Cq.adj))

output # yep!
















