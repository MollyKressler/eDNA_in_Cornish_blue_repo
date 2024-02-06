
## Comparing Metaprobes to Water Bottles 
### using all species dataset 

## created by Molly M Kressler :: February 2024

## Load at Start
pacman::p_load(sf,tidyverse,lubridate,readr,readxl,lubridate,ggplot2,patchwork,cowplot,stringr)

copies <- read_csv("data_downloads_notpublic/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv")
ntc <- read_csv("data_downloads_notpublic/resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv", 
                col_types = cols(...1 = col_skip(), X = col_skip()))%>%
  mutate(Assay.Role = case_when(Assay.Role == 'Negative' ~ 'NTC',
                                Assay.Role == 'Positive' ~ 'Standard'
                                , .default = as.character(Assay.Role)))
 
##########
## - Separate the Field controls from the copies results
##########
## use 'copies', replicateID will include 'WBC' 

field.controls <- copies %>% filter(grepl('WBC', replicateID)) %>%
           mutate('Assay.Role' = 'Field.Control') # This is the field controls for the water bottles (and metaprobes from joint deployments)
field.controls 

samples <- copies %>% filter(!grepl('WBC', replicateID)) %>%
            mutate('Assay.Role' = 'Field.Sample') # This is all field samples, water bottles and metaprobes, without the field controls 


##########
## - Plots
##########

  ## Cq.value/mean by Sampling Event - MP and WB.
    ## geom_stat_summary, with ranges and medians for WBTs and MPs, and a star geom_point for WBC. These plots do NOT include where field samples NOR field controls had NAs. 
  #engraulis
  stat_sum_engraulis_methodtypes_byEVENT <- ggplot(data = samples) + stat_summary(aes(x = eventID, y = engraulis.cq.mean, color = methodtype),
      fun.min = function(z) { quantile(z,0.25) },
      fun.max = function(z) { quantile(z,0.75) },
      fun = median, position = position_dodge(0.6))+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      ylab('Cq (median, IQR)')+
      xlab('Sampling Event')+
      ylim(c(32,42))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      labs(color="Method")
    
  ## Copies over time - continuous, one line for MP, one line for WBT
    # engraulis

  ## Concentration over time - continuous, one line for MP, one line for WBT
    # engraulis


##########
## - Tables
##########
  ## Broad descriptive statistics - one table for whole project, 2 rows, 6 columns: samples sizes of each (sampling events and replicates), average concentration, lowest copy number non-zero, highest copy number non-zero, median copy number non-zero


##########
## - Statistical tests
##########








