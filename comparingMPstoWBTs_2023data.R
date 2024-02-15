
## Comparing Metaprobes to Water Bottles 
### using all species dataset 

## created by Molly M Kressler :: February 2024

##################
##########
## - Load at start 
##########
pacman::p_load(sf,tidyverse,lubridate,readr,readxl,lubridate,ggplot2,patchwork,cowplot,stringr,flextable,rstatix,lme4, modelsummary)

## For local R, not for server RStudio

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/EDNA')

samples <- read_csv('data_edna/compiledspeciesqPCR_SamplingEventCompiledCopiesandCq_withmetadata_2023eDNACornishBlue.csv') # copies and Cq per methodtype per sampling event 
copies <- read_csv("data_edna/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv") # copies and Cq per replicate 
tech.copies <- read_csv("data_edna/compiledspeciesqPCR_copiesTECHNICALreps_withmetadata_2023eDNACornishBlue.csv") # copies and Cq per technical replicate
ntc <- read_csv("data_edna/resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv", 
                col_types = cols(...1 = col_skip(), X = col_skip()))%>%
  mutate(Assay.Role = case_when(Assay.Role == 'Negative' ~ 'NTC',
                                Assay.Role == 'Positive' ~ 'Standard'
                                , .default = as.character(Assay.Role)))
 ## standard curve test assays - per species 
   curves <- read_excel('data_edna/qPCRresults/standard_curve_equations_qPCR_2023species.xlsx')

    engr.standards <- read_csv('data_edna/qPCRresults/2023Engraulisencrasicolus/tiyded_results_En.encras_STANDARDCURVETEST_NOV24_2023.csv',col_types = cols(...1 = col_skip()))
    engr.standards

    prio.standards <- read_csv('data_edna/qPCRresults/2023Prinoaceglauca/Pglauca_STANDARDCURVETEST_JUL312023.csv',col_types = cols(...1 = col_skip()))
    prio.standards

## For server RStudio

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
techreps.field.controls <- tech.copies %>% filter(grepl('WBC', replicateID)) %>%
  mutate('Assay.Role' = 'Field.Control') 
techreps.field.controls 

repsamples <- copies %>% filter(!grepl('WBC', replicateID)) %>%
  mutate('Assay.Role' = 'Field.Sample') # This is all field samples, water bottles and metaprobes, without the field controls 
techreps.samples <- tech.copies %>% filter(!grepl('WBC', replicateID)) %>%
  mutate('Assay.Role' = 'Field.Sample')


##########
## - Taxonomic groups, copies
##########

  by.taxa <- repsamples%>%
    rowwise()%>%
    mutate(allspp.copies = sum(c_across(c(engraulis.copies,sprattus.copies,scombrus.copies,lamna.copies,alopias.copies,prionace.copies)), na.rm = TRUE),
      fish.copies = sum(c_across(c(engraulis.copies,sprattus.copies,scombrus.copies)), na.rm = TRUE),
      shark.copies = sum(c_across(c(lamna.copies,alopias.copies,prionace.copies)), na.rm = TRUE))

##################

##########
## - Plots
##########

#### Concentration over time - continuous, one line for MP, one line for WBT - replicate repsamples data
  # all species
  dnacont_overtime<- ggplot(repsamples, aes(sampDATE,engraulis.dnaCont, group = methodtype, col = methodtype, fill = methodtype))+
    geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
    scale_color_manual(values=c('#DAA507','#8EC7D2'))+
    scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
    xlab('Sampling Date (2023)')+ 
    ylab('Nanodrop Concentration')+
    theme_bw()+
    theme(plot.title = element_text(size=10, face='italic'))+
    labs(fill="Method", col= "Method")
  
#### Standards Curve Test results - present them with Log() of Quantity and then with a table of the slope, intercept, R squared, and efficiency
  ## engraulis
    curve.engr <-  ggplot(data = engr.standards, aes(log(Quantity),Cq)) + 
      geom_point()+
      xlab('Quantity (copies/µL, log)')+
      geom_abline()+
      theme_bw()+
      ggtitle('E. encrasicolus')+
      theme(plot.title = element_text(size=10, face='italic'))

  ## prionace
    curve.prio <-  ggplot(data = prio.standards, aes(log(Quantity),Cq)) + 
      geom_point()+
      xlab('Quantity (copies/µL, log)')+
      geom_abline()+
      xlim(1.5,12)+
      theme_bw()+
      ggtitle('P. glauca')+
      theme(plot.title = element_text(size=10, face='italic'))
  
  ## alopias
    curve.alop <-  ggplot(data = alop.standards, aes(log(Quantity),Cq)) + 
      geom_point()+
      xlab('Quantity (copies/µL, log)')+
      geom_abline()+
      xlim(1.5,12)+
      theme_bw()+
      ggtitle('A. vulpinas')+
      theme(plot.title = element_text(size=10, face='italic'))
  
  ## lamna 
    curve.lamn <-  ggplot(data = lamn.standards, aes(log(Quantity),Cq)) + 
      geom_point()+
      xlab('Quantity (copies/µL, log)')+
      geom_abline()+
      xlim(1.5,12)+
      theme_bw()+
      ggtitle('L. nasus')+
      theme(plot.title = element_text(size=10, face='italic'))
  
  ## scombrus
    curve.scro <-  ggplot(data = scro.standards, aes(log(Quantity),Cq)) + 
      geom_point()+
      xlab('Quantity (copies/µL, log)')+
      geom_abline()+
      xlim(1.5,12)+
      theme_bw()+
      ggtitle('S. scombrus')+
      theme(plot.title = element_text(size=10, face='italic'))
  
  ## sprattus     
    curve.spra <-  ggplot(data = spra.standards, aes(log(Quantity),Cq)) + 
      geom_point()+
      xlab('Quantity (copies/µL, log)')+
      geom_abline()+
      xlim(1.5,12)+
      theme_bw()+
      ggtitle('S. sprattus')+
      theme(plot.title = element_text(size=10, face='italic'))
  

  ## stack/panel them 
      curves.plots <- curve.engr+curve.prio+curve.alop+curve.lamn+curve.scro+curve.spra
      ggsave(curves.plots, file='data_edna/figures_and_tables/comparingmethods/StandardsCurves_allSpecies_2023.png',device=png,units='in',height=11,width=4.5,dpi=900)

#### Cq.value per Technical replicate and Sampling Event - MP and WB. These plots do NOT include where field repsamples NOR field controls had NAs. species specific 
    #engraulis
    stat_sum_engraulis_methodtypes_byEVENT <- ggplot(data = techreps.repsamples)+
      stat_summary(aes(x = eventID, y = engraulis.Cq, color = methodtype), fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median, position = position_dodge(0.6))+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      ylab('Cq')+
      xlab('Sampling Event')+
      ylim(c(32,42))+
      theme_bw()+
      ggtitle('E. encrasicolus')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size=10, face='italic'))+
      labs(color="Method")
    
    # alopias
    stat_sum_alopias_methodtypes_byEVENT <- ggplot(data = repsamples)+
      stat_summary(aes(x = eventID, y = alopias.cq.mean, color = methodtype), fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median, position = position_dodge(0.6))+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      ylab('Cq (median, IQR)')+
      xlab('Sampling Event')+
      ylim(c(32,42))+
      theme_bw()+
      ggtitle('A. vulpinas')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size=10, face='italic'))+
      labs(color="Method")         
    
    # prionace
    stat_sum_prionace_methodtypes_byEVENT <- ggplot(data = repsamples)+
      stat_summary(aes(x = eventID, y = prionace.cq.mean, color = methodtype), fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median, position = position_dodge(0.6))+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      ylab('Cq (median, IQR)')+
      xlab('Sampling Event')+
      ylim(c(32,42))+
      theme_bw()+
      ggtitle('P. glauca')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size=10, face='italic'))+
      labs(color="Method")         
    
    # lamna
    stat_sum_lamna_methodtypes_byEVENT <- ggplot(data = repsamples)+
      stat_summary(aes(x = eventID, y = lamna.cq.mean, color = methodtype), fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median, position = position_dodge(0.6))+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      ylab('Cq (median, IQR)')+
      xlab('Sampling Event')+
      ylim(c(32,42))+
      theme_bw()+
      ggtitle('L. nasus')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size=10, face='italic'))+
      labs(color="Method")         
    
    # scombrus
    stat_sum_scombrus_methodtypes_byEVENT <- ggplot(data = repsamples)+
      stat_summary(aes(x = eventID, y = scombrus.cq.mean, color = methodtype), fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median, position = position_dodge(0.6))+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      ylab('Cq (median, IQR)')+
      xlab('Sampling Event')+
      ylim(c(32,42))+
      theme_bw()+
      ggtitle('S. scombrus')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size=10, face='italic'))+
      labs(color="Method")         
    
    # sprattus 
    stat_sum_sprattus_methodtypes_byEVENT <- ggplot(data = repsamples)+
      stat_summary(aes(x = eventID, y = sprattus.cq.mean, color = methodtype), fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median, position = position_dodge(0.6))+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      ylab('Cq (median, IQR)')+
      xlab('Sampling Event')+
      ylim(c(32,42))+
      theme_bw()+
      ggtitle('S. sprattus')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size=10, face='italic'))+
      labs(color="Method")

#### Copies over time - continuous, one line for MP, one line for WBT. Replicate repsamples, not technical
  ## geom_smooth uses conditional means. standard error bars are calculated using predict(). - all species, fish and shark taxa plots go with the GLMER model table

  ## all species together
     copies_time_allspp <- ggplot(by_taxa, aes(sampDATE,allspp.copies.copies, group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies')+
      ggtitle('All Species')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
         
  ## fish species together 
    copies_time_fish <- ggplot(by_taxa, aes(sampDATE,fish.copies, group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies')+
      ggtitle('Fish Species')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")

  ## shark species together 
    copies_time_shark <- ggplot(by_taxa, aes(sampDATE,shark.copies, group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies')+
      ggtitle('Shark Species')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")

  ## species, one at a time - for supplementary
    # engraulis
    copies_time_engraulis <- ggplot(repsamples, aes(sampDATE,log(mean.engraulis.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies (log)')+
      ggtitle('E. encrasicolus')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # alopias
    copies_time_alopias <- ggplot(repsamples, aes(sampDATE,log(alopias.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies (log)')+
      ggtitle('A. vulpinas')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # prionace
    copies_time_prionace <- ggplot(repsamples, aes(sampDATE,log(prionace.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies (log)')+
      ggtitle('P. glauca')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # lamna
    copies_time_lamna <- ggplot(repsamples, aes(sampDATE,log(lamna.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies (log)')+
      ggtitle('L. nasus')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # scombrus
    copies_time_scombrus <- ggplot(repsamples, aes(sampDATE,log(scombrus.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies (log)')+
      ggtitle('S. scombrus')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # sprattus 
    copies_time_sprattus <- ggplot(repsamples, aes(sampDATE,log(sprattus.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies (log)')+
      ggtitle('S. sprattus')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")



#### stacked/panelled plots - qPCR results
## fish species together 
  fish_methods_comp_together <- (stat_sum_engraulis_methodtypes_byEVENT + copies_time_engraulis) / (stat_sum_scombrus_methodtypes_byEVENT + copies_time_scombrus) / (stat_sum_sprattus_methodtypes_byEVENT + copies_time_sprattus)  + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

## shark species together 
  shark_methods_comp_together <- (stat_sum_alopias_methodtypes_byEVENT + copies_time_alopias) / (stat_sum_prionace_methodtypes_byEVENT + copies_time_prionace) / (stat_sum_lamna_methodtypes_byEVENT + copies_time_lamna)  + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')
  
## by species, one at a time
  # engraulis
  engraulis_methods_comp_together <- stat_sum_engraulis_methodtypes_byEVENT / copies_time_engraulis  + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

  # alopias
  alopias_methods_comp_together <- stat_sum_alopias_methodtypes_byEVENT / copies_time_alopias + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

  # prionace
  prionace_methods_comp_together <- stat_sum_prionace_methodtypes_byEVENT / copies_time_prionace  + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

  # lamna
  lamna_methods_comp_together <- stat_sum_lamna_methodtypes_byEVENT / copies_time_lamna + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

  # scombrus
  scombrus_methods_comp_together <- stat_sum_scombrus_methodtypes_byEVENT / copies_time_scombrus + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

  # sprattus 
  sprattus_methods_comp_together <- stat_sum_sprattus_methodtypes_byEVENT / copies_time_sprattus  + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')




  ## SAVE - local R 
    ## individual plots
    ggsave(stat_sum_engraulis_methodtypes_byEVENT,file='data_edna/figures_and_tables/comparingmethods/engraulis_cqTechReps_bymethod_atSamplingEvents_2023_statsumplot.png',device=png,units='in',height=5,width=6.5,dpi=700)

    ggsave(copies_time_engraulis,file='data_edna/figures_and_tables/comparingmethods/engraulis_copies_bymethod_overtime_2023_geomsmoothplot.png',device=png,units='in',height=5,width=6.5,dpi=700)
    
    ggsave(dnacont_overtime,file='data_edna/figures_and_tables/comparingmethods/concentration_bymethod_overtime_2023_geomsmoothplot.png',device=png,units='in',height=5,width=6.5,dpi=700)

    ## stacked/panelled plots 
    ggsave(engraulis_methods_comp_together,file='data_edna/figures_and_tables/comparingmethods/engraulis_2023_stacked_CQperEvent_CopiesOverTime.png',device=png,units='in',height=8,width=4.5,dpi=900)


##########
## - Statistical tests
##########

  ## glmer with random effect on event ID 
    ## all species
      lmdata1 <- repsamples %>% 
          rowwise()%>%
          mutate(copies = sum(c_across(ends_with('copies')))) # adds up columns that have copy numbers of each species 
      lm1 <- glmer(copies ~ methodtype + (1|eventID), data = lmdata1, family = Gamma(link = 'log'))

    ## fish only
      lmdata.fish <- repsamples %>% 
          dplyr::select(eventID,methodtype,sampleID,replicateID,engraulis.copies,scrombrus.copies,sprattus.copies)%>%
          rowwise()%>%
          mutate(copies = sum(c_across(ends_with('copies')))) # adds up columns that have copy numbers of each species 
      lm2 <- glmer(copies ~ methodtype + (1|eventID), data = lmdata.fish, family = Gamma(link = 'log'))

    ## sharks only
      lmdata.sharks <- repsamples %>% 
          dplyr::select(eventID,methodtype,sampleID,replicateID,alopias.copies,prionace.copies,lamna.copies)%>%
          rowwise()%>%
          mutate(copies = sum(c_across(ends_with('copies')))) # adds up columns that have copy numbers of each species 
      lm3 <- glmer(copies ~ methodtype + (1|eventID), data = lmdata.sharks, family = Gamma(link = 'log'))

    ## code for engraulis only right  now
      lmdata.engr <- repsamples %>% 
          rowwise()%>%
          mutate(copies = sum(c_across(ends_with('copies')))+0.001) # add a constant because the Gamma can't have 0s. 
      lm.engr <- glmer(copies ~ methodtype + (1|eventID), data = lmdata.engr, family = Gamma(link = 'log'))

    ## save model RDS 
      saveRDS(lm1, 'data_edna/modelRDS/glmer_RandomEventID_allspecies_copies_byMethodType.RDS')
      saveRDS(lm2, 'data_edna/modelRDS/glmer_RandomEventID_fish_copies_byMethodType.RDS')
      saveRDS(lm3, 'data_edna/modelRDS/glmer_RandomEventID_sharks_copies_byMethodType.RDS')

  ## Deprecated by the GLMER - Calculate the mean SD between replicates within sampling events of the two methods across all species (as one value), then for fosh species lumped, and then for shark species lumped. Add in a t-test of each column in flextable. 
      sd.allspp <- repsamples%>%
        group_by(sampleID)%>%
        summarise(sd.engr.copies = sd(engraulis.copies, na.rm=TRUE),
          sd.alop.copies = sd(alopias.copies, na.rm=TRUE),
          sd.prio.copies = sd(prionace.copies, na.rm=TRUE),
          sd.lamn.copies = sd(lamna.copies, na.rm=TRUE),
          sd.scro.copies = sd(scrombrus.copies, na.rm=TRUE),
          sd.spra.copies = sd(sprattus.copies, na.rm=TRUE))
      
      sd.allspp.spread <- left_join(sd.allspp, repsamples%>%dplyr::select(sampleID, methodtype))%>%
        group_by(methodtype)%>%
        summarise(method.meanSD.engr.copies = mean(sd.engr.copies),
          method.meanSD.alop.copies = sd(sd.alop.copies, na.rm=TRUE),
          method.meanSD.prio.copies = sd(sd.prio.copies, na.rm=TRUE),
          method.meanSD.lamn.copies = sd(sd.lamn.copies, na.rm=TRUE),
          method.meanSD.scro.copies = sd(sd.scro.copies, na.rm=TRUE),
          method.meanSD.spra.copies = sd(sd.spra.copies, na.rm=TRUE))
     
      sd.allspp.gathered <- sd.allspp.spread%>%
         rowwise()%>%
         mutate(meanSD.allspp.btwrep = mean(c_across(c(method.meanSD.engr.copies,method.meanSD.alop.copies,method.meanSD.prio.copies,method.meanSD.lamn.copies,method.meanSD.scro.copies,method.meanSD.spra.copies)), na.rm = TRUE),
            meanSD.FISH.btwrep = mean(c_across(c(method.meanSD.engr.copies,method.meanSD.scro.copies,method.meanSD.spra.copies)), na.rm = TRUE),
            meanSD.SHARK.btwrep = mean(c_across(c(method.meanSD.alop.copies,method.meanSD.prio.copies,method.meanSD.lamn.copies)), na.rm = TRUE))%>%
         dplyr::select(methodtype,meanSD.allspp.btwrep,meanSD.FISH.btwrep,meanSD.SHARK.btwrep)%>%
          flextable()%>%
          set_header_labels(methodtype = 'Method', n = 'Sample Size', meanSD.allspp.btwrep = 'All Species\n\ Replicate Avg. SD', meanSD.FISH.btwrep = 'Fish Species\n\ Replicate Avg. SD', meanSD.SHARK.btwrep = 'Shark Species\n\ Replicate Avg. SD')%>%
          theme_zebra()%>%
          align(align = 'center', part = 'all')%>%
          font(fontname = 'Arial', part = 'all')%>%
          fontsize(size = 10, part = 'all')%>%
          autofit()

    ## repeat, listing out all species - for Supplementary materials
      sd.allspp2 <- left_join(sd.allspp, repsamples%>%filter()%>%dplyr::select(sampleID, methodtype))%>%
        group_by(methodtype)%>%
        summarise(method.meanSD.engr.copies = mean(sd.engr.copies),
          method.meanSD.scro.copies = mean(sd.scro.copies, na.rm=TRUE),
          method.meanSD.spra.copies = mean(sd.spra.copies, na.rm=TRUE),
          method.meanSD.alop.copies = mean(sd.alop.copies, na.rm=TRUE),
          method.meanSD.prio.copies = mean(sd.prio.copies, na.rm=TRUE),
          method.meanSD.lamn.copies = mean(sd.lamn.copies, na.rm=TRUE))%>%
          flextable()%>%
            set_header_labels(methodtype = 'Method', n = 'Sample Size', method.meanSD.engr.copies = 'E. encrasicolus\n\ Replicate Avg. SD', method.meanSD.scro.copies = 'S. scombrus\n\ Replicate Avg. SD', method.meanSD.spra.copies = 'S. sprattus\n\ Replicate Avg. SD', method.meanSD.alop.copies = 'A. vulpinas\n\ Replicate Avg. SD', method.meanSD.prio.copies = 'P. glauca\n\ Replicate Avg. SD', method.meanSD.lamn.copies = 'L. nasus\n\ Replicate Avg. SD')%>%
            theme_zebra()%>%
            align(align = 'center', part = 'all')%>%
            font(fontname = 'Arial', part = 'all')%>%
            fontsize(size = 10, part = 'all')%>%
            autofit()


    ## For now you can calculate it for engraulis only
      sd.allspp <- repsamples%>%
        group_by(sampleID)%>%
        summarise(sd.engr.copies = sd(engraulis.copies, na.rm=TRUE))
      sd.allspp2 <- left_join(sd.allspp, repsamples%>%dplyr::select(sampleID, methodtype))%>%
        group_by(methodtype)%>%
        summarise(method.meanSD.engr.copies = mean(sd.engr.copies))
      sd.allspp2  

      ## test if this will work with multiple columns (ie multiple species copy numbers) - It does work. Leave this here for later workshopping.
        temp <- repsamples %>% mutate(spB.copies = rnorm(n=nrow(repsamples),mean=.25, sd=6))
        sd.temp <- temp%>%
        group_by(sampleID)%>%
        summarise(sd.engr.copies = sd(engraulis.copies, na.rm=TRUE), sd.spB.copies = sd(spB.copies))
        temp.sd2 <- left_join(sd.temp, temp%>%dplyr::select(sampleID, eventID, methodtype))%>%
           group_by(methodtype)%>%
          summarise(method.meanSD.engr.copies = mean(sd.engr.copies),method.meansd.spB.copies = mean(sd.spB.copies)) 
        temp.sd2  

        temp.eventID.method.meanSD <- left_join(sd.temp, temp%>%dplyr::select(sampleID, eventID, methodtype)) %>% group_by(eventID, methodtype)%>%
          summarise(method.meanSD.engr.copies = mean(sd.engr.copies),method.meansd.spB.copies = mean(sd.spB.copies))

        temp.sd.allspp.gathered <- temp.sd2 %>%
         rowwise()%>%
         mutate(method.meanSD.allspp.btwrep = mean(c_across(c('method.meanSD.engr.copies','method.meansd.spB.copies')), na.rm = TRUE)) ## this adds a column with the information. the new column is the mean of the row. This would also allow the fish and shark data to remain for
        
         stat.test <- t.test(method.meanSD.allspp.btwrep ~ methodtype, data = temp.sd.allspp.gathered)


##########
## - Tables
##########

#### Standards Curve Test results 
  curvelineartestresult.flex <- curves %>%
  dplyr::select(-species, -dateofassay)%>%
  mutate(sp = case_when(sp == 'engraulis' ~ 'E. encrasicolus', sp == 'prionace' ~ 'P. glauca', .default = as.character(sp)))%>%
    flextable()%>%
      set_header_labels(sp = 'Species', intercept = 'Intercept', Slope = 'Slope', r.squared = 'R-squared', efficiency = 'Efficiency')%>%
      theme_zebra()%>%
      align(align = 'center', part = 'all')%>%
      font(fontname = 'Arial', part = 'all')%>%
      fontsize(size = 10, part = 'all')%>%
      autofit()

      save_as_image(curvelineartestresult.flex, 'standardsCurveTest_Table_2023_sharksAndFish.png', webshot = 'webshot2')

      
#### CQ.mean values for all replicates, field controls and lab controls, of both sample types - Supplementary. 15/2/2024: want to add a ratio (e.g. 2/4) describing how many technical reps amplified. But do then only use samples (replicates with mean copy number) for the table
  ## select down to vars for each dataset. rename to get things to match up 
    repsamples2 <- repsamples %>% dplyr::select('Assay.Role', 'date', 'Sample.Name', 'eventID', 'replicateID', 'methodtype', ends_with('.copies'), 'testID')
    repsamples2
    tech.samples2 <- techreps.samples %>% dplyr::select('Assay.Role', 'date', 'Sample.Name', 'replicateID', 'eventID', 'methodtype', ends_with('.copies'), 'testID')
    tech.samples2

    field.controls2 <- field.controls %>% dplyr::select('Assay.Role', 'date', 'Sample.Name', 'replicateID', 'eventID', 'methodtype', ends_with('.copies'), 'testID')
    techreps.field.controls2 <- techreps.field.controls %>% dplyr::select('Assay.Role', 'date', 'Sample.Name', 'eventID', 'replicateID', 'methodtype', ends_with('.copies'), 'testID')
    techreps.field.controls2

    ntc2 <- ntc %>%
      mutate(eventID = 'n/a', methodtype = 'n/a', date = 'n/a', Sample.Name = case_when(is.na(Sample.Name) ~ 'NTC_NA', .default = as.character(Sample.Name))) %>% 
      dplyr::select('Assay.Role', 'date', 'Sample.Name', 'eventID', 'methodtype',  ends_with('.copies'), 'testID', Cq)
    ntc2

  
  ## calculate the number of tech reps that amplified 
    tt <- techreps.samples %>% group_by(Sample.Name) %>% 
        tally(engraulis.Cq != 'NA')
    repsamples3 <- left_join(repsamples2, tt, relationship = 'one-to-one', by='Sample.Name')%>%
      mutate(engraulis.Amplified = paste0(n,'/3'))%>%
      mutate_if(is.numeric, round, digits = 3)
    repsamples3
    ttt <- techreps.field.controls %>% group_by(Sample.Name) %>% 
        tally(engraulis.Cq != 'NA')
    field.controls3 <- left_join(field.controls2, ttt, relationship = 'one-to-one', by='Sample.Name')%>%
      mutate(engraulis.Amplified = paste0(n,'/3'))%>%
      mutate_if(is.numeric, round, digits = 3)
    field.controls3
    tttt <- ntc2 %>% group_by(testID, Sample.Name) %>% 
        tally(Cq != 'NA')
    ntc3 <- left_join(ntc2, tttt, relationship = 'one-to-one', by=c('Sample.Name','testID'))%>%
      mutate(engraulis.Amplified = paste0(n,'/1'))%>%
      dplyr::select(-Cq)%>%
      rename(mean.engraulis.copies = engraulis.copies)%>% 
      mutate_if(is.numeric, round, digits = 3)
    ntc3


  ## make table - supplementary
   results.byrole.replicates <- bind_rows(repsamples3, field.controls3, ntc3) %>%
      dplyr::select(-n, -replicateID)%>%
      flextable()%>%
     set_header_labels(Sample.Name = 'Sample ID', Assay.Role = 'Assay Role', 'eventID' = 'Event ID', 'methodtype' = 'Method Type', testID = 'qPCR Test','date' = 'Date', mean.engraulis.copies = 'E. encrasicolus\n\ copies', engraulis.Amplified = 'Amplified\n\ (E. encrasicolus)')%>%
     theme_zebra()%>%
     align(align = 'center', part = 'all')%>%
     font(fontname = 'Arial', part = 'all')%>%
     fontsize(size = 10, part = 'all')%>%
     autofit()
   results.byrole.replicates

   save_as_image(results.byrole.replicates, 'data_edna/figures_and_tables/results_allsamples_triplicateCQmeans_2023eDNA.png', webshot = 'webshot2')
  
  ## make table shorter - mean copy number PER SAMPLE, not replicate, and amplified ratios  

    samples <- read_csv('data_edna/compiledspeciesqPCR_SamplingEventCompiledCopiesandCq_withmetadata_2023eDNACornishBlue.csv')
    
    ttttt <- techreps.samples %>% group_by(sampleID) %>% 
        tally(engraulis.Cq != 'NA')
    samples2 <- left_join(samples, ttttt, relationship = 'one-to-one', by='sampleID')%>%
      mutate(engraulis.Amplified = paste0(n,'/9'))
    summary(samples2)

    # update the field controls and ntcs 
      samples3 <- samples2 %>% 
        dplyr::select('date', 'eventID', 'methodtype',  ends_with('.copies'), 'testID', ends_with('Amplified'))%>%
        mutate(Assay.Role = 'Field.Sample')%>%
        relocate(Assay.Role,eventID,testID,methodtype, date,  ends_with('Cq'), ends_with('.copies'), ends_with('Amplified'))%>%
        mutate_if(is.numeric, round, digits = 3)
      ntc4 <- ntc3 %>%
        dplyr::select(-n,-Sample.Name) %>%
        rename(sample.mean.engraulis.copies = mean.engraulis.copies)%>%
        relocate(Assay.Role, eventID,testID,methodtype, date,  ends_with('Cq'), ends_with('.copies'), ends_with('Amplified'))
      field.controls4 <- field.controls3 %>% 
        dplyr::select(-n, -replicateID, -Sample.Name)%>%
        rename(sample.mean.engraulis.copies = mean.engraulis.copies)%>%
        relocate(Assay.Role, eventID, testID, methodtype, date,  ends_with('Cq'), ends_with('.copies'), ends_with('Amplified'))

      #### UPDATE WITH MORE SPECIES 
    results.byrole <- bind_rows(samples3, 
      field.controls4, ntc4)%>%
      flextable()%>%
       set_header_labels(Assay.Role = 'Assay Role', 'eventID' = 'Event ID', 'methodtype' = 'Method Type', testID = 'qPCR Test','date' = 'Date', sample.mean.engraulis.copies = 'E. encrasicolus\n\ copies (Sample mean)', engraulis.Amplified = 'Amplified\n\ (E. encrasicolus)')%>%
       theme_zebra()%>%
       align(align = 'center', part = 'all')%>%
       font(fontname = 'Arial', part = 'all')%>%
       fontsize(size = 10, part = 'all')%>%
       autofit()

    results.byrole

   save_as_image(results.byrole, 'data_edna/figures_and_tables/qPCRresults_allspecies_copiesPERevent_AmplifiedRatio_2023eDNA.png', webshot = 'webshot2')









#### sample size and lowest/median/highest copy number per taxa
size <- samples %>%
  group_by(methodtype)%>%
  distinct(eventID)%>%
  tally()
size

by.taxa

#### Soaktime - Supplementary
  soaks <- samples %>%
    filter(methodtype=='metaprobe')%>%
    mutate(timeOUT = case_when(timeOUT == 'TBD' ~ as.character(timeIN), is.na(timeOUT) ~ as.character(timeIN), .default = as.character(paste0(timeOUT,':00'))))%>%
  mutate(timeOUT2 = hms(timeOUT), timeIN2 = hms(timeIN), soaktime = timeOUT2-timeIN2)%>%
    dplyr::select(Sample.Name, replicateID, eventID, soaktime)%>%
    mutate(soaktime = seconds(soaktime))
 soaks.flex <-  soaks %>% flextable() %>%
    set_header_labels(Sample.Name = 'Sample ID', replicateID = 'replicateID', 'eventID' = 'Event ID', soaktime = 'Soak Time (sec)')%>%
    theme_zebra()%>%
    align(align = 'center', part = 'all')%>%
    font(fontname = 'Arial', part = 'all')%>%
    fontsize(size = 10, part = 'all')%>%
    autofit()
    
  save_as_image(soaks.flex, 'soaktimes_metaprobes_seconds_2023.png', webshot = 'webshot2')
  
  
#### Glmer model of copy number between methods type, random effect of event ID
lm1 <- readRDS('data_edna/modelRDS/glmer_RandomEventID_allspecies_copies_byMethodType.RDS') ## all species
lm2 <- readRDS('data_edna/modelRDS/glmer_RandomEventID_fish_copies_byMethodType.RDS') ## fish 
lm3 <- readRDS('data_edna/modelRDS/glmer_RandomEventID_sharks_copies_byMethodType.RDS') ## sharks   

models <- c(lm1, lm2, lm3)   

summary <- modelsummary(models, coef_rename = c('methodtypewaterbottle' = 'Water Bottle', 'methodmetaprobe' = 'Metaprobe'),fmt=3,estimate='estimate', statistic='conf.int',stars=FALSE,conf_level=0.95,output='flextable')

summary_flex <- summary %>%
    set_header_labels('Model 1' = 'All Species', 'Model 2' = 'Fish Species', 'Model 3' = 'Shark Species')%>%
    theme_zebra()%>%
    align(align = 'center', part = 'all')%>%
    font(fontname = 'Arial', part = 'all')%>%
    fontsize(size = 10, part = 'all')%>%
    autofit() # this table goes with the geom_smooth graphs of copies over time by taxon

save_as_image(summary_flex, 'summary_glmers_eventIDrandom_copies_byMethod.png', webshot = 'webshot2')

  













