
## Comparing Metaprobes to Water Bottles 
### using all species dataset 

## created by Molly M Kressler :: February 2024

##################
##########
## - Load at start 
##########
pacman::p_load(sf,tidyverse,lubridate,readr,readxl,lubridate,ggplot2,patchwork,cowplot,stringr,flextable)

## For local R, not for server RStudio

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/EDNA')

copies <- read_csv("data_edna/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv")
ntc <- read_csv("data_edna/resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv", 
                col_types = cols(...1 = col_skip(), X = col_skip()))%>%
  mutate(Assay.Role = case_when(Assay.Role == 'Negative' ~ 'NTC',
                                Assay.Role == 'Positive' ~ 'Standard'
                                , .default = as.character(Assay.Role)))
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

samples <- copies %>% filter(!grepl('WBC', replicateID)) %>%
  mutate('Assay.Role' = 'Field.Sample') # This is all field samples, water bottles and metaprobes, without the field controls 

##################

##########
## - Plots
##########

#### individual plots
###
  ## Cq.value/mean by Sampling Event - MP and WB.
  ## geom_stat_summary, with ranges and medians for WBTs and MPs, and a star geom_point for WBC. These plots do NOT include where field samples NOR field controls had NAs. 
    #engraulis
    stat_sum_engraulis_methodtypes_byEVENT <- ggplot(data = samples)+
      stat_summary(aes(x = eventID, y = engraulis.cq.mean, color = methodtype), fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median, position = position_dodge(0.6))+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      ylab('Cq (median, IQR)')+
      xlab('Sampling Event')+
      ylim(c(32,42))+
      theme_bw()+
      ggtitle('E. encrasicolus')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size=10, face='italic'))+
      labs(color="Method")
    
    # alopias
    stat_sum_alopias_methodtypes_byEVENT <- ggplot(data = samples)+
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
    stat_sum_prionace_methodtypes_byEVENT <- ggplot(data = samples)+
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
    stat_sum_lamna_methodtypes_byEVENT <- ggplot(data = samples)+
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
    stat_sum_scombrus_methodtypes_byEVENT <- ggplot(data = samples)+
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
    stat_sum_sprattus_methodtypes_byEVENT <- ggplot(data = samples)+
      stat_summary(aes(x = eventID, y = sprattus.cq.mean, color = methodtype), fun.min = function(z) { quantile(z,0.25) },fun.max = function(z) { quantile(z,0.75) },fun = median, position = position_dodge(0.6))+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      ylab('Cq (median, IQR)')+
      xlab('Sampling Event')+
      ylim(c(32,42))+
      theme_bw()+
      ggtitle('S. sprattus')+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size=10, face='italic'))+
      labs(color="Method")
  
  ###
  ## Copies over time - continuous, one line for MP, one line for WBT.
  ## geom_smooth uses conditional means. standard error bars are calculated using predict(). 
    # engraulis
    
    copies_time_engraulis <- ggplot(samples, aes(sampDATE,log(engraulis.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
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
    copies_time_alopias <- ggplot(samples, aes(sampDATE,log(alopias.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
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
    copies_time_prionace <- ggplot(samples, aes(sampDATE,log(prionace.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
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
    copies_time_lamna <- ggplot(samples, aes(sampDATE,log(lamna.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
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
    copies_time_scombrus <- ggplot(samples, aes(sampDATE,log(scombrus.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
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
    copies_time_sprattus <- ggplot(samples, aes(sampDATE,log(sprattus.copies+0.00001), group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Copies (log)')+
      ggtitle('S. sprattus')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")

  ###
  ## Concentration over time - continuous, one line for MP, one line for WBT
    # engraulis
    
    dnacont_overtime_engraulis <- ggplot(samples, aes(sampDATE,engraulis.dnaCont, group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Nanodrop Concentration')+
      ggtitle('E. encrasicolus')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # alopias
    dnacont_overtime_alopias <- ggplot(samples, aes(sampDATE,alopias.dnaCont, group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Nanodrop Concentration')+
      ggtitle('A. vulpinas')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # prionace
    dnacont_overtime_prionace <- ggplot(samples, aes(sampDATE,prionace.dnaCont, group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Nanodrop Concentration')+
      ggtitle('P. glauca')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # lamna
    dnacont_overtime_lamna <- ggplot(samples, aes(sampDATE, lamna.dnaCont, group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Nanodrop Concentration')+
      ggtitle('L. nasus')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # scombrus
    dnacont_overtime_scombrus <- ggplot(samples, aes(sampDATE,scombrus.dnaCont, group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Nanodrop Concentration')+
      ggtitle('S. scrombrus')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")
    
    # sprattus 
    dnacont_overtime_sprattus <- ggplot(samples, aes(sampDATE,sprattus.dnaCont, group = methodtype, col = methodtype, fill = methodtype))+
      geom_smooth(span = 3, se = TRUE, level = 0.95, alpha=0.3)+
      scale_color_manual(values=c('#DAA507','#8EC7D2'))+
      scale_fill_manual(values=c('#DAA507','#8EC7D2'))+
      xlab('Sampling Date (2023)')+ 
      ylab('Nanodrop Concentration')+
      ggtitle('S. sprattus')+
      theme_bw()+
      theme(plot.title = element_text(size=10, face='italic'))+
      labs(fill="Method", col= "Method")


#### stacked/panelled plots 
# engraulis
engraulis_methods_comp_together <- stat_sum_engraulis_methodtypes_byEVENT / copies_time_engraulis / dnacont_overtime_engraulis + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

# alopias
alopias_methods_comp_together <- stat_sum_alopias_methodtypes_byEVENT / copies_time_alopias / dnacont_overtime_alopias + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

# prionace
prionace_methods_comp_together <- stat_sum_prionace_methodtypes_byEVENT / copies_time_prionace / dnacont_overtime_prionace + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

# lamna
lamna_methods_comp_together <- stat_sum_lamna_methodtypes_byEVENT / copies_time_lamna / dnacont_overtime_lamna + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

# scombrus
scombrus_methods_comp_together <- stat_sum_scombrus_methodtypes_byEVENT / copies_time_scombrus / dnacont_overtime_scombrus + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

# sprattus 
sprattus_methods_comp_together <- stat_sum_sprattus_methodtypes_byEVENT / copies_time_sprattus / dnacont_overtime_sprattus + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')


## SAVE - local R 
## individual plots
ggsave(stat_sum_engraulis_methodtypes_byEVENT,file='data_edna/figures_and_tables/comparingmethods/engraulis_cqmean_bymethod_atSamplingEvents_2023_statsumplot.png',device=png,units='in',height=5,width=6.5,dpi=700)
ggsave(copies_time_engraulis,file='data_edna/figures_and_tables/comparingmethods/engraulis_copies_bymethod_overtime_2023_geomsmoothplot.png',device=png,units='in',height=5,width=6.5,dpi=700)
ggsave(dnacont_overtime_engraulis,file='data_edna/figures_and_tables/comparingmethods/engraulis_concentration_bymethod_overtime_2023_geomsmoothplot.png',device=png,units='in',height=5,width=6.5,dpi=700)

## stacked/panelled plots 
ggsave(engraulis_methods_comp_together,file='data_edna/figures_and_tables/comparingmethods/engraulis_2023_stacked_ConttOverTime_CQperEvent_CopiesOverTime.png',device=png,units='in',height=11,width=4.5,dpi=900)


##########
## - Statistical tests
##########
  samples
  ## Calculate the mean SD between replicates within sampling events of the two methods across all species
    ## a df with one value for MP and one for WB = the mean SD between replicate samples - across all species 
        ## will need to pivot the data long, so that there is one column for metaprobes         and one for water. bottles, and one column for replicate ID, and one column for         species assay. and then each row will be a copy number. Unique rows will be             identified with species AND replicate ID. 
        # select down to species.copies columns, replicate ID, and methodtype. 
        a <- samples %>% dplyr::select(engraulis.copies, Sample.Name, replicateID, methodtype)
        aa <- a %>% pivot_longer(values_to = 'methodtype', cols = 'engraulis.copies')
      
    ## repeat, for only all fish species
  
    ## repeat, for only all elasmobranch species 
  
    ## For now you can calculate it for engraulis only, without pivoting the data:  
      sd.allspp <- samples%>%
        group_by(sampleID)%>%
        summarise(sd.engr.copies = sd(engraulis.copies, na.rm=TRUE))
      sd.allspp2 <- left_join(sd.allspp, samples%>%dplyr::select(sampleID, methodtype))%>%
        group_by(methodtype)%>%
        summarise(method.meanSD.engr.copies = mean(sd.engr.copies))
      sd.allspp2  



  ##########
  ## - Tables
  ##########
  ## Broad descriptive statistics - one table for whole project, 2 rows, 6 columns: samples sizes of each (sampling events and replicates), average concentration, lowest copy number non-zero, highest copy number non-zero, median copy number non-zero
  # a long data frame, one row per replicate sample, columns for each of the above. then we'll create a summary df. 
  samples # Sample.Name = unique Replicate ID, methodtype = grouping factor, 
  
  ## not species specific. only sample sizes, and mean SD between replicates 
  size <- samples %>%
    group_by(methodtype)%>%
    distinct(eventID)%>%
    tally()
  size
  tbl <- left_join(size, sd.allspp2)
  tbl
  
    tbl.all <- tbl %>% mutate(method.meanSD.engr.copies = round(method.meanSD.engr.copies,3))%>% 
      flextable()%>%
      set_header_labels(methodtype = 'Method', n = 'Sample Size', method.meanSD.engr.copies = 'Replicate Avg. SD')%>%
      theme_zebra()%>%
      align(align = 'center', part = 'all')%>%
      font(fontname = 'Arial', part = 'all')%>%
      fontsize(size = 10, part = 'all')%>%
      autofit()
      
  
  ## species specific with copy number information 
  ## engraulis
  engr.size
  engraulis.comparison.summary.table <- 
    


