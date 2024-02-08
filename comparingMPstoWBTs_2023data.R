
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

#### individual plots
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


# prionace


# lamna


# scombrus


# sprattus 



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


# prionace


# lamna


# scombrus


# sprattus 



#### stacked/panelled plots 
# engraulis
engraulis_methods_comp_together <- stat_sum_engraulis_methodtypes_byEVENT / copies_time_engraulis / dnacont_overtime_engraulis + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', legend.justification = 'centre', legend.direction = 'vertical')

# alopias


# prionace


# lamna


# scombrus


# sprattus 



## SAVE - local R 
## individual plots
ggsave(stat_sum_engraulis_methodtypes_byEVENT,file='data_edna/figures_and_tables/comparingmethods/engraulis_cqmean_bymethod_atSamplingEvents_2023_statsumplot.png',device=png,units='in',height=5,width=6.5,dpi=700)
ggsave(copies_time_engraulis,file='data_edna/figures_and_tables/comparingmethods/engraulis_copies_bymethod_overtime_2023_geomsmoothplot.png',device=png,units='in',height=5,width=6.5,dpi=700)
ggsave(dnacont_overtime_engraulis,file='data_edna/figures_and_tables/comparingmethods/engraulis_concentration_bymethod_overtime_2023_geomsmoothplot.png',device=png,units='in',height=5,width=6.5,dpi=700)

## stacked/panelled plots 
ggsave(engraulis_methods_comp_together,file='data_edna/figures_and_tables/comparingmethods/engraulis_2023_stacked_ConttOverTime_CQperEvent_CopiesOverTime.png',device=png,units='in',height=11,width=4.5,dpi=900)


##########
## - Tables
##########
## Broad descriptive statistics - one table for whole project, 2 rows, 6 columns: samples sizes of each (sampling events and replicates), average concentration, lowest copy number non-zero, highest copy number non-zero, median copy number non-zero


##########
## - Statistical tests
##########








