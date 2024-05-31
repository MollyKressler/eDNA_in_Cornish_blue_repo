## Amplification plots from qPCR 

## Molly M Kressler :: May 2024


#######
## - Load workspace
#######

pacman::p_load(sf,dplyr,lubridate,readr,readxl, ggplot2)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/EDNA/data_edna/')

#######
## - Load data 
#######

	amp1 <- read_excel('qPCRresults/ESI_assays_spring2024/30052024-KRESSLER-STANDARDS-ASSAY1.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal')%>%
		mutate(testID = 'assay1')
	amp2 <- read_excel('qPCRresults/ESI_assays_spring2024/30052024-KRESSLER-STANDARDS-ASSAY2.xlsx', sheet = 'Amplification Data',range = 'A41:E3881', .name_repair = 'universal')%>%
		mutate(testID = 'assay2')
	setup1 <- read_excel('qPCRresults/ESI_assays_spring2024/30052024-KRESSLER-STANDARDS-ASSAY1.xlsx', sheet = 'Sample Setup',range = 'A41:L137', .name_repair = 'universal')%>%
		dplyr::select(-Sample.Name, -Sample.Color, -Biogroup.Name, -Biogroup.Color)
	setup2 <- read_excel('qPCRresults/ESI_assays_spring2024/30052024-KRESSLER-STANDARDS-ASSAY2.xlsx', sheet = 'Sample Setup',range = 'A41:L137', .name_repair = 'universal')%>%
		dplyr::select(-Sample.Name, -Sample.Color, -Biogroup.Name, -Biogroup.Color)

	data1 <- left_join(amp1, setup1%>%select(-Target.Name), relationship = 'many-to-one', by = 'Well')
	data1
	data2 <- left_join(amp2, setup2%>%select(-Target.Name), relationship = 'many-to-one', by = 'Well')
	data2

	data <- rbind(data1, data2) %>% filter(Rn != 'NA') 
	unique(data$Target.Name)

	write.csv(data, 'datafordan_ampplots.csv')

#######
## - Make plots
#######

	# engraulis #799ecb, scomber #003a78, prio #84cb7b, alo #477939, lamna #154200

	ggplot(data%>%filter(Target.Name == 'Lamna'), aes(x = Cycle, y = Rn, group = Quantity))+
		geom_smooth(span = 0.3, method = 'gam', se = FALSE)+
		scale_y_continuous(limits = c(-0.1,4))+
		theme_bw()

		col = c('Engraulis', 'Scomber', 'Prionace', 'Alopias', 'Lamna')
		values = c('#799ecb', '#003a78', '#84cb7b', '#477939', '#154200')

## not working



