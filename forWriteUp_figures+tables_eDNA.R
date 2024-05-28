## - Figures and Tables for eDNA in Cornwall 


## created 29 September 2023 
## by Molly M Kressler

########
## Load data  
########
## cleaned in eDNA_datacleaning.R

pacman::p_load(sf,tidyverse,dplyr,ggplot2, patchwork, cowplot,lubridate,flextable, rnaturalearth)

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/')

## qPCR results, all species, not standards, test assay positive nor negative test controls 
sp1 <- st_as_sf(st_read('EDNA/data_edna/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.shp'),crs='WGS84')
sp2 <- read.csv('EDNA/data_edna/compiledspeciesqPCR_copiesresults_withmetadata_2023eDNACornishBlue.csv')
## qPCR standards, test assay positive and negative test controls 
c <- read.csv('EDNA/data_edna/resultsANDcopies_perStandard_andNegControl_En.encras_TESTS_1to7_NOV2023.csv')# quick fix of label problem
## land shapes
coastalwater <- st_as_sf(st_read('EDNA/data_edna/kernow_ios_coastalwater.shp'), crs='WGS84')
kernow <- st_as_sf(st_read('EDNA/data_edna/kernow.shp'), crs='WGS84')
ios <- st_as_sf(st_read('EDNA/data_edna/islesofscilly.shp'), crs='WGS84')
kernios <- st_as_sf(st_union(kernow, ios), crs='WGS84')
kernios <- st_transform(kernios, crs='WGS84')
## inset globe 
world <- ne_countries(scale='medium', returnclass = 'sf') %>% dplyr::select(sovereignt,region_un, region_wb, subregion,name_sort,abbrev,admin, geometry) %>% filter(sovereignt == 'United Kingdom')
world


#########
## - Locations 
#########

	location_extractions <- ggplot()+
		geom_sf(data=coastalwater,fill='cadetblue4', lwd=0.5, alpha = 0.2)+
		geom_sf(data=kernios,fill='grey82', lwd=0.5)+
		geom_sf(data=st_jitter(st_geometry(sp1), factor=0.01),size=2,pch=20, alpha = 0.5)+
		theme_bw() 
	location_extractions # this is at the sample replicate level
	
	uk <- ggplot()+
		geom_sf(data = world%>%filter(subregion=='Northern Europe'))+
		theme_void()+
		theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(), panel.border = element_rect(fill=NA, colour = 'black'))+
		geom_rect(aes(xmin=-6.5,xmax=-4.1,ymin=49.8,ymax=50.94),col='goldenrod2',fill='goldenrod2',alpha=0.4,lwd=.8)
	uk

	sample_location_extractions <- ggplot()+
		geom_sf(data=coastalwater,fill='cadetblue4', col=NA, lwd=0.5, alpha = 0.5)+
		geom_sf(data=kernios, lwd=0.5)+
		geom_sf(data=st_jitter(st_geometry(sp1%>%distinct(samplID, .keep_all = TRUE)), factor=0.005), col='#000A52',size=3,pch=20, alpha = 0.5)+
		theme_bw() 
	sample_location_extractions # this is at the sample replicate level

	map <- ggdraw()+
		draw_plot(sample_location_extractions)+
		draw_plot(uk, height=0.4, x = -0.29, y = 0.52)
	map

	ggsave(map,file = 'EDNA/data_edna/figures_and_tables/samplinglocations_2023_withUKinset.png',device='png',units='in',dpi=850,height=5,width=5)

	## Distance from shore 

	sp11 <- sp1 %>% mutate(dist2shore.km = st_distance(sp1, kernios)/1000)
	min(sp11$dist2shore.km) # 0.03 km
	median(sp11$dist2shore.km) # 6.149 km
	max(sp11$dist2shore.km) # 50.585 km





#########
## -  
#########






#########
## -  
#########





#########
## -  
#########






#########
## -  
#########







#########
## -  
#########
