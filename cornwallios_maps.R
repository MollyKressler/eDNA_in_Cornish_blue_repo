## Mapping Cornwall and the IOS
## start trying to map out areas for sampling
# initial shapefiles of coutnies and regions downloaded from the Office for National Statistics Open Geography portal. 

setwd('/Users/mollykressler/Documents/EDNA')

pacman::p_load(tidyverse,ggplot2,sf)

# land - from Office for National Statistics Open Geography portal 

ct<-st_as_sf(st_read('englandandwales_countiesUA_offNationalStats/CTYUA_Dec_2015_FCB_in_England_and_Wales.shp'),crs='WGS84')
regions<-st_as_sf(st_read('england_regions_offNationalstats/Regions_(December_2022)_EN_BFC.shp'),crs='WGS84')
ios<-ct%>%filter(ctyua15cd=='E06000053')
kern<-ct%>%filter(ctyua15cd=='E06000052')
sw<-regions%>%filter(RGN22CD=='E12000009')
	ggplot()+geom_sf(data=sw,fill='grey92', alpha=0.3)+geom_sf(data=kern)+geom_sf(data=ios)+theme_bw()

#st_write(ios,'islesofscilly.shp',driver='ESRI Shapefile')
#st_write(kern, 'kernow.shp',driver='ESRI Shapefile')

kernios<-bind_rows(kern,ios)

# coastal water - from DEFRA, pre (rough) cut to region

wet<-st_as_sf(st_read('coastal_waterboundary_DEFRA/WFD_Coastal_Water_Bodies_Cycle_1.shp'),crs='WGS84')
wetcut<-st_as_sf(st_intersection(st_buffer(kernios,13000),wet))

map<-ggplot()+geom_sf(data=wetcut,col='cadetblue2',fill='cadetblue2',alpha=0.7)+geom_sf(data=kernios)+theme_bw()+theme_bw()
ggsave(map,file='cornwallIOSmap_wthcoastalwater.png',device='png',units='in',height=6,width=6,dpi=800)

#st_write(wetcut,'kernow_ios_coastalwater.shp',driver='ESRI Shapefile')

## For Sci Comms - pngs of cornwall

png.cornwall<-ggplot()+geom_sf(data=st_union(kernios),fill='white')+theme_void()
ggsave(png.cornwall,file='png_kernios_whitefill.png',device='png',dpi=1200,unit='in',width=5)

