
### Sampling locations around Cornwall

## created by Molly Kressler on 17 February 2023 

pacman::p_load(sf,tidyverse,ggplot2)
setwd('/Users/mollykressler/Documents/EDNA/data_edna')

# Need to define a bounding box around Cornwall and IOS - do this in Google. Import: 

box<-st_as_sf(st_read('startingboundarybox_forsampling.kml'))%>%st_transform(27700)


ggplot()+geom_sf(data=box,alpha=0.2,fill='goldenrod3')+geom_sf(data=st_union(cw),alpha=0.3,fill='cadetblue4')+geom_sf(data=st_union(kernios),alpha=0.5,fill='goldenrod3')+theme_bw()


# Cut out KernIos and coastal water shapes
	cw<-st_as_sf(st_read('kernow_ios_coastalwater.shp'))
	ios<-st_as_sf(st_read('islesofscilly.shp'))
	kern<-st_as_sf(st_read('kernow.shp'))
	kernios<-bind_rows(kern,ios)
	
	landsea<-st_as_sf(st_union(st_buffer(st_combine(kernios),700),st_combine(cw))) # adding the buffer fills in the cracks between the coastal water and the land, and any lakes within kernow. 

	box.cut<-st_difference(box,landsea)
	ggplot()+geom_sf(data=box.cut,alpha=0.2,fill='goldenrod3')+theme_bw() # good enough

	#st_write(box.cut,'startingboundarybox_forsampling_landRemoved.shp')

# spatially sample 15 random sites - what I should actually do is: divide the box.cut into squares/hexagons, and then randomly sample from within each hexagon. That way the points are more evenly distirbuted around the area. 

	boxhex<-st_make_grid(box.cut,n=4.5,square=FALSE)
	hex<-st_intersection(box.cut,boxhex) # cuts the hexagon grid to only within the area

	hex2<-hex%>%dplyr::select(-Name,-Description)%>%slice(-1)%>%mutate(hexID=paste0('h',c(1:27)),.before='geometry')
	hex2$hexID<-as.factor(hex2$hexID)
		# need to remove the points (only 26 hexagons in study area but showing 28 because of imperfections in land mass shapes)

		ggplot()+geom_sf(data=hex2)+geom_sf(data=hex2[20,],col='violetred4')+theme_bw()

	sites2<-as.data.frame(matrix(ncol=3,nrow=27))
		colnames(sites2)<-c('siteID','geo.x','geo.y')
		sites2$siteID<-paste0('h',c(1:27))
		sites2$geo.x<-1
		sites2$geo.y<-1
		sites3<-st_as_sf(sites2,coords=c('geo.x','geo.y'),crs=27700)		


	for(i in 1:27){
		i<-paste(i)
		h<-paste0('h',i)
		g<-st_geometry(hex2[i,][2])
		p1<-st_sample(hex2[i,],type='random',size=1)
		sites3$geometry[sites3$siteID==h]<-p1
			}
	ggplot()+geom_sf(data=hex2,alpha=0.2,fill='cadetblue2')+geom_sf(data=sites3,col='violetred4',alpha=0.3)+theme_bw() 

	sites3%>%print(n=Inf)

	st_write(hex2,'hexagon_grid_edna_sampling_kernios_nocoastal.shp')
	st_write(sites3,'sampling_sites_randomgen_onePERhexagon_kernios.shp')





# Specifically within the Fal region? Define a boundng box that stretches up to 20 miles from Fal Harbour entrance, cut out land and coastal water, and randomly sample two locations 





