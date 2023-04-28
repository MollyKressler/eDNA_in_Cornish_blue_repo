
### Sampling locations around Cornwall

## created by Molly Kressler on 17 February 2023 

pacman::p_load(sf,tidyverse,ggplot2,ggsn)
setwd('/Users/mollykressler/Documents/EDNA/data_edna')


############## Load at start 

	area<-st_as_sf(st_read('startingboundarybox_forsampling_landRemoved.shp'))
	cw<-st_as_sf(st_read('kernow_ios_coastalwater.shp'))
	ios<-st_as_sf(st_read('islesofscilly.shp'))
	kern<-st_as_sf(st_read('kernow.shp'))
	kernios<-bind_rows(kern,ios)

	locations<-st_as_sf(st_read('kernow_locations.kml'))%>%st_transform(27700)

	Falharbourlimit<-st_as_sf(st_read('FalHarbourShapefiles/harbour_limits_FHC.shp')) # Fal Harbour Shapefles 
		Falharbourlimit$id<-'FalH' 
		
		# majority is within the coastal waters. 


	hex<-st_as_sf(st_read('hexgrid_edna/hexagon_grid_edna_sampling_kernios_nocoastal.shp')) # proposed hexagon grid
		# update the shapefle as partners agree to areas and remove hexagons to out to sea
			 st_write(hex,'hexgrid_edna/hexagon_grid_edna_sampling_kernios_nocoastal.shp',append=FALSE)

	area_with_land<-ggplot()+geom_sf(data=area,alpha=0.2,fill='cadetblue4')+geom_sf(data=kernios,alpha=0.8,fill='grey82')+theme_bw() 
	hex_with_land<-ggplot()+geom_sf(data=hex,alpha=0.2,fill='cadetblue4')+geom_sf(data=kernios,alpha=0.8,fill='grey82')+theme_bw() 
		# ggsave(area_with_land,file='comms_figures/samplingarea_nogrd_withland.png',device='png',units='in',height=7,width=7,dpi=900)


############## Hexagon map with hex-ID labels
	# want a map with labels of hexagons, labels of locations, land shape, cw shape, and hexagon grid.

	hexlabelled<-ggplot()+geom_sf(data=st_union(kernios),alpha=0.5,fill='grey82')+geom_sf(data=locations,pch=19,size=2,col='black')+geom_sf(data=hex,alpha=0.05,fill='cadetblue4')+geom_sf_text(data=hex,aes(label=hexID),col='black',size=4,vjust=.5)+geom_sf_text(data=locations,aes(label=Name),col='black',hjust=-0.2,vjust=0.1,size=2.5)+north(data=hex,location='topleft',scale=0.12,symbol=12)+xlab(NULL)+ylab(NULL)+theme_bw() 
	ggsave(hexlabelled,file='comms_figures/labelled_hexagons_kernios_locations.png',device='png',units='in',height=7,width=7,dpi=1500)
 

############## Update hex shapefile with informaiton on which parnters have commited to which hexagons

	hex%>%print(n=Inf)
	# from emails from 2-3 March 2023:
		# Padstow Sea life safars  can confirm h20
		# Atlantic diving (chris lowe, newquay) can do Boscastle to st ives, offshore more difficult but does run shark cage diving in the summer. summer is irregular and would need to pay for fuel. or go as commercial crew - need qualifs. --- hexagon h16
		# Falmouth h19
		# marine discovery - penzance, H13 confirmed

	hex[20,2]<-as.character('Padstow, Sea Life Safaris')
	hex[16,2]<-as.character('Newquay, CLowe at Atlantic Diving')
	hex[19,2]<-as.character('FalHarbour, VSpooner')
	hex[13,2]<-as.character('Penzance, Marine Discovery')


	# simplify the hexagon map - remove the edge ones that are very far away and we cant sample - 2,5,11,17,,4,8,1,3,6,12,18,21,26

	hexs2remove<-c('h2','h5','h11','h17','h4','h8','h1','h3','h6','h12','h18','h21','h26')
	hex<-hex%>%filter(!(hexID%in%hexs2remove))
	# re-label the hexagons 1-n
	hex<-hex%>%mutate(hexID=paste0('h',c(1:14)),.before='geometry')

	hex[1,2]<-as.character('IOS, OExeter')

	
############## Working out hexagon areas for partners in regions of operation
#### Fal Harbour SV Killigrew

	# plot the limit(s) 

		ggplot()+geom_sf(data=cw,alpha=0.5,fill='goldenrod2')+geom_sf(data=Falharbourlimit,alpha=0.5,fill='violetred4')+geom_sf(data=hex,alpha=0.5,fill='cadetblue4')+theme_bw() 
	# identify the hexagon that Fal limits are within - can see from the map it is only contained within one. 

		st_intersects(hex,Falharbourlimit)%>%print(n=Inf) # Hexagon 19

	# Add a column in hex, and add notes about whos authority/capacty is there.

		hex<-hex%>%mutate(partners='tbd',.before='geometry')
		hex[19,2]<-'FalHarbour'
		hex%>%print(n=Inf)

	# Killigrew can go out 20 nautical miles from the coast. 

		falhex<-hex[19,]
		falmax<-st_buffer(st_centroid(Falharbourlimit),32000)

		ggplot()+geom_sf(data=falmax,alpha=0.5,fill='goldenrod2')+geom_sf(data=Falharbourlimit,alpha=0.5,fill='violetred4')+geom_sf(data=hex[19,],alpha=0.5,fill='cadetblue4')+theme_bw() 
		# Killigrew can reach most of the hexagon

		#st_write(falhex,'hexgrid_edna/FalmouthHarbour_hexgrid_from_hexagon_grid_edna_sampling_kernios_nocoastal.shp')

	# make a plot to send to Vicki. 

		forfhc_land<-st_intersection(falmax,kernios)
		forFHC<-ggplot()+geom_sf(data=forfhc_land,alpha=0.5,fill='goldenrod2')+geom_sf(data=Falharbourlimit,alpha=0.5,fill='violetred4')+geom_sf(data=hex[19,],alpha=0.5,fill='cadetblue4')+north(data=forfhc_land,location='topleft',scale=0.11,symbol=12)+theme_bw() 
		ggsave(forFHC,file='comms_figures/samplingarea_offFHC.png',device='png',units='in',height=5,width=4,dpi=900)
 
#### Collin Trundle and IFCA, also out of Fal. Offered to help with hexagon east of FalHarbour Hex #19. - Hexagon # 22

	st_touches(hex,hex[19,])%>%print(n=Inf) # Hexagons:13,15,21,22,24
		
	hex[22,2]<-'CTrundleIFCA'
	hex%>%print(n=Inf)
	forCTrundle_land<-st_intersection(st_buffer(hex[22,],20000),kernios)
	#st_write(hex[22,],'hexgrid_edna/CTrundleIFCA_hexgrid_from_hexagon_grid_edna_sampling_kernios_nocoastal.shp')

	forCTatIFCA<-ggplot()+geom_sf(data=forCTrundle_land,alpha=0.5,fill='goldenrod2')+geom_sf(data=hex[22,],alpha=0.5,fill='cadetblue4')+north(data=forCTrundle_land,location='topleft',scale=0.11,symbol=12)+theme_bw() 

 	ggsave(forCTatIFCA,file='comms_figures/samplingarea_CTrundle_IFCA.png',device='png',units='in',height=5,width=7,dpi=900)


#### Penzance - Hexagon 13 (and 10 and 15 at a reach) 

 	a<-st_nearest_feature(locations[1,],hex) #13 
	hex[13,2]<-'Penzance'# change these names later if/when partners come on board
	hex[10,2]<-'Penzance adjacent'
	hex[15,2]<-'Penzance adjacent'
	hex%>%print(n=Inf)
	penzance_hexagons<-ggplot()+geom_sf(data=st_union(kernios),alpha=0.3,fill='grey82')+geom_sf(data=locations[1,],pch=19,size=2,col='black')+geom_sf(data=hex[15,],alpha=0.5,fill='cadetblue4')+geom_sf(data=hex[13,],alpha=0.5,fill='cadetblue4')+geom_sf(data=hex[10,],alpha=0.5,fill='cadetblue4')+geom_sf_text(data=hex[c(15,13,10),],aes(label=c('A','B','C')),col='black',size=6)+geom_sf_text(data=locations[1,],aes(label=Name),col='black',hjust=-0.5)+north(data=kernios,location='topleft',scale=0.12,symbol=12)+xlab(NULL)+ylab(NULL)+theme_bw() 
	# 13 is the direct hexagon
	# To the north and south there are two more which some tours might operate within. include? hex 10 and 
 
 	ggsave(penzance_hexagons,file='comms_figures/samplingarea_Penzance.png',device='png',units='in',height=5,width=7,dpi=900)


##### Padstow & Newquay 

 	st_nearest_feature(hex[c(16,20),],hex) #16 and 20
	hex[20,2]<-'Padstow' # change these names later if/when partners come on board
	hex[16,2]<-'Newquay'
	hex%>%print(n=Inf)

	padstow_and_newquay_hexagons<-ggplot()+geom_sf(data=st_union(kernios),alpha=0.3,fill='grey82')+geom_sf(data=locations[c(2:3),],pch=19,size=2,col='black')+geom_sf(data=locations[c(2:3),],pch=19,size=2,col='black')+geom_sf(data=hex[20,],alpha=0.5,fill='cadetblue4')+geom_sf(data=hex[16,],alpha=0.5,fill='cadetblue4')+geom_sf(data=hex[17,],alpha=0.5,fill='cadetblue4')+geom_sf(data=hex[25,],alpha=0.5,fill='cadetblue4')+geom_sf(data=hex[14,],alpha=0.5,fill='cadetblue4')+geom_sf_text(data=hex[c(16,14,17,20,25),],aes(label=c('A','B','C','D','E')),col='black',size=6)+geom_sf_text(data=locations[c(2:3),],aes(label=Name),col='black',hjust=-0.15,vjust=0.01)+north(data=kernios,location='topleft',scale=0.12,symbol=12)+xlab(NULL)+ylab(NULL)+theme_bw() 

 	ggsave(padstow_and_newquay_hexagons,file='comms_figures/samplingarea_padstow_and_newquay_tentative.png',device='png',units='in',height=5,width=7,dpi=900)










############## MAKING SHAPEFILES, do not do every open. 
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
		ggplot()+geom_sf(data=landsea,alpha=0.2,fill='goldenrod3')+theme_bw() # good enough

		#st_write(box.cut,'startingboundarybox_forsampling_landRemoved.shp')

	# spatially sample 15 random sites - what I should actually do is: divide the box.cut into squares/hexagons, and then randomly sample from within each hexagon. That way the points are more evenly distirbuted around the area. 

		boxhex<-st_make_grid(box.cut,n=4.5,square=FALSE)
		hex<-st_intersection(box.cut,boxhex) # cuts the hexagon grid to only within the area

		hex<-hex%>%dplyr::select(-Name,-Description)%>%slice(-1)%>%mutate(hexID=paste0('h',c(1:27)),.before='geometry')
		hex$hexID<-as.factor(hex$hexID)
			# need to remove the points (only 26 hexagons in study area but showing 28 because of imperfections in land mass shapes)

			ggplot()+geom_sf(data=hex)+geom_sf(data=hex[20,],col='violetred4')+theme_bw()

		sites2<-as.data.frame(matrix(ncol=3,nrow=27))
			colnames(sites2)<-c('siteID','geo.x','geo.y')
			sites2$siteID<-paste0('h',c(1:27))
			sites2$geo.x<-1
			sites2$geo.y<-1
			sites3<-st_as_sf(sites2,coords=c('geo.x','geo.y'),crs=27700)		


		for(i in 1:27){
			i<-paste(i)
			h<-paste0('h',i)
			g<-st_geometry(hex[i,][2])
			p1<-st_sample(hex[i,],type='random',size=1)
			sites3$geometry[sites3$siteID==h]<-p1
				}
		ggplot()+geom_sf(data=hex,alpha=0.2,fill='cadetblue2')+geom_sf(data=sites3,col='violetred4',alpha=0.3)+theme_bw() 

		sites3%>%print(n=Inf)

		#st_write(hex,'hexagon_grid_edna_sampling_kernios_nocoastal.shp')
		#st_write(sites3,'sampling_sites_randomgen_onePERhexagon_kernios.shp')


	## Making a map of hexagon area for communicating with partners about location support.

	grid<-st_as_sf(st_read('hexagon_grid_edna_sampling_kernios_nocoastal.shp'),crs=27700)
	spots<-st_as_sf(st_read('sampling_sites_randomgen_onePERhexagon_kernios.shp'),crs=27700)
	cw<-st_as_sf(st_read('kernow_ios_coastalwater.shp'))
	ios<-st_as_sf(st_read('islesofscilly.shp'))
	kern<-st_as_sf(st_read('kernow.shp'))
	kernios<-bind_rows(kern,ios)
		
		map<-ggplot()+geom_sf(data=grid,alpha=0.2,fill='cadetblue2')+geom_sf(data=kernios,fill='grey72',alpha=0.5)+geom_sf(data=cw,fill='dodgerblue2',alpha=0.2)+theme_bw() 
		ggsave(map,file='map_with_grid_NO_sampling_spots.png',device='png',units='in',height=6,width=6,dpi=800)














