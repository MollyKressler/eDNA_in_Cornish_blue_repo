
### Sampling locations around Cornwall

## created by Molly Kressler on 17 February 2023 

pacman::p_load(sf,tidyverse,ggplot2,ggsn,lubridate,viridis,flextable,patchwork)
setwd('/Users/mollykressler/Documents/EDNA/data_edna')

############## Load at start 

	area<-st_as_sf(st_read('startingboundarybox_forsampling_landRemoved.shp'),crs='WGS84')
	cw<-st_as_sf(st_read('kernow_ios_coastalwater.shp'),crs='WGS84')
	ios<-st_as_sf(st_read('islesofscilly.shp'),crs='WGS84')
	kern<-st_as_sf(st_read('kernow.shp'),crs='WGS84')
	kernios<-bind_rows(kern,ios)
	kernios2<-st_combine(kernios)

	locations<-st_as_sf(st_read('kernow_locations.kml'))%>%st_transform('WGS84')

	Falharbourlimit<-st_as_sf(st_read('FalHarbourShapefiles/harbour_limits_FHC.shp')) # Fal Harbour Shapefles 
		Falharbourlimit$id<-'FalH' 
		
		# majority is within the coastal waters. 


	hex<-st_as_sf(st_read('hexgrid_edna/hexagon_grid_edna_sampling_kernios_nocoastal.shp')) # proposed hexagon grid
		# update the shapefle as partners agree to areas and remove hexagons to out to sea
		#	 st_write(hex,'hexgrid_edna/hexagon_grid_edna_sampling_kernios_nocoastal.shp',append=FALSE)

	area_with_land<-ggplot()+geom_sf(data=area,alpha=0.2,fill='cadetblue4')+geom_sf(data=kernios,alpha=0.8,fill='grey82')+theme_bw() 
	hex_with_land<-ggplot()+geom_sf(data=hex,alpha=0.2,fill='cadetblue4')+geom_sf(data=kernios,alpha=0.8,fill='grey82')+theme_bw() 
		# ggsave(area_with_land,file='comms_figures/samplingarea_nogrd_withland.png',device='png',units='in',height=7,width=7,dpi=900)




#################################
## - eDNA sampling locations 
#################################
##  REPEAT WITH NEW DATA DOWNLOADS
#################################
## created by Molly Kressler :: September 2023
#################################

## Load at Start
	 ext <- read.csv('data_downloads_notpublic/extraction_data_edna_asofSEPT2023.csv')
	 samp <- read.csv('data_downloads_notpublic/samplingdata_edna_asofSEPT2023.csv')%>%rename(Latitude=longitude,Longitude=latitude)%>%
	 	mutate(Latitude = as.character(Latitude),Longitude = as.character(Longitude))
	 head(samp)

## Clean and add location and eventID to each eppendorf row in extraction.csv
	 e2 <- ext %>%
	 	rename(extraction.date=date)%>%
	 	mutate(sampledat = str_split(sampleID,'-',simplify=TRUE)[,1],
	 		type = str_split(sampleID,'-',simplify=TRUE)[,2],
	 		AorB = str_split(sampleID,'-',simplify=TRUE)[,3], 
	 		method.subsample = str_split(sampleID,'-',simplify=TRUE)[,4])%>%
	 	mutate(eventID = paste0(sampledat,"-",type,"-",AorB), .before=sampleID)%>%
	 	mutate(sampling.date=dmy(sampledat),month=month(sampling.date))%>%
	 	mutate(across(contains('.quantity'),~replace(.,is.na(.),0)))%>%
	 	mutate(across(contains('.presence'),~replace(.,is.na(.),0)))%>%
	 	filter(!is.na(dnaCont)) # empty last row
	 	head(e2)
	 	sapply(e2,class)

## 	Left_join locations from samp to e2 using the event ID

	d1 <- left_join(e2, samp%>%dplyr::select(eventID,Latitude,Longitude,recorded.by,vessel,passengers,time.timeIN,timeOUT), by = 'eventID', multiple = 'all')%>%
		mutate(timeOUT = ifelse(is.na(timeOUT),time.timeIN,timeOUT))%>%
	 	mutate(method.type = case_when(startsWith(as.character(method.subsample),'W') ~ 'waterbottle', TRUE ~ 'metaprobe'))
	 	## mutate for timeOUT fills in the timeout column for water bottles to match the time.timeIN value. 
	head(d1)
	summary(d1)



## Save as CSV and SHP
	
	write.csv(d1,'eDNA_data_meta_and_qPCR_ASOF_15september2023.csv') # extraction + meta
	
	sf1 <- st_as_sf(d1,coords=c('Longitude','Latitude'),crs='WGS84') #extraction + meta
		unique(sf1$geometry) # will be FEWER than in samp/sampsf1, but equal to how many locations you've done extractions and qPCR for 

	sampsf1 <- 	st_as_sf(samp,coords=c('Longitude','Latitude'),crs='WGS84')%>%
	 	mutate(type = str_split(eventID,'-',simplify=TRUE)[,2],
	 		sampledat = str_split(eventID,'-',simplify=TRUE)[,1],.before='geometry')%>%
	 	mutate(sampling.date=dmy(sampledat), month = as.character(month(sampling.date)))%>%
	 	mutate(method.type = case_when(startsWith(as.character(sampleID),'W') ~ 'waterbottle', TRUE ~ 'metaprobe'))		# sampling data
	 	sampsf1
		unique(sampsf1$geometry)

 ########
 ## Plots 
 ########

	month.palette <- c('#461582','#6a0078','#85006b','#9a005c','#a9004b','#b20039','#b82711','#ac5600','#a06800','#7f8800','#49a306','#1e0e9e')
	sixmonths.palette <- c('#1e0e9e','#6a0078','#a9004b','#b20039','#b82711','#7f8800')

  ### Descriptive Stats
	## Samples per month - creates a bar chart of samples collected per month with a table showign the number of each type per month 
		samplespermonth <- ggplot(data=sampsf1,aes(x=month,col=month,fill=month,alpha=0.8))+
			geom_histogram(stat='count')+
			scale_colour_manual(values=sixmonths.palette)+
			scale_fill_manual(values=sixmonths.palette)+
			scale_x_discrete(labels=c('April','May','June','July','August','September'))+
			xlab('Month (2023)')+ ylab('Samples (count)')+
			guides(col='none',fill='none',alpha='none')+
			scale_y_continuous(breaks=c(3,6,9,12),limits=c(0,14))+
			theme_bw()
		methodtypepermonth_controlledformultiplerowsperday <- sampsf1%>%group_by(method.type,month)%>%dplyr::select(method.type,month)%>%tally()%>%ungroup()%>%as.tibble()%>%dplyr::select(-geometry)%>%pivot_wider(names_from=month,values_from=n)
		methodtypepermonth.table<-methodtypepermonth_controlledformultiplerowsperday%>%
			flextable()%>%
			set_header_labels('method.type'='Method','4'='April','5'='May','6'='June','7'='July','8'='August','9'='September')%>%
			bold(bold=TRUE,part='header')%>%
			align(align = 'center', part = 'all')%>%
			font(fontname = 'Arial', part = 'all')%>%
			color(color='black',part='all')%>%
			fontsize(size = 10, part = 'all')%>%
			autofit()
		met <- flextable::gen_grob(methodtypepermonth.table)
		outreachANDvessels_permonth_apriltosept2023 <- cowplot::plot_grid(samplespermonth,met,nrow=2,ncol=1,rel_heights=c(6,1),scale=c(1,.75))
		ggsave(outreachANDvessels_permonth_apriltosept2023,file='samplesPERmonth_sampling_april2sept2023.png',device=png,units='in',height=5,width=5,dpi=600)
	
	## Outreach per month - creates a bar chart with number of passengers (unique) exposed to the project per month, and the number of participating (unique) vessel per month 

		passengerspermonth_controlledformultiplerowsperday <- sampsf1%>%group_by(eventID)%>%slice_head(n=1)%>%ungroup()%>%group_by(month)%>%tally(passengers)

		outreachpermonth<-ggplot(data=passengerspermonth_controlledformultiplerowsperday,aes(x=month,y=n,col=month,fill=month,alpha=0.8))+
			geom_col()+
			scale_colour_manual(values=sixmonths.palette)+
			scale_fill_manual(values=sixmonths.palette)+
			scale_x_discrete(labels=c('April','May','June','July','August','September'))+
			xlab('Month (2023)')+ ylab('Outreach (count)')+
			guides(col='none',fill='none',alpha='none')+
			theme_bw()
			# have to divide by two because the data contain 
		vesselspermonth_controlledformultiplerowsperday <- sampsf1%>%group_by(eventID)%>%slice_head(n=1)%>%ungroup()%>%group_by(vessel)%>%group_by(month)%>%tally()
		vesselspermonth.table<-vesselspermonth_controlledformultiplerowsperday%>%
			as.tibble()%>%
			dplyr::select(month,n)%>%
			pivot_wider(names_from=month,values_from=n)%>%
			flextable()%>%
			set_header_labels('4'='April','5'='May','6'='June','7'='July','8'='August','9'='September')%>%
			bold(bold=TRUE,part='header')%>%
			align(align = 'center', part = 'all')%>%
			font(fontname = 'Arial', part = 'all')%>%
			color(color='black',part='all')%>%
			fontsize(size = 10, part = 'all')%>%
			autofit()
		out <- flextable::gen_grob(vesselspermonth.table)
		outreachANDvessels_permonth_apriltosept2023 <- cowplot::plot_grid(outreachpermonth,out,nrow=2,ncol=1,rel_heights=c(6,1),scale=c(1,.75))
		ggsave(outreachANDvessels_permonth_apriltosept2023,file='outreach_and_vessels_permonth_april2september2023.png',device=png,units='in',height=5,width=5,dpi=600)

	## DNA concentration over time

		dnacont.overtime<-ggplot(data=d1,aes(x=sampling.date,y=dnaCont,group=method.type,col=method.type,fill=method.type))+
				geom_smooth()+
				scale_fill_manual(values=c('cadetblue2','violetred4'))+
				scale_color_manual(values=c('cadetblue2','violetred4'))+
				scale_x_continuous(labels=c('April','May','June','July','August','September'))+
				scale_y_continuous(breaks=c(0,20,40,60,80,100))+
				xlab('Sampling Date (2023)')+ ylab('DNA Concentraton (double-stranded)')+
				theme_bw()
		ggsave(dnacont.overtime,file='dnaCont_oversamplingperiod_april2september2023.png',device=png,units='in',height=5,width=6.5,dpi=700)

  ### Maps 
	## LOCATIONS
		##  using sampling data only (more samples collected than processed at most times)
		locations_sampling <- ggplot()+
			geom_sf(data=kernios2,alpha=0.8,fill='grey72')+
			geom_sf(data=sampsf1, col='cadetblue4',size=2,pch=18,alpha=0.45)+
			theme_bw() 
			ggsave(locations_sampling,file='locations_sampling_april2sept2023.png',device=png,units='in',height=5,width=5,dpi=600)
		
		locations_sampling_CSorR <- ggplot()+
			geom_sf(data=kernios2,alpha=0.8,fill='grey72')+
			geom_sf(data=sampsf1%>%filter(type=='R'), col='cadetblue4',size=2,pch=18,alpha=0.45)+
			geom_sf(data=sampsf1%>%filter(type=='CS'), col='violetred4',size=2,pch=16)+
			theme_bw() 

		locations_sampling_bymonth <- ggplot()+
			geom_sf(data=kernios2,alpha=0.8,fill='grey72')+
			geom_sf(data=sampsf1, aes(colour=month),size=3,pch=18,alpha=0.75)+
			scale_colour_manual(values=sixmonths.palette)+
			theme_bw() 

		##  using extraction + metadata sf
		location_extractions <- ggplot()+
			geom_sf(data=kernios2,col='grey72')+
			geom_sf(data=sf1, col='cadetblue4',size=2,pch=18,alpha=0.45)+
			theme_bw() 


	## DENSITY OF DETECTIONS, method specific then species specfic (so 12 plots in total) - creates a map plot with points of sampling locations (grouped by sampleID), where the size of the point ~ the relative abundance of a species  
		## data = extraction data frame.

### YOU LEFT OFF HERE



	## SPECIES PROPORTION OF DETECTIONS, method specific (2 plots in total) - creates a map plot with points of sampling locations (grouped by sampleID), where the point is a pe chart, and each species is a different colour, and the slice width is ~ the relative abundance of a species at THAT location. 






#################################
#################################
## - Partner areas 
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





#################################
#################################
## - MAKING SHAPEFILES, do not do every open. 
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














