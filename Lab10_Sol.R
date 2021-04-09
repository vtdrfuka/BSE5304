pacman::p_del(httr,EcoHydRology,curl,elevatr,raster,rgdal,
               data.table,foreign,maptools,dataRetrieval,gdistance)
pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
               data.table,foreign,maptools,dataRetrieval,gdistance)

dir.create("~/Week10/")
setwd("~/Week10/")

#
# Note we have a new library to access USGS Waterdata
# https://owi.usgs.gov/R/dataRetrieval.html
# https://owi.usgs.gov/R/training-curriculum/usgs-packages/dataRetrieval-readNWIS/
#
?dataRetrieval  # Review the man page for this package
?readNWISuv
?readNWISdv
?readNWISdata
#
# Need to figure out which data to download. 
# https://nwis.waterdata.usgs.gov/nwis/pmcodes?radio_pm_search=param_group&pm_group=All+--+include+all+parameter+groups&pm_search=&casrn_search=&srsname_search=&format=html_table&show=parameter_group_nm&show=parameter_nm&show=casrn&show=srsname&show=parameter_units
# 
##############################################
# 0205551460 LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
##############################################
make_usgs_gage_list=function(siteNo = "0205551460",
   parameterCd = c("00060","00065"),
   start.date = "2017-05-01",  # Not frozen to not frozen
   end.date = "2017-11-01") {  # to still not frozen

  el_listo_amazing=list()   # Organize the data in a nice list as in previous labs
  el_listo_amazing[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
#  agency_cd	site_no        	dateTime X_00060_00000 X_00060_00000_cd
#1  	USGS 0205551460 2017-05-01 04:00:00      	6.38            	A
#2  	USGS 0205551460 2017-05-01 04:05:00      	6.38            	A
#  X_00065_00000 X_00065_00000_cd tz_cd
#1      	2.74            	A   UTC
#2      	2.74            	A   UTC
#
# And of course we want to work in SI units so:
  el_listo_amazing$flowdata$depth_m=el_listo_amazing$flowdata$X_00065_00000*0.3048
# m/ft depth
  el_listo_amazing$flowdata$cms=el_listo_amazing$flowdata$X_00060_00000*.02832
# m3/ft3 flow
#
# Let's add in the USGS gage site information to the list and inspect
  el_listo_amazing[["site"]]=readNWISsite(siteNo)
  class(el_listo_amazing$site$dec_lat_va)
#
# Set the Manning Coefficient in the USGS Gage's Site Table
#
  el_listo_amazing$site$man_n=.035/1.49
#
# Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
  coordinates(el_listo_amazing$site)=~dec_long_va+dec_lat_va
#
#Might be nice to have a function that gets these data, no?
  return(el_listo_amazing)
}
USGS02056000=make_usgs_gage_list("02056000")
USGS0205551460=make_usgs_gage_list("0205551460")
USGS02055100=make_usgs_gage_list("02055100")
USGS02055000=make_usgs_gage_list("02055000")
USGS02054530=make_usgs_gage_list("02054530")

# Data from DEM:
  ab_ll=rbind(USGS02056000$site,
                USGS0205551460$site,
                USGS02055100$site,
                USGS02055000$site,
                USGS02054530$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                     trunc((180+coordinates(USGS02055000$site)[1])/6+1), 
                     " +datum=WGS84 +units=m +no_defs")
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords
mydem=get_aws_terrain(locations=ab_utm@coords, 
                        z = 11, prj = proj4_utm,expand=1)
#
# Lets plot the DEM and the gage locations so we can guess 
# what gages connect with what gages
#
plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
# From Lab02, I know I can get an overview of streams with the 
# USGS H
url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_03010101_HU8_Shape.zip"
curl_download(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
streams_utm=spTransform(streams,crs_utm)

#
# Breaking things  for educational purposes
USGS02056000$flowdata=USGS02056000$flowdata[,c(1,2,3,4,5,8,10)]
# Assume in the inventory link that for this gage, our Gage height is missing. 

# rating curve is: https://en.wikipedia.org/wiki/Rating_curve
# and use the readNWISrating() function to grab it for this gage
USGS02056000[["rating"]]=readNWISrating(USGS02056000$site$site_no)
plot(USGS02056000$rating$DEP,USGS02056000$rating$INDEP,xlab="DEP",
     ylab="INDEP",ylim=c(0,5),xlim=c(0,3000))
#
# Note that this is very similar to what we saw in the previous gage's results
# and as it turns out, we can use it to estimate a 00065 measurement as 
# we did for the previous gage.
USGS02056000$flowdata$X_00065_00000=approx(USGS02056000$rating$DEP,
  USGS02056000$rating$INDEP, xout = USGS02056000$flowdata$X_00060_00000, ties = min)$y
points(USGS02056000$flowdata$X_00060_00000,USGS02056000$flowdata$X_00065_00000,
  col="red")
min(USGS02056000$flowdata$X_00065_00000,na.rm=T)*0.3048
#
# HW1 Answer
# 

USGS02056000$flowdata$depth_m=USGS02056000$flowdata$X_00065_00000*0.3048-
                (min(USGS02056000$flowdata$X_00065_00000,na.rm=T)*0.3048)
USGS0205551460$flowdata$depth_m=USGS0205551460$flowdata$X_00065_00000*0.3048-
            (min(USGS0205551460$flowdata$X_00065_00000,na.rm=T)*0.3048)
# m/ft depth
#
# A quick readthrough of the Example 1: Hiking around Maunga Whau
# in the package vignette. 
# vignette("Overview", package = "gdistance")
# Set the starting and ending locations
# determine the river reach length and slope using the gdistance package.
#
lazydan=function(siteA=USGS0205551460,siteB=USGS02056000){
#  siteA=USGS0205551460;siteB=USGS02056000
  A=SpatialPoints(siteA$site)# Up gradient site Lick Run
  B=SpatialPoints(siteB$site) # Down gradient site ROA River atNiagara
  proj4string(A)=proj4_ll
  proj4string(B)=proj4_ll
  A_utm=spTransform(A,crs_utm)
  B_utm=spTransform(B,crs_utm)
  # Cut the DEM down to a more manageable size
  cropmydem=crop(mydem,extend(extent(ab_utm),600))
  cropmydem=trim(cropmydem)
  cropmydem=cropmydem*1000.0
  plot(cropmydem)
  plot(ab_utm,add=T)
  # Set up the weighting functions
  altDiff <- function(x){x[2] - x[1]}
  hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
  slope <- geoCorrection(hd)
  adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
  speed <- slope
  speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
  Conductance <- geoCorrection(speed)
  # Find and plot the flow path
  AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
  plot(AtoB,add=T)
  cropstreams_utm=crop(streams_utm,cropmydem)
  plot(cropstreams_utm,col="blue",add=T)
  plot(AtoB,add=T)
#  SpatialLinesLengths(AtoB)
  siteA$site$L=SpatialLinesLengths(AtoB) # km to m
  siteA$site$L # reach length in m
  #
  #
  # Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
  # estimate of slope than taking the point slopes at the gage site
  #
  siteA$site$slope=(extract(mydem,A_utm)-
       extract(mydem,B_utm))/siteA$site$L
  siteA$site$slope
  # So now we have flow depth (y "$depth_m"), manning's n ("$man_n"), Q ("$cms"), and slope ("$slope") rearrange to solve for B
  # B=(n*Q)/(y^(5/3)*sqrt(So))
  siteA$flowdata$B=(siteA$site$man_n*
         siteA$flowdata$cms)/(siteA$flowdata$depth_m^(5/3)*
         sqrt(siteA$site$slope))
  #head(siteA$flowdata)
  #  agency_cd	site_no        	dateTime X_00060_00000 X_00060_00000_cd
  #1  	USGS 05267000 2017-05-01 04:00:00      	6.38            	A
  #2  	USGS 05267000 2017-05-01 04:05:00      	6.38            	A
  #  X_00065_00000 X_00065_00000_cd tz_cd   	cms  depth_m    	B
  #1      	2.74            	A   UTC 0.1806816 0.835152 0.103032
  #2      	2.74            	A   UTC 0.1806816 0.835152 0.103032
  #
  # Lets look at how B changes with flow. 
  plottitle=paste0(siteA$site$station_nm," to ",siteB$site$station_nm)
  plot(siteA$flowdata$dateTime,siteA$flowdata$B, main=plottitle)
  # Does this seem reasonable (...like order of magnitude reasonable)? You can 
  # perform a quick and dirty check using google earth and measuring the channel 
  # width in a few places.
  #
  plot(siteA$flowdata$cms,siteA$flowdata$depth_m, main=plottitle)

  siteA$flowdata$ck =
    5/3*sqrt(siteA$site$slope)/siteA$site$man_n*
    (siteA$flowdata$depth_m^(2/3))
  #
  siteA$flowdata$dt =
    siteA$site$L/siteA$flowdata$ck
  
  plot(siteA$flowdata$dateTime,siteA$flowdata$dt)
  siteA$flowdata$outTime=siteA$flowdata$dateTime+
    siteA$flowdata$dt
  
  return(siteA)
}
USGS0205551460=lazydan()
USGS02055100=lazydan(siteA = USGS02055100)
USGS02054530=lazydan(siteA = USGS02054530)
USGS02055000=lazydan(siteA = USGS02055000)


lazydan2=function(site=USGS0205551460,wavegrowth=1.1){
  site$flowdata$newwave=
  site$flowdata$cms *wavegrowth <
  data.table::shift(site$flowdata$cms)
# Add plot of the point found
len=length(site$flowdata$newwave)
site$flowdata$newwave[is.na(site$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  if(site$flowdata$newwave[i]==T &
     site$flowdata$newwave[i-1]==T){
    site$flowdata$newwave[i]=F
  }
}
summary(site$flowdata$newwave)
which(site$flowdata$newwave == TRUE)


}

# Find the time locations where waves begin
lazydan2(USGS0205551460,wavegrowth = 1.16)
location=15113
maxcms=max(USGS0205551460$flowdata$cms[location:(location+600)],na.rm = T)
plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$cms/maxcms,
       type="l",xlim=c(USGS0205551460$flowdata$dateTime[location],
                       USGS0205551460$flowdata$dateTime[location+600]),ylim=c(0,1.01))
points(USGS0205551460$flowdata$dateTime[USGS0205551460$flowdata$newwave],
       USGS0205551460$flowdata$cms[USGS0205551460$flowdata$newwave]/max(USGS0205551460$flowdata$cms,na.rm = T),col=2)
lines(USGS0205551460$flowdata$outTime,USGS0205551460$flowdata$cms/maxcms,col=2)


lazydan2(USGS02055100,wavegrowth = 1.08)
maxcms=max(USGS02055100$flowdata$cms[location:(location+600)],na.rm = T)
lines(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms/maxcms,
     type="l",xlim=c(USGS02055100$flowdata$dateTime[location],
                     USGS02055100$flowdata$dateTime[location+600]),ylim=c(0,1.01))
lines(USGS02055100$flowdata$outTime,USGS02055100$flowdata$cms/maxcms,col=2)

lazydan2(USGS02054530,wavegrowth = 1.03)
maxcms=max(USGS02054530$flowdata$cms[location:(location+600)],na.rm = T)
lines(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms/maxcms,
      type="l",xlim=c(USGS02054530$flowdata$dateTime[location],
                      USGS02054530$flowdata$dateTime[location+600]))
lines(USGS02054530$flowdata$outTime,USGS02054530$flowdata$cms/maxcms,col=2)

lazydan2(USGS02055000,wavegrowth = 1.08)
maxcms=max(USGS02055000$flowdata$cms[location:(location+600)],na.rm = T)
lines(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms/maxcms,
      type="l",xlim=c(USGS02055000$flowdata$dateTime[location],
                      USGS02055000$flowdata$dateTime[location+600]))
lines(USGS02055000$flowdata$outTime,USGS02055000$flowdata$cms/maxcms,col=2)




