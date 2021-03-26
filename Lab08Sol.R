#
objects()
rm(list=objects())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,soilDB,rgeos,data.table)
options("download.file.extra"="-L -k")
options("download.file.method"="curl")

# 
#  Always nice to start with hydro and DEM data to warm ourselves up
# 
# Get our gold standard flow data from USGS 0205551460 
# LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = "2019-01-01")

# We want Q in mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
# Get our LatLon to UTM requirements settled
proj4_utm = paste0("+proj=utm +zone=", trunc((180+myflowgage$declon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
latlon <- cbind(myflowgage$declon,myflowgage$declat)
myflowgage$gagepoint_ll <- SpatialPoints(latlon)

proj4string(myflowgage$gagepoint_ll)=proj4_ll
myflowgage$gagepoint_utm=spTransform(myflowgage$gagepoint_ll,crs_utm)
#
# Finding your basin can really help save time in processing
#
url=paste0("https://www.google.com/maps/@",
           myflowgage$declat,",",myflowgage$declon,",18z")
browseURL(url)

searchlength=sqrt(myflowgage$area*6)*1000  # meters per km
searchlength=12000 # meters
bboxpts=myflowgage$gagepoint_utm@coords
bboxpts=rbind(bboxpts,bboxpts+c(-searchlength,searchlength))
# Outlet in Lower Right, so from outlet, -x and +y
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# grabbing and projecting latlon DEM to UTM
mydem=get_aws_terrain(locations=bboxpts@coords, 
                      z = 12, prj = proj4_utm)
mydem=projectRaster(mydem,crs=crs_utm)
writeRaster(mydem,filename = "mydem.tif",overwrite=T)
plot(mydem)
plot(bboxpts,add=T)

### And that was all to get an updated background DEM, since we are
#   importing previous lab's data
url="https://github.com/vtdrfuka/BSE5304/raw/main/Lab08Data.zip"
download.file(url,"Lab08Data.zip")
unzip("Lab08Data.zip")
load("Lab08.RData")
plot(TIC)
# Just as a reminder as to how to get your soils data, you created a shape
# file from your delineated basin, and uploaded it as an AOI template to
url="https://github.com/vtdrfuka/BSE5304/raw/main/mydemw_poly_gdal.zip"
download.file(url,"mydemw_poly_gdal.zip")
unzip("mydemw_poly_gdal.zip")
# Open the WebSoilSurvey site to: 
browseURL("https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx")
# "Create AOI from a zipped shapefile" as we did in previous weeks
# Open "Download Soils Data" Tab
# "Create Download Link" in lower right hand corner
# Right-Click on download link and "Copy Link Address" and 
# paste into a url object:
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/emdfriehdwyjort5pohresl3/wss_aoi_2021-03-25_22-10-28.zip"
download.file(url,"wss_aoi.zip")
unzip("wss_aoi.zip")
mysoil_ll=readOGR("wss_aoi_2021-03-25_22-10-28/spatial/soilmu_a_aoi.shp")
View(mysoil_ll@data)
mysoil_utm = spTransform(mysoil_ll,crs_utm)
plot(mysoil_utm,add=T)

chorizon_tab=read.csv("wss_aoi_2021-03-25_22-10-28/tabular/chorizon.txt",sep="|",header=F)
component_tab=read.csv("wss_aoi_2021-03-25_22-10-28/tabular/comp.txt",sep="|",header=F)
View(chorizon_tab)
View(component_tab)
#soildb_US_2003_mdb=list()
#for (tblname in mdb_tables("wss_aoi_2021-03-25_22-10-28/soildb_US_2003.mdb")){
#  soildb_US_2003_mdb[[tblname]]=read_mdb("wss_aoi_2021-03-25_22-10-28/soildb_US_2003.mdb",table = tblname)
#}
# save(soildb_US_2003_mdb, file="soildb_US_2003_mdb.RData")
# Download Soils MDB based Template
url="https://github.com/vtdrfuka/BSE5304/raw/main/soildb_US_2003_mdb.RData"
download.file(url,"soildb_US_2003_mdb.RData")
load("soildb_US_2003_mdb.RData")
View(soildb_US_2003_mdb[["chorizon"]])
colnames(chorizon_tab)=colnames(soildb_US_2003_mdb[["chorizon"]])
colnames(component_tab)=colnames(soildb_US_2003_mdb[["component"]])
View(hru_table)
hru_table_soils=merge(hru_table,component_tab)
hru_table_soils=merge(hru_table_soils,chorizon_tab)
View(hru_table_soils)

MUSLE=aggregate(hru_table_soils,list(hru_table_soils$TIclass),mean,na.rm=T)
#
# Easiest first! Eq. 4:1.1.15 Course Fragment Factor
MUSLE$CFRG=exp(-0.053*MUSLE$frag3to10_r)
MUSLE
#
# LSusle is calculated using eq. 4.1.12
MUSLE$alpha=atan(MUSLE$slp/100)
MUSLE$LSm=.6*(1-exp(-35.835*MUSLE$slp/100))
MUSLE$LS=(MUSLE$slopelenusle_r/22.1)^MUSLE$LSm * (65.41*sin(MUSLE$alpha)^2+4.56*sin(MUSLE$alpha)+0.065)
#
# Pusle Homework 2
?predict
Pusle_df=data.frame(slope=c(0,1.5,4,7,10.5,14.5,18.5,23,100),pusle=c(.6,.6,.5,.5,.6,.7,.8,.9,.9),maxlen=c(122,122,91,61,37,24,18,15,15))
new=data.frame(slope=tan(MUSLE$slp)*100)
attach(Pusle_df)
predict(lm(pusle~slope),new)
MUSLE$Pusle=predict(lm(pusle~slope),new)
detach(Pusle_df)

#MUSLE$Pusle=.50
#
# Cusle
MUSLE$Cusle=.20
#
# Kusle
MUSLE$Kusle=0.28
#
# Build a constant for those we are not changing day to day
attach(MUSLE)
MUSLE$KCPLSCFRG118=11.8*Kusle*Cusle*Pusle*LS*CFRG
detach(MUSLE)
MUSLE # Make sure values look correct, Pusle, Cusle, Kusle
#
# Now we need to use each of the TIClass Q solutions from Lab06 to calculate
# peak flows (qpeak) and complete the MUSLE Sediment Loss for each class.
# Run Model
#
# Now we need to use each of the TIClass Q solutions from Lab07 to calculate
# peak flows (qpeak) and complete the MUSLE Sediment Loss for each class.
# Run Model
#download.file("https://github.com/vtdrfuka/BSE5304/raw/main/grabdata.R",
#              "grabdata.R")
#download.file("https://github.com/vtdrfuka/BSE5304/raw/main/functions.R",
#              "functions.R")
# Use the estimated S for our watershed (Lab06)
Sest = 157
# We will split into 5 VSA areas represented by 5 TI Classes
nTIclass=5
VSAsol=data.table(TIClass=seq(from=nTIclass,to=1),
                  As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]
#
VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1
# Calculate TI Class localized sigma and Curve Number
VSAsol[,sigma:=Sest*sSratio]
VSAsol[,CN:=25400/(sigma+254)]
VSAsol

VSAParams=merge(VSAsol,MUSLE,by.x="TIClass",by.y="TIclass")
View(VSAParams)
modeldata$HillslopeAboveExcess=0
TIC01=modeldata
TIC02=modeldata
TIC03=modeldata
TIC04=modeldata
TIC05=modeldata
# We now need our watershed CN model which just happens to be in 
# our github repository!
url="https://github.com/vtdrfuka/BSE5304/raw/main/functions.R"
download.file(url,"functions.R")
source("functions.R")
# For TIC01 CNavg=VSAParams$CN[1] but confirm
TIC01 = CN_Model(fnc_CNModel = TIC01, CNavg=VSAParams$CN[1])
TIC01$qpeak=TIC01$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec
TIC01$sed=(TIC01$Qpred*TIC01$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[1]    # Eq. 4:1.1.1 SWAT Theory
# Route the water and continue
TIC02$HillslopeAboveExcess=TIC01$Qpred
TIC02 = CN_Model(fnc_CNModel = TIC02, CNavg=VSAParams$CN[2])
TIC02$qpeak=TIC02$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec
TIC02$sed=(TIC02$Qpred*TIC02$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[2]  

TIC03$HillslopeAboveExcess=TIC02$Qpred
TIC03 = CN_Model(fnc_CNModel = TIC03, CNavg=VSAParams$CN[3])
TIC03$qpeak=TIC03$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec
TIC03$sed=(TIC03$Qpred*TIC03$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[3]  

TIC04$HillslopeAboveExcess=TIC03$Qpred
TIC04 = CN_Model(fnc_CNModel = TIC04, CNavg=VSAParams$CN[4])
TIC04$qpeak=TIC04$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec
TIC04$sed=(TIC04$Qpred*TIC04$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[4]  

TIC05$HillslopeAboveExcess=TIC04$Qpred
TIC05 = CN_Model(fnc_CNModel = TIC05, CNavg=VSAParams$CN[5])
TIC05$qpeak=TIC05$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec
TIC05$sed=(TIC05$Qpred*TIC05$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[5]  


plot(TIC05$date,cumsum(TIC05$sed),type="l",lty=1,col=1,ylab="Sed Tons",xlab="date")
lines(TIC04$date,cumsum(TIC04$sed),type="l",lty=1,col=2,ylab="Sed Tons",xlab="date")
lines(TIC03$date,cumsum(TIC03$sed),type="l",lty=1,col=3,ylab="Sed Tons",xlab="date")
lines(TIC02$date,cumsum(TIC02$sed),type="l",lty=1,col=4,ylab="Sed Tons",xlab="date")
lines(TIC01$date,cumsum(TIC01$sed),type="l",lty=1,col=5,ylab="Sed Tons",xlab="date")
legend("topleft",legend=c("TIC05","TIC04","TIC03","TIC02","TIC01"),col=c(1,2,3,4,5),lty=1,cex=.5)
lines(TIC05$date,cumsum(TIC05$sed),type="l",lty=2,col=1,ylab="Sed Tons",xlab="date")
lines(TIC04$date,cumsum(TIC04$sed),type="l",lty=2,col=2,ylab="Sed Tons",xlab="date")
lines(TIC03$date,cumsum(TIC03$sed),type="l",lty=2,col=3,ylab="Sed Tons",xlab="date")
lines(TIC02$date,cumsum(TIC02$sed),type="l",lty=2,col=4,ylab="Sed Tons",xlab="date")
lines(TIC01$date,cumsum(TIC01$sed),type="l",lty=2,col=5,ylab="Sed Tons",xlab="date")

# This is the point where we saved our data for Week09's Lab
#save(list=c("TIC01","TIC02","TIC03","TIC04","TIC05","TIC","modeldata","hru_table"),file = "Lab08.RData")

# Graduate Homework which is of value
?reclassify
# reclassify the values into three groups 
# all values > 0 and <= 0.25 become 1, etc.
m <- c(1, 365*mean(TIC01$sed),2, 365*mean(TIC02$sed),3, 
       365*mean(TIC03$sed),4, 365*mean(TIC04$sed),5, 365*mean(TIC05$sed))
rclmat <- matrix(m, ncol=2, byrow=TRUE)
RTIC <- reclassify(TIC, rclmat)

plot(RTIC)

#
# Ready for Lab09!
#



