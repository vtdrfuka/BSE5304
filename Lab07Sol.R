if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,soilDB,rgeos)
options("download.file.extra"="-L -k")
options("download.file.method"="curl")
# Get our gold standard flow data from USGS 0205551460 
# LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = "2019-01-01")

# We want Q in mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

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

searchlength=sqrt(myflowgage$area*6)*1000  # meters per km
bboxpts=myflowgage$gagepoint_utm@coords
bboxpts=rbind(bboxpts,bboxpts+c(-searchlength,searchlength))
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# grabbing and projecting latlon DEM to UTM
mydem=get_aws_terrain(locations=bboxpts@coords, 
                      z = 12, prj = proj4_utm,src ="aws")

writeRaster(mydem,filename = "mydem.tif",overwrite=T)
plot(mydem)
plot(bboxpts,add=T)
# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
system("pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")

# DInf flow directions
system("mpiexec -n 8 dinfflowdir -ang mydemang.tif -slp mydemslp.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
ang=raster("mydemang.tif")
slp=raster("mydemslp.tif")

# Dinf contributing area
system("mpiexec -n 8 areadinf -ang mydemang.tif -sca mydemsca.tif")
sca=raster("mydemsca.tif")
# Threshold 1/overguesstimate
# mythresh=.01*(raster::cellStats(sca,stat = max))
#
mythresh=200
system(paste0("mpiexec -n 8 threshold -ssa mydemsca.tif -src mydemsrc.tif -thresh ",mythresh))
src=raster("mydemsrc.tif")
plot(src)
plot(myflowgage$gagepoint_utm,add=T)
zoomext=myflowgage$gagepoint_utm@coords
zoomext=rbind(zoomext,zoomext+res(src)*10)
zoomext=rbind(zoomext,zoomext-res(src)*10)
zoomext=SpatialPoints(zoomext,proj4string = crs_utm)
zoom(src,ext=zoomext)

# D8 flow directions
system("mpiexec -n 8 d8flowdir -p mydemp.tif -sd8 mydemsd8.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)

## a quick R function to write a shapefile
#makeshape.r=function(sname="shape",n=1)
#{
#  xy=locator(n=n)
#  points(xy)
#  
  #Point
#  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
#  ddTable <- data.frame(Id=c(1),Name=paste("outlet",1:n,sep=""))
#  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
#  write.shapefile(ddShapefile, sname, arcgis=T)
#}

#makeshape.r("approxoutlets")


outlet=SpatialPointsDataFrame(myflowgage$gagepoint_utm,data.frame(Id=c(1),outlet=paste("outlet",1,sep="")))
writeOGR(outlet,dsn=".",layer="approxoutlets",driver="ESRI Shapefile", overwrite_layer=TRUE)

# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o approxoutlets.shp -om outlet.shp")
outpt=readOGR("outlet.shp")
approxpt=readOGR("approxoutlets.shp")
zoom(src,ext=zoomext)
plot(approxpt, add=T,col="red")
plot(outpt, add=T,col="blue")

# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 

# Threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh 5000")
src1=raster("mydemsrc1.tif")
plot(src1)
plot(outpt, add=T,col="blue")
# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemssa.tif -src mydemsrc1.tif -o outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")

# Build a mask, trim, and crop the basin files
mydemw=raster("mydemw.tif")
plot(mydemw)
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
slp=raster("mydemslp.tif")
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)

mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)

TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)
zoom(TI,ext=zoomext)

pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
plot(TIC)
# A series of plots to show all of the components
#
# DEM - Filled DEM
mybasinfel=crop(fel,mybasinmask)
mybasinfel=mask(mybasinfel,mybasinmask)

par(mfrow = c(2, 2))
plot(mybasinfel)
plot(mybasinfel-mybasindem)
plot(TI)
plot(TIC,col=rainbow(5))

# Lets look closer at our "subbasins"!
mybasindemw=crop(mydemw,mybasinmask)
mybasindemw=mask(mybasindemw,mybasinmask)
plot(mybasindemw)



