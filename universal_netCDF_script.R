rm(list=ls())
## The only potential things that should need changed for use with other netCDF files
## 1. Tree ring data file (Line 16)
## 2. Beginning and ending years for tree ring data (Line 19, numbers)
## 3. netCDF file name and/or var.name (Line 23, quotes)
## 4. Spatial extent (Line 28)
## 5. Line 127-128, naming/order of the raster brick layers; currently for southern hemisphere
## 6. Color ramp (Line 136)
## 7. Coastline file, although if packaged with this file it shouldn't need changed.

library(ncdf4)
library(raster)

##Tree ring data
trDat <- read.csv("mt_read_mean_hp_kbp.csv", header = TRUE)

#Decide start year and end year. Crop data to it
trDat <-trDat[which(trDat$year >= 1979 & trDat$year<= 2011),]

#Bring in mean sea level pressure data
#This file has two variables slp and msl. If only one variable, varname isn't required.
dat <- brick("ERAinterim.nc", varname = "sst")

## set spatial extent for area interested in. Longitude min - max then latitude min - max
## Creating this way allows for other created rasters to recognize as an extent with 
## the proper spatial extent, otherwise have to set xmin, xmax, ymin, ymax individually.
ext <- extent(75.5, 175.5, -80.5, 0.5)

#Spatial crop using extent
datC <- crop(dat, ext)

# Subset rasters based on the tree ring data. Won't exclude partial seasons.
# second subset for that.
# substr(names(dat)....) reads raster name and selects the characters at specified points
# and compares to values in tree ring file
datC <- datC[[which(as.numeric(substr(names(dat), 2, 5)) >= trDat$year[1] & 
                     as.numeric(substr(names(dat), 2, 5)) <= trDat$year[nrow(trDat)])]]
datC <- datC[[-c(1:2, nlayers(datC))]] #subsets based on years in tree ring data file

##set up identifier for seasons to use when doing raster math
starts = rep(seq(1,nlayers(datC)/3, 1), each=3)

##Get Mean of seasonal sea level pressure
datM <- stackApply(datC, starts, mean)

##remove unnecessary files
rm(dat, starts)

## Rename the raster layers to something more readable
## First year not a complete set of seasons, must rename them individually then paste together years and seasons
## for the rest.
names(datM) <- c(paste("MAM", substr(names(datC[[1]]), 2, 5)),
                 paste("JJA", substr(names(datC[[1]]), 2, 5)),
                 paste("SON", substr(names(datC[[1]]), 2, 5)), 
                 paste( rep( c("DJF", "MAM", "JJA", "SON"), (nlayers(datM)-3)/4), 
                        rep(seq(as.numeric(substr(names(datC[[12]]), 2, 5)), 
                                as.numeric(substr(names(datC[[nlayers(datC)]]), 2, 5))
                                , 1), each=4)))
rm(datC)

# If data is not continuous you must
# replace all NA with a number outside of the bounds of data in order to do correlations.
# For temperature, 0 will not work since values do go to 0 or below. 
# Done at this stage because data can be too large.
# Subset Seasons
SON <- values(subset(datM, grep("SON", names(datM), value=T)))
SON[is.na(SON[])] <- -10
DJF <- values(subset(datM, grep("DJF", names(datM), value=T)))
DJF[is.na(DJF[])] <- -10
JJA <- values(subset(datM, grep("JJA", names(datM), value=T)))
JJA[is.na(JJA[])] <- -10
MAM <- values(subset(datM, grep("MAM", names(datM), value=T)))
MAM[is.na(MAM[])] <- -10

#create rasters to ingest the spatial correlations
CorT <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
Cor <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)

## create correlation based on tree ring indices
## Create raster of p values used to create the cropped confidence intervals
for(i in 1:dim(SON)[1]){
  Cor[i] <- cor(x=SON[i,], y = trDat$mr_kbp, method = 'pearson') 
  CorT[i] <- cor.test(x=SON[i,], y = trDat$mr_kbp, method = 'pearson')$p.value 
}

## Create the clipped 95% Confidence intervals
## SpCorT contains p-values. Replace any values >0.05 with NA
CorT[CorT > 0.05] <- NA
## Using the altered data frame, create the cropped field
son1 <- mask(Cor, CorT)
## Remove unnecessary data
rm(SON)

#One less DJF than the rest of the seasons because of split-year season
trDatS <- trDat[-1,]
for(i in 1:dim(DJF)[1]){
  CorT[i] <- cor.test(x=DJF[i,], y = trDatS$mr_kbp, method = 'pearson')$p.value 
  Cor[i] <- cor(x=DJF[i,], y = trDatS$mr_kbp, method = 'pearson') 
}

CorT[CorT > 0.05] <- NA
djf1 <- mask(Cor, CorT)
rm(DJF)

for(i in 1:dim(JJA)[1]){
  Cor[i] <- cor(x=JJA[i,], y = trDat$mr_kbp, method = 'pearson') 
  CorT[i] <- cor.test(x=JJA[i,], y = trDat$mr_kbp, method = 'pearson')$p.value 
}

CorT[CorT > 0.05] <- NA
jja1 <- mask(Cor, CorT)
rm(JJA)

for(i in 1:dim(MAM)[1]){
  Cor[i] <- cor(x=MAM[i,], y = trDat$mr_kbp, method = 'pearson') 
  CorT[i] <- cor.test(x=MAM[i,], y = trDat$mr_kbp, method = 'pearson')$p.value 
}

CorT[CorT > 0.05] <- NA
mam1 <- mask(Cor, CorT)
rm(MAM, trDat, trDatS, datM, i) 

#Stack the new rasters to be displayed.
Seasons <- brick(jja1, son1, djf1, mam1)
names(Seasons) <- c("Winter", "Spring", "Summer", "Autumn")
rm(mam1, jja1, son1, djf1, Cor, CorT)

#Load in a shapefile and crop for the region of interest 
library(rgdal)
library(rgeos)
coast_shapefile <- crop(readOGR("ne_10m_coastline.shp"), ext)

#Create color ramps for mapping and number of colors to use
library(colorRamps)
col5 <- colorRampPalette(c('#08519c', 'gray96', "#fee0d2", "firebrick3"))

#Load the lattice packages to display the maps
library(rasterVis)
library(gridExtra)

#Plot correlation map using the color ramp and levels using lattice
#Layout dictates how many plots in row, col format
#col.region is the color ramp previously created, pretty makes the color breaks happen logically
#colorkey is information about the legend, wanted bottom so have to give it the space
#par.settings is various graphical settings outside of plots
# +layer(...) adds the coastlines onto the map
levelplot(Seasons, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Seasonal Sea Surface Temperature 1979 - 2009",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))
