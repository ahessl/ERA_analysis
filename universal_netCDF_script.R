rm(list=ls())
#setwd("C:/Users/S/Desktop/NetCDFRPlay/ERA")
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
trDat <- read.table("mt_read_mean_hp_kbp.csv", header = TRUE, sep=',')

#Decide start year and end year. Crop data to it
trDat <-trDat[which(trDat$year >= 1979 & trDat$year<= 2011),]
rownames(trDat) <- trDat$year

#Bring in mean sea level pressure data
#This file has two variables slp and msl. If only one variable, varname isn't required.
dat <- brick("ERAinterim.nc", varname = "sst")

## set spatial extent for area interested in. Longitude min - max then latitude min - max
## Creating this way allows for other created rasters to recognize as an extent with 
## the proper spatial extent, otherwise have to set xmin, xmax, ymin, ymax individually.
ext <- extent(80.25, 180, -80.25, -19.50)

#Spatial crop using extent
datC <- crop(dat, ext)

# Subset rasters based on the tree ring data. Won't exclude partial seasons.
# second subset for that.
# substr(names(dat)....) reads raster name and selects the characters (2 through 5) 
# in the names (time in this case) to extract years and compares to values in tree ring file
datC <- datC[[which(as.numeric(substr(names(dat), 2, 5)) >= trDat$year[1] & 
                     as.numeric(substr(names(dat), 2, 5)) <= trDat$year[nrow(trDat)])]]
datC <- datC[[-c(1:2, nlayers(datC))]] #subsets based on years in tree ring data file

##set up identifier for seasons to use when doing raster math
#starts <- rep(seq(1,nlayers(datC)/3, 1), each=3) #3 months for every season
#AMY's alternative to STARTS
library(chron)
#d <- as.Date(gsub(".", "/", substr(names(datC), 2, 11), fixed=T)) #pull out the yr, month, day and fix formatting 

## SC: This doesn't work, and I'm not sure how to fix yet. December is problematic. Working on it.....
yr_mo_dy <- substr(names(datC), 2, 11)
d <- as.Date(gsub(".", '/', yr_mo_dy, fixed = T)) #fix the format by replacing "." with "/"

seasons <- c('DJF', 'MAM', 'JJA', 'SON')[ # select from character vector with numeric vector
  1+((as.POSIXlt(d)$mon+1) %/% 3)%%4]

paste( 1900 + # this is the base year for POSIXlt year numbering 
         as.POSIXlt( d )$year + 
         1*(as.POSIXlt( d )$mon==11) ,   # offset needed for December
       c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
         1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
       , sep="-")

#year + season
yr_season <- paste(seasons, substr(d, 1, 4), sep = ".")

##remove unnecessary files
rm(yr_mo_dy, d)

# ##Get Mean of seasonal sea level pressure using seasons
datM <- stackApply(datC, yr_season, mean) #raster with mean for each season
names(datM) <- unique(yr_season) #Is this more efficient since it doesn't have to call on itself?
  
  #substr(names(datM),7,14) #cleaner name for each seasonal mean


# If data is not continuous you must
# replace all NA with a number outside of the bounds of data in order to do correlations.
# For temperature, 0 will not work since values do go to 0 or below. 
# Done at this stage because data can be too large.
# Also subsets seasons at same time
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
temp <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)

## create correlation based on tree ring indices
## Create raster of p values used to create the cropped confidence intervals

## Presenting the treeCorr Function!!!
## This compresses multitudes of previous coding into a very small package
## Correlation, masking all done in one go. Plot ready after this.
treeCorr <- function(x, y){
  for(i in 1:dim(x)[1]){
      Cor[i] <- cor(x=x[i,], y = y, method = 'pearson') 
      CorT[i] <- cor.test(x=x[i,], y = y, method = 'pearson')$p.value
  }
  CorT[CorT > 0.05] <- NA
  Cor <- brick(Cor)
  temp <- mask(Cor, CorT)
  return(temp)
   
}

## Simplified creating rasters. This isn't completed yet. Making it universal for HP is going to take more work since it's shorter
## Doesn't work quite right yet because of the December issue on lines 49/50

DJF_c <- treeCorr(DJF[,2:ncol(DJF)], trDat$mr_kbp[-c(nrow(trDat))]) #### tg-1
MAM_c <- treeCorr(MAM, trDat$mr_kbp) ##### t
JJA_c <- treeCorr(JJA[,-c(ncol(JJA))], trDat$mr_kbp[-c(1)]) #### t+1
SON_c <- treeCorr(SON[,-c(ncol(SON))], trDat$mr_kbp[-c(1)]) #### t+1


rm(SON, DJF, JJA, MAM, trDat, trDatS, datM, i) 

#Stack the new rasters to be displayed.
Seasons <- stack(DJF_c[[1]], MAM_c[[1]], JJA_c[[1]], SON_c[[1]])
names(Seasons) <- c("DJF", "MAM", "JJA", "SON")
rm(MAM_c, JJA_c, DJF_c, SON_c, Cor, CorT)

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
levelplot(Seasons, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Seasonal SST 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))


# Amy Playing.
# make a mean raster of each season
SON_mean_r <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
SON_mean_r <- setValues(SON_mean_r, rowMeans(SON))

DJF_mean_r <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
DJF_mean_r <- setValues(DJF_mean_r, rowMeans(DJF))

MAM_mean_r <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
MAM_mean_r <- setValues(MAM_mean_r, rowMeans(MAM))

JJA_mean_r <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
JJA_mean_r <- setValues(JJA_mean_r, rowMeans(JJA))

#Stack the new rasters to be displayed.
Season_m <- brick(JJA_mean_r, SON_mean_r, DJF_mean_r, MAM_mean_r)
names(Season_m) <- c("JJAm", "SONm", "DJFm", "MAMm")
rm(JJA_mean_r, SON_mean_r, DJF_mean_r, MAM_mean_r)

#make plot of mean of each season

levelplot(Season_m, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Seasonal Mean SST 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) +
  layer(sp.lines(coast_shapefile))


# End Amy playing.
