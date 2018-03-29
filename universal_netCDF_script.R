
## The only potential things that should need changed for use with other netCDF files
## 1. Tree ring data file (Line 16)
## 2. Beginning and ending years for tree ring data (Line 19, numbers)
## 3. netCDF file name and/or var.name (Line 23, quotes)
## 4. Spatial extent (Line 28)
## 5. Line 127-128, naming/order of the raster brick layers; currently for southern hemisphere
## 6. Color ramp (Line 136)
## 7. Coastline file, although if packaged with this file it shouldn't need changed.
rm(list=ls())
#setwd("C:/Users/S/Desktop/NetCDFRPlay/ERA")

library(ncdf4)
library(raster)

#Bring in mean sea level pressure data
#This file has two variables slp and msl. If only one variable, varname isn't required.
dat <- brick("Reanal20thc/pres.sfc.mon.mean.nc")#, varname = "tmp")

## set spatial extent for area interested in. Longitude min - max then latitude min - max
## Creating this way allows for other created rasters to recognize as an extent with 
## the proper spatial extent, otherwise have to set xmin, xmax, ymin, ymax individually.

#ext <- extent(60, 180, -80, -4)
ext <- extent(-180, 180, -80, 0)

#Spatial crop using extent
datC <- crop(dat, ext)

## Tree ring data
trDat <- read.table("../KBP_South/KBPS_cull_gap.rwl_tabs.txt", header = TRUE)

## Decide start year and end year based on target and tree ring data
#F_yr <- min(as.numeric(substr(names(datC), 2, 5)))
F_yr <- 1979
L_yr <- 1998
#L_yr <- as.numeric(max(trDat$year))


#Crop data to it
trDat <-trDat[which(trDat$year >= F_yr-1 & trDat$year<= L_yr),]

# Subset rasters based on the tree ring data. Won't exclude partial seasons.
# second subset for that.
# substr(names(dat)....) reads raster name and selects the characters (2 through 5) 
# in the names (time in this case) to extract years and compares to values in tree ring file
datC <- datC[[which(as.numeric(substr(names(dat), 2, 5)) >= F_yr & 
                      as.numeric(substr(names(dat), 2, 5)) <= L_yr)]]
datC <- datC[[-c(1:2, (nlayers(datC)-3):nlayers(datC))]] #removes first incomplete season JF and last SON from year

library(chron)

yr_mo_dy <- substr(names(datC), 2, 11)
d <- as.Date(gsub(".", '/', yr_mo_dy, fixed = T)) #fix the format by replacing "." with "/"
##AH: changed $year==12 to $mon<8 (POSIXlt indexes months from 0
##AH: and growing season starts in SEP) 

yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                      as.POSIXlt( d )$year - 
                      1*(as.POSIXlt( d )$mon<8) ,   # offset needed for growing season in SH
                    c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                      1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                    , sep="-")

##remove unnecessary files
rm(yr_mo_dy, d)

# ##Get Mean of seasonal sea level pressure using seasons
datM <- stackApply(datC, yr_season, mean) #raster with mean for each season
names(datM) <- unique(yr_season) #Is this more efficient since it doesn't have to call on itself?

# If data is not continuous you must
# replace all NA with a number outside of the bounds of data in order to do correlations.
# For temperature, 0 will not work since values do go to 0 or below. 
# Done at this stage because data can be too large.
# Also subsets seasons at same time
# AH: PRETTY!

for (i in unique(substring(yr_season, 6))){
  assign(paste0(i), values(subset(datM, grep(i, names(datM), value=T))))
  }

SON[is.na(SON[])] <- -10
DJF[is.na(DJF[])] <- -10
JJA[is.na(JJA[])] <- -10
MAM[is.na(MAM[])] <- -10

#AH tries to remove trend from each cell in raster
#Need a sequence of years the right length for each season
#Need to iterate over each seasonal matrix
#could be rolled into your loop above, when it works!
#for DJF
lm_x <- seq(1:dim(DJF)[2]) #regress over time, exact years not important
DJF_r <- t(resid(lm(t(DJF)~ lm_x)))

#amy tries a loop and fails!
for (i in unique(substring(yr_season, 6))){
  lm_x <- seq(1:dim(get(i))[2])
  r <- t(resid(lm(t(get(i)) ~ lm_x))) #something is stuffing up here,got it working on DJF alone
  assign(paste0(i), r) #and maybe here with the names vs matrices
  rm(r)
}
###

#create rasters to ingest the spatial correlations
CorT <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
Cor <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
temp <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)


## Presenting the treeCorr Function!!!
## This compresses multitudes of previous coding into a very small package
## Correlation, masking all done in one go. Plot ready after this.

## SC: The following is the combination of previous functions.
fullCorr <- function(x, y){
  rng <- range(as.numeric(substr(grep(
    unique(substr(as.character(colnames(x)), 7, 9)), 
    colnames(x), value=T),2, 5)))
  trYr <-trDat[which(trDat$year >= rng[1] & trDat$year<= rng[2]),]
  for(i in 1:dim(x)[1]){
    Cor[i] <- cor(x=x[i,], y = trYr[,y], method = 'pearson') ## create correlation based on tree ring
    CorT[i] <- cor.test(x=x[i,], y = trYr[,y], method = 'pearson')$p.value ## p values used to create the cropped confidence intervals
  }
  CorT[CorT > 0.05] <- NA
  Cor <- brick(Cor)
  temp <- mask(Cor, CorT)
  return(temp)
}

DJF_c <- fullCorr(DJF, "ars")
MAM_c <- fullCorr(MAM, "ars")
JJA_c <- fullCorr(JJA, "ars")
SON_c <- fullCorr(SON, "ars")

rm(SON, DJF, JJA, MAM, trDat, datM) 

#Stack the new rasters to be displayed.
Seasons <- stack(DJF_c[[1]], MAM_c[[1]], JJA_c[[1]], SON_c[[1]])
names(Seasons) <- c("DJF", "MAM", "JJA", "SON")

rm(Cor, CorT, temp)
#Load in a shapefile and crop for the region of interest 
library(rgdal)
library(rgeos)
coast_shapefile <- crop(readOGR("GISdata/ne_10m_coastline.shp"), ext)

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
levelplot(Seasons, layout=c(2,2), col.regions = col5, pretty=TRUE, 
          main= paste("Pearson R Reanal SFCP", F_yr, "-", L_yr),
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))


