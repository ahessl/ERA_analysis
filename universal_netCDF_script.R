
## The only potential things that should need changed for use with other netCDF files
## 1. Tree ring data file (Line 16)
## 2. Beginning and ending years for tree ring data (F_yr, L_yr)
## 3. netCDF file name and/or var.name
## 4. Spatial extent (ext)
## 6. Color ramp (col5)
## 7. Coastline file, although if packaged with this file it shouldn't need changed.
rm(list=ls())

### Installing Shawn's fancy package
library(devtools)
devtools::install_github("ahessl/ERA_analysis/SpatCor")

#setwd("C:/Users/S/Dropbox/ATSE Thesis Workspace/ERAclimateExploration")

library(ncdf4)
library(raster)
library(SpatCor)
library(dplR)

#Name the files for convenience:

netcdf.file <- "../NOAA20thc/hgt500.nc"
treering.file <- "../../KBP_South/KBPS_cull_gap.rwl_tabs.txt"

# Bring in netcdf file via SpatCor
dat <- ncdfRead(netcdf.file)


# Get a list of the variables to choose from and confirm that the time 
#origin and units are appropriate......
nc <- nc_open(netcdf.file)

#select the variable
print(names(nc[['var']]))

var.name <- names(nc[['var']])[1] #MIGHT NEED ADJUSTING


t <- ncvar_get(nc, "time")
tunits <- ncatt_get(nc, "time", "units")
print(tunits)
tustr <- strsplit(tunits$value, " ")

#extract lon dimension to check its range
lon <- ncvar_get(nc, "lon")
mlon <- max(lon)
print(mlon)

nc_close(nc)

dat <- rotate(brick(netcdf.file, varname= var.name)) #rotate converts lons from 0-360 to -180-180

#for some reason, brick cannot deal with "months since..." ???
if (unlist(tustr)[1]=="months") {
      mons <- length(t)
      org <- as.Date(unlist(tustr)[3])
      dates_f <- seq.Date(org, by="month", length.out=mons)
      names(dat) <- dates_f
    } else {
      names(dat) <- names(dat)
    }
## set spatial extent for area interested in. Longitude min - max then latitude min - max
#ext <- extent(144, 149, -44, -40) #awap micro extent
ext <- extent(0, 180, -80, -4)
#ext <- extent(-180, 180, -80, 0)
datC <- crop(dat, ext) #Spatial crop using extent

## Tree ring data
trDat_tab <- read.table(treering.file, header = TRUE)


## Decide start year and end year based on target and tree ring data
F_yr <- min(as.numeric(substr(names(datC), 2, 5)))
#F_yr <- 1959
#L_yr <- 1998
L_yr <- as.numeric(max(trDat_tab$year))

#Subset Tree Ring data to F_yr and L_Yr
trDat_raw <- subset(trDat_tab, year>=F_yr & year<=L_yr )

#Detrend Tree Ring Index Linear or FD
tr_i <- trDat_raw$ars
## First Differences Method
#r <- diff(tr_i), 1))
# Linear Model Method
# lm_x <- seq(1:length(tr_i))
# r <- resid(lm(tr_i ~ lm_x))
#Difference from spline
r <- tr_i - ffcsaps(trDat_raw$ars, trDat_raw$year, nyrs=30) #subtract detrended series from original series
trDat <- as.data.frame(cbind(year=trDat_raw$year, ars=r)) #put dataframe back together


# Subset rasters based on the F_yr, L_yr
# substr(names(dat)....) reads raster name and selects the characters (2 through 5) 
# in the names (time in this case) to extract years and compares to values in tree ring file
datC <- datC[[which(as.numeric(substr(names(dat), 2, 5)) >= F_yr & 
                      as.numeric(substr(names(dat), 2, 5)) <= L_yr)]]
datC <- datC[[-c(1:2, (nlayers(datC)-3):nlayers(datC))]] #removes first incomplete season JF and last SON from year


### Use seasNm from Shawn - data, SchulmanShift TRUE or FALSE for seasonal offset (SH), lag (0,1), function (e.g. sum, mean)
datM <- seasNm(datC, F, 0, mean)

# Subset season, replace NAs with -9999 for correlation/regression.
# Can run linear model extracting residuals or first differences.

# SC: First Differences added!

for (i in unique(substring(names(datM), 7))){ 
  d <- values(subset(datM, grep(i, names(datM), value=T)))
  d[is.na(d[])] <- -9999
  assign(paste0(i), d)
  ## First Differences Method
  #r <- t(diff(t(d), 1))
  # Linear Model Method
  lm_x <- seq(1:dim(get(i))[2])
  r <- t(resid(lm(t(get(i)) ~ lm_x)))
  rm(lm_x)
  assign(paste0(i), r)
  rm(r, d)
}

#This doesn't work yet, trying to figure it out.....Make changes easier....
#detrCL(datM, "fd")

#create rasters to store the spatial correlations
CorT <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
Cor <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
temp <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)

## Correlation, masking all done in one go. Plot ready after this.

fullCorr <- function(ClimMatr, TRInd, TrYrs){ # ClimMatr=climate data matrix, TRInd=column name from trDat in (df$Indices), TrYrs = years (df$years)
  #Creates date range with tree and climate data max/min years
  rng <- c(max(range(as.numeric(substr(colnames(ClimMatr),2, 5)))[1], range((TrYrs))[1]),
           min(range(as.numeric(substr(colnames(ClimMatr),2, 5)))[2], range((TrYrs))[2]))
  trYr <- setNames(data.frame(TrYrs, TRInd), c("year", "data")) #internal function use -- speeds things up a little more for cor tests
  trYr <-trYr[which(trYr$year >= rng[1] & trYr$year<= rng[2]),] #crops the tree ring data based on the range; rings can be variable (stupid SH lag thing)
  x <- ClimMatr[,which(substring(colnames(ClimMatr), 2, 5) >= rng[1] & substring(colnames(ClimMatr), 2, 5) <= rng[2])]
  for(i in 1:dim(x)[1]){
    Cor[i] <- cor(x=x[i,], y = trYr$data, method = 'pearson') ## create correlation based on tree ring
    CorT[i] <- cor.test(x=x[i,], y = trYr$data, method = 'pearson')$p.value ## p values used to create the cropped confidence intervals
  }
  CorT[CorT > 0.05] <- NA
  Cor <- raster::brick(Cor)
  temp <- raster::mask(Cor, CorT)
  return(temp)
}

# part of the SpatCor package
# Working on how to integrate the previous 3 things...
SON_c <- fullCorr(SON, trDat$ars, trDat$year)
DJF_c <- fullCorr(DJF, trDat$ars, trDat$year)
MAM_c <- fullCorr(MAM, trDat$ars, trDat$year)
JJA_c <- fullCorr(JJA, trDat$ars, trDat$year)


rm(SON, DJF, JJA, MAM, trDat, datM) 

#Stack the new rasters to be displayed.
Seasons <- stack(SON_c[[1]], DJF_c[[1]], MAM_c[[1]], JJA_c[[1]])
names(Seasons) <- c("SON", "DJF", "MAM", "JJA")

rm(Cor, CorT, temp)
#Load in a shapefile and crop for the region of interest 
library(rgdal)
library(rgeos)
coast_shapefile <- crop(readOGR("../GISdata/ne_10m_coastline.shp"), ext)

#Create color ramps for mapping and number of colors to use
library(colorRamps)
col5 <- colorRampPalette(c('#08519c', 'gray96', "#fee0d2", "firebrick3"))

#Load the lattice packages to display the maps
library(rasterVis)
library(gridExtra)

title.txt <- basename(netcdf.file) #b/c I am losing track of what's what;
product <- unlist(strsplit (netcdf.file, "[/]"))[2]

#Plot correlation map using the color ramp and levels using lattice
#Layout dictates how many plots in row, col format
#col.region is the color ramp previously created, pretty makes the color breaks happen logically
#colorkey is information about the legend, wanted bottom so have to give it the space
#par.settings is various graphical settings outside of plots
# +layer(...) adds the coastlines onto the map
pdf(paste("../SpatialCorr/", product, "_", title.txt, "_", F_yr, "_", L_yr, ".pdf", sep=""), width=6, height=4.3, pointsize=7, family="Helvetica")
levelplot(Seasons, layout=c(2,2), col.regions = col5, pretty=TRUE, 
          main= paste(title.txt, F_yr, "-", L_yr),
          colorkey=list(space="right"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))
dev.off()

