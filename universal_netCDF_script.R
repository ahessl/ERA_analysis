
## The only potential things that should need changed for use with other netCDF files
## 1. Tree ring data file (Line 16)
## 2. Beginning and ending years for tree ring data (F_yr, L_yr)
## 3. netCDF file name and/or var.name
## 4. Spatial extent (ext)
## 6. Color ramp (col5)
## 7. Coastline file, although if packaged with this file it shouldn't need changed.
rm(list=ls())
#setwd("C:/Users/S/Dropbox/ATSE Thesis Workspace/ERAclimateExploration")

library(ncdf4)
library(raster)

#Name the files for convenience:
netcdf.file <- "../CRUT/cru_ts4.01.1901.2016.tmp.dat.nc"
treering.file <- "../../KBP_South/KBPS_cull_gap.rwl_tabs.txt"

# If the netcdf file has one layer, varname isn't required. 
# If it has more than one, the first will be loaded and give names for all. Use varname="" and 
# put the variable name abbreviation in "" to load it.

# Get a list of the variables to choose from and confirm that the time 
#origin and units are appropriate......
nc <- nc_open(netcdf.file)

#select the variable
print(names(nc[['var']]))
var.name <- names(nc[['var']])[1]

t <- ncvar_get(nc, "time")
tunits <- ncatt_get(nc, "time", "units")
print(tunits)
tustr <- strsplit(tunits$value, " ")

nc_close(nc)

dat <- brick(netcdf.file, varname= var.name)

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
## Creating this way allows for other created rasters to recognize as an extent with 
## the proper spatial extent, otherwise have to set xmin, xmax, ymin, ymax individually.

#ext <- extent(60, 180, -80, -4)
ext <- extent(-180, 180, -80, 0)

#Spatial crop using extent
datC <- crop(dat, ext)

## Tree ring data
trDat <- read.table(treering.file, header = TRUE)


## Decide start year and end year based on target and tree ring data
#F_yr <- min(as.numeric(substr(names(datC), 2, 5)))
F_yr <- 1959
#L_yr <- 1998
L_yr <- as.numeric(max(trDat$year))


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

### Current Growing Year (as determined by tree dates)
yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                      as.POSIXlt( d )$year - 
                      1*(as.POSIXlt( d )$mon<8) ,   # offset needed for growing season in SH
                    c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                      1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                    , sep="-")

### One Year Climate Lag (e.g. tree year 1980 will be associated with climate data for tree year 1979)
#yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
#                      as.POSIXlt( d )$year + 
#                      1*(as.POSIXlt( d )$mon>7) ,   # offset needed for lagged season in SH
#                    c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
#                      1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
#                    , sep="-")

##remove unnecessary files
rm(yr_mo_dy, d)

# ##Get Mean of seasonal sea level pressure using seasons
datM <- stackApply(datC, yr_season, mean) #raster with mean for each season
names(datM) <- unique(yr_season) #Is this more efficient since it doesn't have to call on itself?

# Subset season, replace NAs with -9999 for correlation/regression.
# Can run linear model extracting residuals or first differences.

# SC: First Differences added!

for (i in unique(substring(yr_season, 6))){
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

#create rasters to store the spatial correlations
CorT <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
Cor <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)
temp <- setExtent(raster(nrow = nrow(datM), ncol = ncol(datM)),ext)

## Correlation, masking all done in one go. Plot ready after this.

fullCorr <- function(x, y){ # x=climate data, y=column name from trDat in ""
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
#add wide or narrow!

#Plot correlation map using the color ramp and levels using lattice
#Layout dictates how many plots in row, col format
#col.region is the color ramp previously created, pretty makes the color breaks happen logically
#colorkey is information about the legend, wanted bottom so have to give it the space
#par.settings is various graphical settings outside of plots
# +layer(...) adds the coastlines onto the map
levelplot(Seasons, layout=c(2,2), col.regions = col5, pretty=TRUE, 
          main= paste(title.txt, F_yr, "-", L_yr),
          colorkey=list(space="right"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))


