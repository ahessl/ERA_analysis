
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
treering.file <- "../../Law_Dome/DSS_Oct-March_average_d18O_2012-1900 CE_20180530.csv"

# Bring in netcdf file via SpatCor
dat <- ncdfRead(netcdf.file)




#for some reason, ncdfRead is not able to read in files that index time in months 
# if (unlist(tustr)[1]=="months") {
#       mons <- length(t)
#       org <- as.Date(unlist(tustr)[3])
#       dates_f <- seq.Date(org, by="month", length.out=mons)
#       names(dat) <- dates_f
#     } else {
#       names(dat) <- names(dat)
#     }

#set spatial extent for area interested in. Longitude min - max then latitude min - max
ext <- extent(-180, 180, -80, 0) #use the full extent to be sure that rotate is working
datC <- crop(dat, ext) #Spatial crop using extent

## Tree ring data
trDat_tab <- read.csv(treering.file, header = TRUE, skip=3)


## Decide start year and end year based on target and tree ring data
yrs_rng <- range(intersect(as.numeric(substr(names(datC), 2, 5)), 
                           as.numeric(trDat_tab$year)))
#F_yr <- yrs_rng[1]
L_yr <- yrs_rng[2]
#OR manually set
F_yr <- 1959
#L_yr <- 1998

#Subset Tree Ring data to F_yr and L_Yr
trDat_raw <- subset(trDat_tab, year>=F_yr & year<=L_yr )

#Detrend Tree Ring Index Linear, FD or Spline
tr_i <- trDat_raw$d18O_ONDJFM
## First Differences Method
#r <- diff(tr_i), 1))
# Linear Model Method
# lm_x <- seq(1:length(tr_i))
# r <- resid(lm(tr_i ~ lm_x))
#Difference from spline
r <- tr_i - ffcsaps(trDat_raw$d18O_ONDJFM, trDat_raw$year, nyrs=30) #subtract detrended series from original series
trDat <- as.data.frame(cbind(year=trDat_raw$year, d18O_ONDJFM=r)) #put dataframe back together


# Subset rasters based on the F_yr, L_yr
# substr(names(dat)....) reads raster name and selects the characters (2 through 5) 
# in the names (time in this case) to extract years and compares to values in tree ring file
datC <- datC[[which(as.numeric(substr(names(dat), 2, 5)) >= F_yr & 
                      as.numeric(substr(names(dat), 2, 5)) <= L_yr)]]
datC <- datC[[-c(1:2, (nlayers(datC)-3):nlayers(datC))]] #removes first incomplete season JF and last SON from year

seasNm <- function(climDat, SchulmanShift = FALSE, lg = 0, FUN){
  yr_mo_dy <- substr(names(climDat), 2, 11)
  d <- as.Date(gsub(".", '/', yr_mo_dy, fixed = T)) #fix the format by replacing "." with "/"
  if(SchulmanShift == TRUE) {
    if(lg == 0) {
      ### Current Growing Year (as determined by tree dates)
      yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                            as.POSIXlt( d )$year - 
                            1*(as.POSIXlt( d )$mon<8) ,   # offset needed for growing season in SH
                          c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                            1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                          , sep="-")
      datM <- stackApply(climDat, yr_season, match.fun(FUN)) #raster with mean for each season
      names(datM) <- unique(yr_season) 
      return(datM)
    } else {
      
      ### One Year Climate Lag (e.g. tree year 1980 will be associated with climate data for tree year 1979)
      yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                            as.POSIXlt( d )$year + 
                            1*(as.POSIXlt( d )$mon>7) ,   # offset needed for lagged season in SH
                          c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                            1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                          , sep="-")
      datM <- stackApply(climDat, yr_season, match.fun(FUN)) #raster with mean for each season
      names(datM) <- unique(yr_season) 
    }
  }else{
### No Schulman shift - December as year change; December included with Jan and Feb
    if(lg ==0 ){
      yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                            as.POSIXlt( d )$year - 
                            1*(as.POSIXlt( d )$mon<2) ,
                          c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                            1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                          , sep="-")
      datM <- stackApply(climDat, yr_season, match.fun(FUN))
      names(datM) <- unique(yr_season)
      return(datM)
    }else{
      ##1 year lag
      yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                            as.POSIXlt( d )$year - 
                            1*(as.POSIXlt( d )$mon<2) ,
                          c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                            1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                          , sep="-")
      datM <- stackApply(climDat, yr_season, match.fun(FUN))
      names(datM) <- unique(yr_season)
      return(datM)
    }
  }
}

### Use seasNm from Shawn - data, SchulmanShift TRUE or FALSE for seasonal offset (SH), lag (0,1), function (e.g. sum, mean)
datM <- seasNm(datC, F, 0, mean)


###### TEST CODE - DO NOT USE YET ######

## SC: This code has 3 menus: What month the year ends; How many months included; and apply Schulman
## SC: I have tested it with my data - needs looked at by other eyes. Will also require cleaning.

library(dplyr)
yr_mo_dy <- substr(names(datC), 2, 11)
d <- as.Date(gsub(".", '/', yr_mo_dy, fixed = T)) #fix the format by replacing "." with "/"

mo <- menu(month.name, title = "What month does the year end?")

season <- switch(menu(c(2,3,4,6), title = "How many months in each season?"), 2,3,4,6)
s.s <- menu(c("No", "Yes"), title = "Should a Schulman shift be applied?") - 1

mo <- ifelse (mo == 12, mo - 1, mo)

# Dataframe of first letter of the month names, season (chosen above), and the POSIXlt value for 
# season name generation and indexing later
indx <- data.frame(cbind
                   (MonName = substr(months(seq.Date(as.Date(paste("1999", mo+1, "01", sep = "/")),, "month", length = 12)), 1,1),
                     ssn = rep(1:(12/season), each = season),
                     POSmon = as.POSIXlt(seq.Date(as.Date(paste("1999", mo+1, "01", sep = "/")),, "month", 
                                                  length = 12))$mon))

# Only way I could figure out how to generate a list of seasonal names that could change based on previous parameters.
ssn.nms <- indx %>% group_by(ssn) %>% summarize (ssnNm = paste(MonName, collapse = ""))

yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                      as.POSIXlt( d )$year - 
                      s.s*(as.POSIXlt( d )$mon < mo) ,   # offset applied with Schulman shift question
                    ssn.nms[[2]][          # indexing from created dataframe and name list
                      ifelse(a<-match(as.POSIXlt(d)$mon, indx$POSmon), indx$ssn[a], NA)
                      ]
                    , sep="-")

datM <- stackApply(datC, yr_season, mean) #raster with mean for each season
names(datM) <- unique(yr_season)

###### END TEST CODE ######


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

