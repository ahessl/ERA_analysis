#setwd("C:/Users/S/Dropbox/ATSE Thesis Workspace/ERA_Download")
rm(list=ls())
library(raster)
library(ncdf4)

# Get a list of the stupid variables to choose from and confirm that the time 
#origin and units are appropriate......
nc <- nc_open('../ERA_Download/era_interim_moda_2mT.nc')

#select the variable
print(names(nc[['var']]))
var.name <- names(nc[['var']])

t <- ncvar_get(nc, "time")
tunits <- ncatt_get(nc, "time", "units")
print(tunits)

nc_close(nc)


dat <- brick("../ERA_Download/era_interim_moda_2mT.nc", varname=var.name) #adjusted for selection above
trDat <- read.table("../../KBP_South/KBPS_cull_gap.rwl_tabs.txt", header = TRUE)

## Select start year and end year based on target and tree ring data
F_yr <- min(as.numeric(substr(names(dat), 2, 5)))
L_yr <- as.numeric(max(trDat$year))

#Crop data to it
trDat <-trDat[which(trDat$year >= F_yr-1 & trDat$year<= L_yr),]

## set spatial extent for area interested in. Longitude min - max then latitude min - max
## Creating this way allows for other created rasters to recognize as an extent with 
## the proper spatial extent, otherwise have to set xmin, xmax, ymin, ymax individually.

ext <- extent(60, 200.25, -80.25, -4.50)

#Spatial crop using extent
datC <- crop(dat, ext)
datC <- datC[[which(as.numeric(substr(names(datC), 2, 5)) >= F_yr & 
                      as.numeric(substr(names(datC), 2, 5)) <= L_yr)]]
datC <- datC[[-c(1:2, (nlayers(datC)-3):nlayers(datC))]] #removes first incomplete season JF and last SON from year


##### Seasonal Indexing ######
library(chron)

yr_mo_dy <- substr(names(datC), 2, 11)
d <- as.Date(gsub(".", '/', yr_mo_dy, fixed = T)) #fix the format by replacing "." with "/"
## changed $year==12 to $mon<8 (POSIXlt indexes months from 0
## and growing season starts in SEP) 


yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                      as.POSIXlt( d )$year - 
                      1*(as.POSIXlt( d )$mon<8) ,   # offset needed for grwoing season in SH
                    c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                      1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                    , sep="-")

## Get yearly seasonal means
datM <- stackApply(datC, yr_season, mean) #raster with mean for each season
names(datM) <- unique(yr_season)

########### Creating Composite ############

#### Create seasonal mean of all climate data ####

## SC: This needs cleaned a bit more.
datMs <- stackApply(datM,substring(names(datM), 7), mean) #create seasonal mean of entire dataset
names(datMs) <- unique(substring(names(datM), 7)) #give them meaningful names
datMs <- subset(datMs, c(3,4,1,2)) #reorder because they are setup as MAM, JJA, SON, DJF


#### Create seasonal mean of high years and low years (whatever that means for your purposes) ####

### This section uses climate extremes to extract specific years
### created c.comp.csv with columns STR w/ values h and l for high and low years; column year is the year it occurred
### combined high and low in one place
# c.cm <- read.csv("C:/Users/S/Dropbox/ATSE Thesis Workspace/ERA_Download/c.comp.csv")
# 
# c.yr <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in%  c.cm$year[c.cm$STR == "h"] )]] #High value year
# #c.yr <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in%  c.cm$year[c.cm$STR == "l"] )]] #Low value year
# 
# dat.c <- stackApply(c.yr, substring(names(c.yr), 7), mean)#seasonal mean for high years
# names(dat.c) <- unique(substring(names(c.yr), 7)) #meaningful names
# com.c <- dat.c - datMs #composite difference high years

#dat.c <- stackApply(l.yr, substring(names(c.yr), 7), mean) #seasonal mean for low years
#names(dat.c) <- unique(substring(names(c.yr), 7)) #meaningful names
#com.c <- dat.c - datMs #composite difference low years


#### Tree ring guiding the analysis....for whatever it's worth. ####

## There has been very little work done to this. 
#Better to select based on quantiles:
## Smallest x% years, Largest x% years
quants <- quantile(trDat$ars, probs = c(0.10, 0.90))
lq_yrs <- trDat$year[which(trDat$ars<quants[1])]
uq_yrs <- trDat$year[which(trDat$ars>quants[2])]

## Extract wide and narrow years from TR indices. h for high growth; l for low growth
#tr.h <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in%  trDat$year[order(trDat$mr_kbp)[1:5]])]]
#tr.l <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in%  trDat$year[order(trDat$mr_kbp, decreasing = T)[1:5]])]]
tr.h <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in% uq_yrs) ]]
tr.l <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in% lq_yrs) ]]

dat.h <- stackApply(tr.h, substring(names(tr.h), 7), mean) #seasonal mean for wide years
names(dat.h) <- unique(substring(names(tr.h), 7)) #meaningful names

dat.l <- stackApply(tr.l, substring(names(tr.l), 7), mean)
names(dat.l) <- unique(substring(names(tr.l), 7))

com.h <- dat.h - datMs #composite difference
com.l <- dat.l - datMs #composite difference

########### Plotting the composites ###########

library(rgdal)
library(rgeos)
coast_shapefile <- crop(readOGR("../GISData/ne_10m_coastline.shp"), ext)

# Create color ramps for mapping and number of colors to use
library(colorRamps)
col5 <- colorRampPalette(c('#08519c', 'lightblue3','gray96', "#fee0d2", "firebrick3"))

# Load the lattice packages to display the maps
library(rasterVis)
library(gridExtra)

# Make a plot!
levelplot(com.h, layout=c(2,2), col.regions = col5, pretty=TRUE, main= paste("Sea Level Pressure Composite:", F_yr, "-", L_yr),
          colorkey=list(space="right"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))

#### Correlation with tree ring data ####

# Pull High Value Years from tree ring data using c.cm from before
tr.yr <- trDat[which(trDat$year %in%  c.cm$year[c.cm$STR == "h"] ),]


# Low Value Years from tree ring data using c.cm from before
#tr.yr <- trDat[which(trDat$year %in%  c.cm$year[c.cm$STR == "l"] ),]

for (i in unique(substring(yr_season, 6))){
  d <- values(subset(c.yr, grep(i, names(c.yr), value=T)))
  d[is.na(d[])] <- -9999
  r <- d
  #  assign(paste0(i), d)
  ## First Differences Method
#  r <- t(diff(t(d), 1)) ## Don't use this right now. Have to figure out how to deal with it.
  ## Linear Model Method  
  #lm_x <- seq(1:dim(get(i))[2])
  #r <- t(resid(lm(t(get(i)) ~ lm_x))) 
  #rm(lm_x)
  assign(paste0(i), r)
  rm(r, d)
}

#create rasters to store the spatial correlations
CorT <- setExtent(raster(nrow = nrow(c.yr), ncol = ncol(c.yr)),ext)
Cor <- setExtent(raster(nrow = nrow(c.yr), ncol = ncol(c.yr)),ext)
temp <- setExtent(raster(nrow = nrow(c.yr), ncol = ncol(c.yr)),ext)

## Correlation, masking all done in one go. Plot ready after this.

fullCorr <- function(x, y){ # x=climate data, y=column name from trDat in ""
  trYr <-tr.yr
  for(i in 1:dim(x)[1]){
    Cor[i] <- cor(x=x[i,], y = trYr[,y], method = 'pearson') ## create correlation based on tree ring
    CorT[i] <- cor.test(x=x[i,], y = trYr[,y], method = 'pearson')$p.value ## p values used to create the cropped confidence intervals
  }
  CorT[CorT > 0.05] <- NA
  Cor <- brick(Cor)
  temp <- mask(Cor, CorT)
  return(temp)
}

DJF_c <- fullCorr(DJF, "mr_kbp")
MAM_c <- fullCorr(MAM, "mr_kbp")
JJA_c <- fullCorr(JJA, "mr_kbp")
SON_c <- fullCorr(SON, "mr_kbp")

rm(SON, DJF, JJA, MAM, trDat, datM) 

#Stack the new rasters to be displayed.
Seasons <- stack(SON_c[[1]], DJF_c[[1]], MAM_c[[1]], JJA_c[[1]])
names(Seasons) <- c("SON", "DJF", "MAM", "JJA")

rm(Cor, CorT, temp)

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