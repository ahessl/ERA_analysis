setwd("C:/Users/S/Dropbox/ATSE Thesis Workspace/ERA_Download")
rm(list=ls())
library(raster)
library(ncdf4)

# Get a list of the stupid variables to choose from......
nc = nc_open('era_interim_moda_All.nc')
var.name = names(nc[['var']])
nc_close(nc)


dat <- brick("era_interim_moda_All.nc", varname=var.name[7])
trDat <- read.csv("../ERAclimateExploration/chronos.csv", head=T)

## Decide start year and end year based on target and tree ring data
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
### created csv with columns STR w/ values h and l for high and low years; column year is the year it occurred
### combined high and low in one place
c.cm <- read.csv("c.comp.csv")
h.yr <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in%  c.cm$year[c.cm$STR == "h"] )]]
l.yr <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in%  c.cm$year[c.cm$STR == "l"] )]]

dat.h <- stackApply(h.yr, substring(names(h.yr), 7), mean)#seasonal mean for high years
names(dat.h) <- unique(substring(names(h.yr), 7)) #meaningful names

dat.l <- stackApply(l.yr, substring(names(l.yr), 7), mean) #seasonal mean for high years
names(dat.l) <- unique(substring(names(l.yr), 7)) #meaningful names

com.h <- dat.h - datMs #composite difference
com.l <- dat.l - datMs #composite difference

rm(c.cm, h.yr, l.yr, dat.h, dat.l, datM, datC)
#### Tree ring guiding the analysis....for whatever it's worth. ####

## Extract wide and narrow years from TR indices. h for high growth; l for low growth
#tr.h <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in%  trDat$year[order(trDat$mr_kbp)[1:5]])]]
#tr.l <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in%  trDat$year[order(trDat$mr_kbp, decreasing = T)[1:5]])]]

#dat.h <- stackApply(tr.h, substring(names(tr.h), 7), mean) #seasonal mean for wide years
#names(dat.w) <- unique(substring(names(tr.h), 7)) #meaningful names

#dat.l <- stackApply(tr.l, substring(names(tr.l), 7), mean)
#names(dat.n) <- unique(substring(names(tr.l), 7))

#com.h <- dat.h - datMs #composite difference
#com.l <- dat.l - datMs #composite difference

########### Plotting the composites ###########

library(rgdal)
library(rgeos)
coast_shapefile <- crop(readOGR("../GISData/ne_10m_coastline.shp"), ext)

#Create color ramps for mapping and number of colors to use
library(colorRamps)
col5 <- colorRampPalette(c('#08519c', 'lightblue3','gray96', "#fee0d2", "firebrick3"))

#Load the lattice packages to display the maps
library(rasterVis)
library(gridExtra)


#### Low value years ####
levelplot(com.l, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Narrow Composite Mean SLP Diff 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))


#### High value years ####
levelplot(com.h, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Wide Composite Mean SLP Diff 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))

#### Low minus High ####
com.lh <- com.l - com.h

levelplot(com.lh, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Narrow-Wide Composite Mean SLP Diff 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))
