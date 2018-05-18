#setwd("C:/Users/S/Dropbox/ATSE Thesis Workspace/ERA_Download")
rm(list=ls())

### Installing Shawn's fancy package
library(devtools)
devtools::install_github("ahessl/ERA_analysis/SpatCor")

library(raster)
library(ncdf4)
library(SpatCor)

netcdf.file <- "../NOAA20thc/hgt500.nc"
treering.file <- "../../KBP_South/KBPS_cull_gap.rwl_tabs.txt"

## New function....gets rid of the messing with the file here.
## Also has interactive functionality!!!
#dat <- ncdfRead(netcdf.file) #not rotatiing!!!


# Get a list of the stupid variables to choose from and confirm that the time 
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

dat <- rotate(brick(netcdf.file, varname=var.name)) #adjusted for selection above
#rotate converts lons from 0-360 to -180-180



#for some reason, brick cannot deal with "months since..." ???
#  if (unlist(tustr)[1]=="months") {
#      mons <- length(t)
#      org <- as.Date(unlist(tustr)[3])
#      dates_f <- seq.Date(org, by="month", length.out=mons)
#      names(dat) <- dates_f
#    } else {
#      names(dat) <- names(dat)
#    }

trDat <- read.table(treering.file, header = TRUE)

## Select start year and end year based on target and tree ring data
F_yr <- min(as.numeric(substr(names(dat), 2, 5)))
L_yr <- as.numeric(max(trDat$year))

#Crop data to it
trDat <-trDat[which(trDat$year >= F_yr-1 & trDat$year<= L_yr),]

## set spatial extent for area interested in. Longitude min - max then latitude min - max
## Creating this way allows for other created rasters to recognize as an extent with 
## the proper spatial extent, otherwise have to set xmin, xmax, ymin, ymax individually.

#ext <- extent(60, 200.25, -80.25, -4.50)
ext <- extent(-180, 180, -80, 0)

#Spatial crop using extent
datC <- crop(dat, ext)

#Temporal limits using fyr, lyr
datC <- datC[[which(as.numeric(substr(names(datC), 2, 5)) >= F_yr & 
                      as.numeric(substr(names(datC), 2, 5)) <= L_yr)]]
datC <- datC[[-c(1:2, (nlayers(datC)-3):nlayers(datC))]] #removes first incomplete season JF and last SON from year

#library(SpatCor)
datM <- seasNm(datC, "s", 0, mean)
##### Seasonal Indexing ######
#library(chron)

#yr_mo_dy <- substr(names(datC), 2, 11)
#d <- as.Date(gsub(".", '/', yr_mo_dy, fixed = T)) #fix the format by replacing "." with "/"
## changed $year==12 to $mon<8 (POSIXlt indexes months from 0
## and growing season starts in SEP) 


#yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                      #as.POSIXlt( d )$year - 
                      #1*(as.POSIXlt( d )$mon<8) ,   # offset needed for grwoing season in SH
                    #c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                     # 1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                    #, sep="-")

## Get yearly seasonal means
#datM <- stackApply(datC, yr_season, mean) 
#names(datM) <- unique(yr_season)

##Need to detrend data####

########### Creating Composite ############

#### Create seasonal mean of all climate data ####

## Added in some SpatCor fanciness to get things done a little faster with fewer commands.
## comCalc has 4 arguments tree ring data, numeric quantile, data, 
## and up (default is upper) - for math purposes; replace up with any value and you get lower quantile.
## 
## Error Message:
## If quantile is not a number, stops the whole thing
## Need to add this error message
## If the seasons aren't even, stops the process and provides the data to determine the season.
##
## Reproduce error codes: 
##    com.h <- compCalc(trDat$mr_ars, "upper", datM)

com.h <- compCalc(trDat$ars, .85, datM)
com.l <- compCalc(trDat$ars, .15, datM, 5)

## SC: This needs cleaned a bit more.
#datMs <- stackApply(datM,substring(names(datM), 7), mean) #create seasonal mean of entire dataset
#names(datMs) <- unique(substring(names(datM), 7)) #give them meaningful names
#datMs <- subset(datMs, c(3,4,1,2)) #reorder because they are setup as MAM, JJA, SON, DJF


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


#### Select years to composite based on tree rings.
# based on quantiles:
## Smallest x% years, Largest x% years

#quants <- quantile(trDat$ars, probs = c(0.15, 0.85))
#lq_yrs <- trDat$year[which(trDat$ars<quants[1])]
#uq_yrs <- trDat$year[which(trDat$ars>quants[2])]


## Extract wide and narrow years from TR indices. h for high growth; l for low growth
#tr.h <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in% uq_yrs) ]]
#tr.l <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in% lq_yrs) ]]

#dat.h <- stackApply(tr.h, substring(names(tr.h), 7), mean) #seasonal mean for wide years
#names(dat.h) <- unique(substring(names(tr.h), 7)) #meaningful names
#dat.h <- subset(dat.h, c("SON", "DJF", "MAM", "JJA")) #keeps order correct

#dat.l <- stackApply(tr.l, substring(names(tr.l), 7), mean)
#names(dat.l) <- unique(substring(names(tr.l), 7))
#dat.l <- subset(dat.l, c("SON", "DJF", "MAM", "JJA"))


#com.h <- dat.h - datMs #composite difference to get anomalies
#com.l <- dat.l - datMs #composite difference to get anomalies

########### Plotting the composites ###########

library(rgdal)
coast_shapefile <- crop(readOGR("../GISData/ne_10m_coastline.shp"), ext)

# Create color ramps for mapping and number of colors to use
library(colorRamps)
col5 <- colorRampPalette(c('#08519c', 'lightblue3','gray96', "#fee0d2", "firebrick3"))

# Load the lattice packages to display the maps
library(rasterVis)

title.txt <- basename(netcdf.file) #b/c I am losing track of what's what;
product <- unlist(strsplit (netcdf.file, "[/]"))[2]

compNms <-c("com.h", "com.l")
for(i in compNms){
  rwdth <- toupper(substring(i, 5))
  pdf(paste(product, "_", title.txt, "_", F_yr, "_", L_yr, rwdth,".pdf", sep=""),
      width=6, height=4.3, pointsize=7, family="Helvetica")
  print(levelplot(get(i), layout=c(2,2), col.regions = col5, pretty=TRUE, main= paste(title.txt, F_yr, "-", L_yr),
                  colorkey=list(space="right"),
                  par.settings = list(layout.heights=list(xlab.key.padding=1),
                                      strip.background=list(col="lightgrey")), par.strip.text = list(font="bold")) + 
          layer(sp.lines(coast_shapefile)))
  dev.off()
}

# Make two plots and save them in directory called "Composites"
#Can you make a loop to do this?  I struggled and gave up.
#pdf(paste("../Composites/", product, "_", title.txt, "_", F_yr, "_", L_yr, "L", ".pdf", sep=""), width=6, height=4.3, pointsize=7, family="Helvetica")
# levelplot(com.l, layout=c(2,2), col.regions = col5, pretty=TRUE, main= paste(title.txt, F_yr, "-", L_yr),
#          colorkey=list(space="right"),
#          par.settings = list(layout.heights=list(xlab.key.padding=1),
#                              strip.background=list(col="lightgrey")
#          ), par.strip.text = list(font="bold")) + 
#  layer(sp.lines(coast_shapefile))
# dev.off()
 
# pdf(paste("../Composites/", product, "_", title.txt, "_", F_yr, "_", L_yr, "H",".pdf", sep=""), width=6, height=4.3, pointsize=7, family="Helvetica")
# levelplot(as.name(try[1]), layout=c(2,2), col.regions = col5, pretty=TRUE, main= paste(title.txt, F_yr, "-", L_yr),
#           colorkey=list(space="right"),
#           par.settings = list(layout.heights=list(xlab.key.padding=1),
#                               strip.background=list(col="lightgrey")
#           ), par.strip.text = list(font="bold")) + 
#   layer(sp.lines(coast_shapefile))
# dev.off()

################################################################### 
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