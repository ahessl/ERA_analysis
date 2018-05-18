#' Composite Correlation
#'
#' This function allows you to create composite spatial field correlations 
#' @param trData tree ring data - probably will appear something like trDat$...
#' @param quantile selection of quantile numerically
#' @param fullMean the full dataset of seasonal means from year 1..n
#' @param up deciding whether upper or lower quantile for calculations, default at up
#' @return com.c based on choice in quantile.
#' @export
compCalc <- function(trData, quantile, fullMean, up = "u"){
  if(is.numeric(quantile) == FALSE ){
    stop("The quantile argument must be a numeric value from 0 to 1")
  }else{
    
    datMs <- raster::stackApply(fullMean,substring(names(fullMean), 7), mean) #create seasonal mean of entire dataset
    names(datMs) <- unique(substring(names(fullMean), 7)) #give them meaningful names
    datMs <- raster::subset(datMs, c("SON", "DJF", "MAM", "JJA")) #reorder because they are setup as MAM, JJA, SON, DJF
      
    quants <- quantile(trData, probs = quantile)
    if(up == "up"){
      uq_yrs <- trDat$year[which(trData > quants)]
      a <- devtools::menu(c("Yes", "No"), "Would you to store the upper quantile data?")
      if (a == 1){
        assign(UpperYears, uq_yrs, enviro=.GlobalEnv)
      }
      tr.c <- fullMean[[which(as.numeric(substr(names(fullMean), 2, 5)) %in% uq_yrs) ]]
      
      dat.c <- raster::stackApply(tr.c, substring(names(tr.c), 7), mean) #seasonal mean for wide years
      names(dat.c) <- unique(substring(names(tr.c), 7)) #meaningful names
      dat.c <- raster::subset(dat.c, c("SON", "DJF", "MAM", "JJA")) #keeps order correct
      com.c <- dat.c - datMs #composite difference
      return(com.c)
    }else{
      uq_yrs <- trDat$year[which(trData < quants)]
      
      tr.c <- fullMean[[which(as.numeric(substr(names(fullMean), 2, 5)) %in% uq_yrs) ]]
      
      dat.c <- raster::stackApply(tr.c, substring(names(tr.c), 7), mean) #seasonal mean for wide years
      names(dat.c) <- unique(substring(names(tr.c), 7)) #meaningful names
      dat.c <- raster::subset(dat.c, c("SON", "DJF", "MAM", "JJA")) #keeps order correct
      com.c <- dat.c - datMs #composite difference
      return(com.c)
    }
    }
  }
  