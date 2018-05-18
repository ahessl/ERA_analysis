#' Composite Correlation
#'
#' This function allows you to create composite spatial field correlations 
#' @param trData tree ring data - probably will appear something like trDat$...
#' @param quantile selection of quantile numerically single value or list of c()
#' @param fullMean the full dataset of seasonal means from year 1..n
#' @param up deciding whether upper or lower quantile for calculations, default not set.
#' @return com.c based on choice in quantile.
#' @export
compCalc <- function(trData, quantile, fullMean, up){
  if(is.numeric(quantile) == FALSE ){
    stop("The quantile argument must be a numeric value from 0 to 1")
  }else{
    
    datMs <- raster::stackApply(fullMean,substring(names(fullMean), 7), mean) #create seasonal mean of entire dataset
    names(datMs) <- unique(substring(names(fullMean), 7)) #give them meaningful names
    datMs <- raster::subset(datMs, c("SON", "DJF", "MAM", "JJA")) #reorder because they are setup as MAM, JJA, SON, DJF
    crMn <- function(fullMean, q_yrs){
      tr.c <- fullMean[[which(as.numeric(substr(names(fullMean), 2, 5)) %in% q_yrs) ]]
      dat.c <- raster::stackApply(tr.c, substring(names(tr.c), 7), mean) #seasonal mean for wide years
      names(dat.c) <- unique(substring(names(tr.c), 7)) #meaningful names
      dat.c <- raster::subset(dat.c, c("SON", "DJF", "MAM", "JJA")) #keeps order correct
      com.c <- dat.c - datMs
    }
    quants <- quantile(trData, probs = quantile)
    if(up == "b"){
      uq_yrs <- trDat$year[which(trData > quants[2])]
      lq_yrs <- trDat$year[which(trData < quants[1])]
      a <- menu(c("Yes", "No"), "Would you to store the quantile years?")
      if (a == 1){
        bth <- list(uq_yrs, lq_yrs)
        mxl <- max(sapply(bth, length))
        Quantiles <<- setNames(data.frame(sapply(bth, function(x){c(x, rep(NA, mxl - length(x)))})), 
                              c("Upper.Quant", "Lower.Quant"))
        tr.u <<- crMn(fullMean, uq_yrs)
        tr.l <<- crMn(fullMean, lq_yrs)
        return()
      }else{
        tr.u <<- crMn(fullMean, uq_yrs)
        tr.l <<- crMn(fullMean, lq_yrs)
      }
      if(up == "u"){
        q_yrs <- trDat$year[which(trData > quants)]
        a <- menu(c("Yes", "No"), "Would you to store the upper quantile years?")
        if (a == 1){
        UpperQuart <<- data.frame(Upper.Quart = q_yrs)
        com.h <<- crMn(fullMean, q_yrs)
        return()
        }else{
          if(up == "l"){
            q_yrs <- trDat$year[which(trData < quants)]
            a <- menu(c("Yes", "No"), "Would you to store the lower quantile years?")
            if (a == 1){
            LowerQuart <<- data.frame(Upper.Quart = q_yrs)
            crMn(fullMean, q_yrs)
            com.l <<- crMn(fullMean, q_yrs)
            return()
        }
      }
      }
    }
    }
  }
}
  