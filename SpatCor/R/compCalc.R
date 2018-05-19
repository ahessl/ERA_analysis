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
    crMn <- function(fullMn, qyrs){
      tr.c <- fullMn[[which(as.numeric(substr(names(fullMn), 2, 5)) %in% qyrs) ]]
      dat.c <- raster::stackApply(tr.c, substring(names(tr.c), 7), mean) #seasonal mean for wide years
      names(dat.c) <- unique(substring(names(tr.c), 7)) #meaningful names
      dat.c <- raster::subset(dat.c, c("SON", "DJF", "MAM", "JJA")) #keeps order correct
      com.c <- dat.c - datMs
    }
    sumStatb <- function(uyrs, lyrs, fllMn){
      t1 <- list(u=list(), l=list(), m=list())
      for(i in unique(substr(names(fllMn), 7,9))){
        t1$u[[i]]<-paste(length(names(fllMn[[which(as.numeric(substr(names(fllMn), 2, 5)) %in% uq_yrs & 
                                                    substr(names(fllMn),7,9) %in% i)]])), i, sep=" ")
        t1$l[[i]]<-paste(length(names(fllMn[[which(as.numeric(substr(names(fllMn), 2, 5)) %in% lq_yrs & 
                                                    substr(names(fllMn),7,9) %in% i)]])), i, sep=" ")
        t1$m[[i]]<-paste(length(names(fllMn[[which(substr(names(fllMn),7,9) %in% i)]])), i, sep=" ")
      }
      
      Stats<- 
        list(Upper = uq_yrs, Lower = lq_yrs,Mean = paste(min(as.numeric(substr(names(fllMn), 2, 5))), "-",
                          max(as.numeric(substr(names(fllMn), 2, 5) )), sep = " "),
             UpperSsn = paste(t1$u), LowerSsn = paste(t1$l), MeanSsn = paste(t1$m)
        )
      return(Stats)
    }
    sumStat <- function(qyrs, fllMn){
      t1 <- list(q=list(), m=list())
      for(i in unique(substr(names(fllMn), 7,9))){
        t1$q[[i]]<-paste(length(names(fllMn[[which(as.numeric(substr(names(fllMn), 2, 5)) %in% qyrs & 
                                                     substr(names(fllMn),7,9) %in% i)]])), i, sep=" ")
        t1$m[[i]]<-paste(length(names(fllMn[[which(substr(names(fllMn),7,9) %in% i)]])), i, sep=" ")
      }
      
      Stats<- 
        list(QuantileYrs = qyrs, Mean = paste(min(as.numeric(substr(names(fllMn), 2, 5))), "-",
                                                         max(as.numeric(substr(names(fllMn), 2, 5) )), sep = " "),
              QuantSsn = paste(t1$q), MeanSsn = paste(t1$m)
        )
      return(Stats)
    }
    
    if(up == "b"){
      qs <- quantile
      quant <- quantile(trData, probs = qs)
      uq_yrs <- trDat$year[which(trData > quant[2])]
      lq_yrs <- trDat$year[which(trData < quant[1])]
      a <- menu(c("Yes", "No"), title = "Would you to store composite summary?")
      if (a == 1){
        #bth <- list(uq_yrs, lq_yrs)
        #mxl <- max(sapply(bth, length))
        #Quantiles <<- setNames(data.frame(sapply(bth, function(x){c(x, rep(NA, mxl - length(x)))})), 
                              #c("Upper.Quant", "Lower.Quant"))
        ComSummary <<- sumStatb(uq_yrs, lq_yrs, datM) 
        com.h <<- crMn(fullMean, qyrs = uq_yrs)
        com.l <<- crMn(fullMean, qyrs = lq_yrs)
      }else{
        com.h <<- crMn(fullMean, uq_yrs)
        com.l <<- crMn(fullMean, lq_yrs)
      }
    }else{
      quants <- quantile(trData, probs = quantile)
      if(up == "u"){
        q_yrs <- trDat$year[which(trData > quants)]
        a <- menu(c("Yes", "No"), title = "Would you to store the composite summary?")
        if (a == 1){
          #UpperQuant <<- data.frame(Upper.Quant = q_yrs)
          ComSummary <<- sumStat(q_yrs, datM)
          com.h <<- crMn(fullMean, q_yrs)
        }else{
          com.h <<- crMn(fullMean, q_yrs)
        }
      }else{
        if(up == "l"){
          q_yrs <- trDat$year[which(trData < quants)]
          a <- menu(c("Yes", "No"), title = "Would you to store the composite summary?")
          if (a == 1){
            #LowerQuant <<- data.frame(Upper.Quant = q_yrs)
            ComSummary <<- sumStat(q_yrs, datM)
            com.l <<- crMn(fullMean, q_yrs)
          }else{
            com.l <<- crMn(fullMean, q_yrs)
          }
        }else{
        stop("Please use l for lower, u for upper, or b for both")
        }
      }
    }
  }
}
  