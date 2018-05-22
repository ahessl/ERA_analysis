#' Spatial Field Correlation Function
#'
#' This function allows you to create spatial field correlations with tree rings and
#' gridded climate data
#' @param x climate data
#' @param y tree ring data
#' @return results from temp
#' @export

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
