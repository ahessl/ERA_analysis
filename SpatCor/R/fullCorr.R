#' Spatial Field Correlation Function
#'
#' This function allows you to create spatial field correlations with tree rings and
#' gridded climate data
#' @param x climate data
#' @param y tree ring data
#' @return results from temp
#' @export

## Correlation, masking all done in one go. Plot ready after this.

fullCorr <- function(x, y){ # x=climate data, y=column name from trDat in ""
  rng <- range(as.numeric(substr(grep( #Creates date range based on the climate data
    unique(substr(as.character(colnames(x)), 7, 9)),
    colnames(x), value=T),2, 5)))
  trYr <- setNames(data.frame(trDat$year, y), c("year", "data")) #internal function use -- speeds things up a little more for cor tests
  trYr <-trYr[which(trYr$year >= rng[1] & trYr$year<= rng[2]),] #crops the tree ring data based on the range; rings can be variable (stupid SH lag thing)
  for(i in 1:dim(x)[1]){
    Cor[i] <- cor(x=x[i,], y = trYr$data, method = 'pearson') ## create correlation based on tree ring
    CorT[i] <- cor.test(x=x[i,], y = trYr$data, method = 'pearson')$p.value ## p values used to create the cropped confidence intervals
  }
  CorT[CorT > 0.05] <- NA
  Cor <- raster::brick(Cor)
  temp <- raster::mask(Cor, CorT)
  return(temp)
}
