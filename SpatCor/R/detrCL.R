#' Dating of netcdf Files and creating a seasonal compilation
#'
#' This function allows you to create spatial field correlations with tree rings and
#' gridded climate data
#' @param climData climate data
#' @param mod mode of detrending
#' @return results from climData
#' @export
#' 

detrCL <- function(climData, mod){ #mod is fd or lm
  for (i in unique(substring(names(climData), 7))){
    d <- values(subset(climData, grep(i, names(climData), value=T)))
    d[is.na(d[])] <- -9999
    if(mod == "fd"){
    ## First Differences Method
    r <- t(diff(t(d), 1))
    do.call("<<-", list(i, r))
    }else{
      if(mod == "lm"){
    # Linear Model Method
        lm_x <- seq(1:dim(get(i))[2])
        r <- t(resid(lm(t(get(i)) ~ lm_x)))
        do.call("<<-", list(i, r))
      }else{
        do.call("<<-", list(i, d))
      }
    }
  }
}
    