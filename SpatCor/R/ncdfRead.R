#' Reading in netCDF Files
#'
#' This allows users to choose climate variables interactively and creates a raster brick
#' @param ncdf tree ring data - probably will appear something like trDat$...
#' @return dat based on choice in quantile.
#' @export

ncdfRead <- function(ncdf){
  nc <- ncdf4::nc_open(ncdf)
  
  #select the variable
  idx <- devtools::menu(names(nc[['var']]), "What variable would you like?")
  var.name <<- names(nc[['var']])[idx]
  
  t <- ncdf4::ncvar_get(nc, "time")
  tunits <- ncdf4::ncatt_get(nc, "time", "units")
  tustr <- strsplit(tunits$value, " ")[[1]][1]
  
  lon <- ncvar_get(nc, "lon")
  mlon <- max(lon)
  
  if (mlon < 181){
    dat <- raster::brick(ncdf, varname = var.name)
  }else{
    dat <- raster::rotate(raster::brick(ncdf, varname=var.name)) #adjusted for selection above
  }
  ncdf4::nc_close(nc)
  
  #for some reason, brick cannot deal with "months since..." ???
  if (tustr=="months") {
    mons <- length(t)
    org <- as.Date(unlist(tustr)[3])
    dates_f <- seq.Date(org, by="month", length.out=mons)
    names(dat) <- dates_f
  } else {
    names(dat) <- names(dat)
  }
  return(dat)
}
