#MonthsSince(netcdf)
#Function to get time units for netcdf files, check if units
#are in 'months since', create a dates_f that is a months 
#since origin sequence of values
#requires 
library(ncdf4)
MonthsSince <- function(netcdf.file) {
    nc <- nc_open(netcdf.file)
    t <- ncvar_get(nc, "time")
    tunits <- ncatt_get(nc, "time", "units")
    print(tunits)
    tustr <- strsplit(tunits$value, " ")
    nc_close(nc)
        if (unlist(tustr)[1]=="months") {
          mons <- length(t)
          org <- as.Date(unlist(tustr)[3])
          dates_f <- seq.Date(org, by="month", length.out=mons)
          names(dat) <- dates_f
        } else {
          names(dat) <- names(dat)
        }
}