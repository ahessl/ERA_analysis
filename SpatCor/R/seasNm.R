#' Dating of netcdf Files and creating a seasonal compilation
#'
#' This function allows you to create spatial field correlations with tree rings and
#' gridded climate data
#' @param climDat climate data
#' @param hem hemisphere of tree ring data
#' @param lg lag of seasonal effect
#' @param FUN function to be used in stackApply
#' @return results from datM
#' @export

##climDat should be pretreated, Hemisphere is for north/south naming ("n", "s"), lg is a lag (0, 1)

seasNm <- function(climDat, hem, lg, FUN){
    yr_mo_dy <- substr(names(climDat), 2, 11)
    d <- as.Date(gsub(".", '/', yr_mo_dy, fixed = T)) #fix the format by replacing "." with "/"
    if(hem == "s") {
      if(lg == 0) {
    ### Current Growing Year (as determined by tree dates)
      yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                      as.POSIXlt( d )$year - 
                      1*(as.POSIXlt( d )$mon<8) ,   # offset needed for growing season in SH
                    c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                      1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                    , sep="-")
      datM <- stackApply(climDat, yr_season, match.fun(FUN)) #raster with mean for each season
      names(datM) <- unique(yr_season) 
      return(datM)
      } else {

### One Year Climate Lag (e.g. tree year 1980 will be associated with climate data for tree year 1979)
        yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                      as.POSIXlt( d )$year + 
                      1*(as.POSIXlt( d )$mon>7) ,   # offset needed for lagged season in SH
                    c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                      1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                    , sep="-")
        datM <- stackApply(climDat, yr_season, match.fun(FUN)) #raster with mean for each season
        names(datM) <- unique(yr_season) 
      }
      }else{
        if(lg ==0 ){
          yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                              as.POSIXlt( d )$year + 
                              1*(as.POSIXlt( d )$mon>10) ,   # offset needed for growing season in SH
                            c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                              1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                            , sep="-")
        datM <- stackApply(climDat, yr_season, match.fun(FUN))
        return(datM)
        }else{
          yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                                as.POSIXlt( d )$year + 
                                1*(as.POSIXlt( d )$mon>10) ,   # offset needed for growing season in SH
                              c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                                1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                              , sep="-")
          datM <- stackApply(climDat, yr_season, match.fun(FUN))
          return(datM)
        }
      }
    }
