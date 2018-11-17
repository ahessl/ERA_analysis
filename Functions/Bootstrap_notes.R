library(boot)
 
   x <- DJF[1,]
   y <- trYr$ars
   
     cor(cbind(x, y))

    myCor <- function(data, index){
       cor(data[index, ])[1, 2]
       }
   
  b <- boot(data=cbind(x, y), statistic=myCor, R=1000)
  boot.ci(b,type="basic")$basic[4:5]
  