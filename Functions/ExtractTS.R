#ExtractTS
#Extract a time series from a grid box of netcdf data

#Extract a single cell 
gbox <- extent (140, 150, -42, -40)
gbox_m <- extract(datC, grid_box, fun="mean", method="simple")


