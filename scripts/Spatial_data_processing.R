
# spatial data 
#

#------elevation -------
#added elevation to points worked 10/15/2023 but apparently the 

library(raster)
library(rgdal)
sts<-read.csv('data_analysis/sites_coordinates.csv')
dem<-raster('data_analysis/elev/USGS_13_n44w114_20130911.tif')
coordinates<-data.frame(lon=sts$lon, lat=sts$lat)
coordinates_sp <- SpatialPoints(coordinates, proj4string = CRS(proj4string(dem)))
sts$elevation <- extract(dem, coordinates_sp)

# we add elevation to the bat data.

# bat2021_v2<-read.csv('data_analysis/bat2021_v2.csv')
bat2021_v2<-left_join(bat2021_v2, sts, by="site")

write.csv(bat2021_v2,file = 'data_analysis/bat2021_v3.csv') # we update the data base







# percentage of riparian area?
# We will use NDVI data.

ndvi_2021<-raster('data_analysis/NDVI/MYD13Q1.A2021185.h09v04.061.2021202231848.hdf',crs=)
coordinates_sp <- SpatialPoints(coordinates, proj4string = CRS(proj4string(ndvi_2021)))

buffer_distance <-100

# Create a buffer around the point
buffer <- raster::buffer(coordinates_sp, width = buffer_distance)

# Extract the NDVI values for the buffer
ndvi_values <- extract(ndvi_2021, buffer)

dvi_stats <- lapply(ndvi_values, function(x) {
  if (length(x) > 0) {
    # Calculate mean NDVI within the buffer
    mean_ndvi <- mean(x, na.rm = TRUE)
    return(mean_ndvi)
  } else {
    # If no NDVI values were found within the buffer, return NA
    return(NA)
  }
})



