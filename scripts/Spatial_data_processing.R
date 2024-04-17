
# spatial data 
#

#------elevation -------
#added elevation to points worked 10/15/2023 but apparently the 

library(raster)
library(rgdal)
library(rgeos)# creates buffers

sts<-read.csv('data_for_analysis/sites_coordinates.csv')
dem<-raster('data_for_analysis/elev/USGS_13_n44w114_20130911.tif')
coordinates<-data.frame(lon=sts$lon, lat=sts$lat)
coordinates_sp <- SpatialPoints(coordinates, proj4string = CRS(proj4string(dem)))
sts$elevation <- extract(dem, coordinates_sp)

# we add elevation to the bat data.

# bat2021_v2<-read.csv('data_analysis/bat2021_v2.csv')
bat2021_v2<-left_join(bat2021_v2, sts, by="site")

write.csv(sts,file = 'data_for_analysis/elev/elev.csv') # write the elevation







# percentage of riparian area?
# We will use NDVI data.
# dowloaded data from :https://search.earthdata.nasa.gov/search/granules?p=C2307290656-LPCLOUD!C2307290656-LPCLOUD&pg[1][v]=t&pg[1][qt]=2021-06-01T00%3A00%3A00.000Z%2C&pg[1][gsk]=start_date&pg[1][m]=download&pg[1][cd]=f&q=NDVI&sb[0]=-115.68164%2C42.80419%2C-112.83398%2C45.33567&as[science_keywords][0]=Biosphere%3AVegetation%3AVegetation%20Index%3ANormalized%20Difference%20Vegetation%20Index%20(Ndvi)&tl=1680503012!3!!&fst0=land%20surface&fst1=Biosphere&fsm1=Vegetation&fs11=Vegetation%20Index&fs21=Normalized%20Difference%20Vegetation%20Index%20(Ndvi)&lat=20.834681452138085&long=-171.7310116952234

ndvi_2021<-raster('data_for_analysis/NDVI/MYD13Q1.A2021185.h09v04.061.2021202231848.hdf') # Modis vegetation index product

#check projection
pj<-crs(ndvi_2021)
# field sites re-projected to raster file size.

site.pj <- SpatialPoints(sts[-1], proj4string = crs(ndvi_2021))
# Extract longitude and latitude coordinates from site.pj
lon <- coordinates(site.pj)[, 1]
lat <- coordinates(site.pj)[, 2]

# Create a new SpatialPoints object with corrected coordinates order
site_corrected <- SpatialPoints(coords = cbind(lon, lat), proj4string = CRS(proj4string(site.pj)))

# Create a buffer around the points with a radius of 100 meters
buffer <- buffer(site_corrected, width = 100)

# Extract NDVI values within the buffer area
ndvi_within_buffer <- extract(ndvi_2021, buffer) # this part is not working. I don't know why.

# Calculate mean NDVI value within the buffer
mean_ndvi <- sapply(ndvi_within_buffer, mean, na.rm = TRUE)

# Print mean NDVI values
print(mean_ndvi)


plot(ndvi_2021, col = rainbow(100))
legend("topright", legend = "MODISVegetation Index")
title(main = "Vegetation Index Map")
points(coordinates_sp$x, coordinates_sp$y, col = "red", pch = 20)




