
# buildings database

rm(list=ls())
source('scripts/utils.R')

spain <- getpeninsularspain()

# load files in terra format
eubucco <- vect('spatial_data/buildings/v0_1-ESP.gpkg')

# # try to add 1 column tu sum
# eu_sf <- st_read(eubucco)
# eu_sf$num <- 1
# eubucco <- vect(eu_sf)
# count <- terra::rasterize(eubucco, base, background=0, touches=F, field='num', fun='sum')
# writeRaster(raster(count), 'spatial_data/buildings/count.asc', format='ascii', overwrite=T)

base <- raster('spatial_data/climate/present/CHELSA_tas_01_2016_V.2.1.tif') %>%
  mask_raster() %>%
  rast() %>%
  project("EPSG:3035")

# rasterize polygons calculating different variables: surface, count and mean age

fraction <- terra::rasterize(eubucco, base, background=0, touches=F, cover=T) %>%
  mask_raster()
# with cover=T it estimates the fraction of cell that is covered by polygons
writeRaster(raster(fraction), 'spatial_data/buildings/fraction.asc', format='ascii', overwrite=T)

count <- terra::rasterize(eubucco, base, background=0, touches=F, fun='count', sum=T) %>%
  mask_raster()
# need to add sum=T because otherwise it creates a binary raster
writeRaster(raster(count), 'spatial_data/buildings/count.asc', format='ascii', overwrite=T)

eubucco_f <- eubucco[!is.na(eubucco$age), ]
age <- terra::rasterize(eubucco_f, base, background=0, touches=F, field='age', fun='mean') %>%
  mask_raster()
writeRaster(raster(age), 'spatial_data/buildings/age.asc', format='ascii', overwrite=T)

fraction <- raster('spatial_data/buildings/fraction.asc')
crs(fraction) <- crs(base)
fraction <- mask_raster(fraction)
writeRaster(fraction, 'spatial_data/buildings/fraction.asc', format='ascii', overwrite=T)

# calculate fragmentation index between buildings with landscapemetrics package

# clean env
rm(list=ls())
source('scripts/utils.R')

# add new variable (for now just fraction) to modeling df
# we must calculate buffer

# load dataset
dataset <- read.csv2('data/modeling_data_v3.csv')[,1:5]
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

spain <- getpeninsularspain()

clc_reclass <- raster('spatial_data/clc_2018/clc_reclass.asc')
fraction <- raster('spatial_data/buildings/fraction.asc')
clc_urban <- calc(clc_reclass, function(x) {ifelse(x==1, 1, 0)})

for (i in 1:nrow(dataset_sf)) {
  buffer <- st_buffer(dataset_sf[i,], dist = 2000)
  values <- na.omit(as.data.frame(raster::extract(stack2016, buffer)))
}

