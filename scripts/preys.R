
# generate preys raster

rm(list=ls())
source('scripts/utils.R')

base <- raster('spatial_data/env2016/tmean.asc')
base[] <- 1:ncell(base) # unique values

files <- list.files('spatial_data/preys', recursive=T, full.names = T, pattern='.shp')

preys <- base
preys[] <- 0

for (file in files) {
  shp <- st_read(file) %>% st_set_crs(25830) %>% st_transform(crs(base))
  values <- na.omit(unlist(raster::extract(base, shp))) # extract pixels in each cell
  binary <- calc(base, function(x) { ifelse(x %in% values, 1, 0) }) # generate binary raster
  preys <- preys + binary # accumulate binaries in the raster
}

writeRaster(preys, 'spatial_data/preys/preys.asc', format='ascii', overwrite=T)
