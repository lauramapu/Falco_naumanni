
# buildings database

rm(list=ls())
#source('scripts/utils.R')

library(sf)
library(raster)
library(terra)
library(dplyr)
library(mapSpain)
library(MazamaSpatialUtils)

getpeninsularspain <- function() {
  provinces <- mapSpain::esp_get_prov()[!mapSpain::esp_get_prov()$iso2.prov.name.es %in%
                                          c("Las Palmas", "Santa Cruz de Tenerife", "Baleares", "Melilla", "Ceuta"), ]
  provinces$campo <- 1
  spain <- provinces %>%
    dissolve(field='campo') %>%
    st_transform(crs=4326)
}

# crop+mask function for all variables
mask_raster <- function(raster_layer) {
  maskk <- st_transform(spain, st_crs(raster_layer))
  cropped <- raster::crop(raster_layer, extent(maskk))
  masked <- raster::mask(cropped, maskk)
  return(masked)
}

spain <- getpeninsularspain()

lustre_dir <- '/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/'

# load files in terra format (absolute paths)
eubucco <- vect(paste0(lustre_dir, 'spatial_data/buildings/v0_1-ESP.gpkg')) %>%
  crop(project(vect(spain), "EPSG:3035"))

base <- raster(paste0(lustre_dir, 'spatial_data/climate/present/CHELSA_tas_01_2016_V.2.1.tif')) %>%
  mask_raster() %>%
  rast() %>%
  project("EPSG:3035")

# rasterize polygons
fraction <- terra::rasterize(eubucco, base, background=0, touches=F, cover=T)
# with cover=T it estimates the fraction of cell that is covered by polygons
count <- terra::rasterize(eubucco, base, background=0, touches=T, fun='count')
age <- terra::rasterize(eubucco, base, background=0, touches=T, field='age', fun='mean')

writeRaster(raster(fraction), paste0(lustre_dir, 'spatial_data/buildings/fraction.asc'), format='ascii', overwrite=T)
writeRaster(raster(count), paste0(lustre_dir, 'spatial_data/buildings/count.asc'), format='ascii', overwrite=T)
writeRaster(raster(age), paste0(lustre_dir, 'spatial_data/buildings/age.asc'), format='ascii', overwrite=T)
