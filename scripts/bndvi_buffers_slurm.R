# calculate buffer for BNDVI
# values for each year

rm(list=ls())

library(sf)
library(raster)
library(dplyr)
library(mapSpain)
library(MazamaSpatialUtils)
library(doParallel)

# base <- raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/climate/present/CHELSA_tas_01_2016_V.2.1.tif')
# 
# getpeninsularspain <- function() {
#   provinces <- mapSpain::esp_get_prov()[!mapSpain::esp_get_prov()$iso2.prov.name.es %in%
#                                           c("Las Palmas", "Santa Cruz de Tenerife", "Baleares", "Melilla", "Ceuta"), ]
#   provinces$campo <- 1
#   spain <- provinces %>%
#     dissolve(field='campo') %>%
#     st_transform(crs=4326)
# }
# 
# spain <- getpeninsularspain()
# 
# merge_rasters <- function(pattern) {
#   paths <- list.files('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi', pattern=pattern, full.names=T)
#   list <- lapply(paths, raster)
#   merged <- do.call(merge, list)
#   return(merged)
# }
# 
# # crop+mask function for all variables
# mask_raster <- function(raster_layer) {
#   maskk <- st_transform(spain, st_crs(raster_layer))
#   cropped <- raster::crop(raster_layer, extent(maskk))
#   masked <- raster::mask(cropped, maskk)
#   return(masked)
# }
# 
# years <- c('Mean_2016', 'StdDev_2016', 'Mean_2017', 'StdDev_2017', 'Mean_2018', 'StdDev_2018')
# 
# for (i in years) {
#   options(scipen = 999) # turn off scientific notation
#   r <- merge_rasters(i)
#   r@data@names<-i
#   print(paste0(i, ' merged'))
#   r <- mask_raster(r)
#   r <- projectRaster(r, base)
#   writeRaster(r, paste0('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_', i, '.asc'))
#   rm(r)
#   print(paste0(i, ' saved'))
# }

# again calculate median of shannon index in small buffer
# (calculated in script clc_preparation)
landuses <- read.csv('data/landuses_12km.csv')
median_shannon <- median(landuses[,'shannon_index']) # 0.6836354

# load last modeling data (clc_preparation)
mod_data <- read.csv('data/modeling_data_v3.csv')
# we erase -frozen, -tmin, -tmax, -elev because of high corr
# convert to sf
dataset_sf <- st_as_sf(mod_data, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

# detect cores
if (Sys.getenv("SLURM_JOB_ID") != "") {
  n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK"))  # SLURM
} else {
  n.cores <- parallel::detectCores() - 1 # local
}
# create cluster
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

# calculate buffers per year and shannon index

bndvi <- foreach(i=1:nrow(dataset_sf),
                 .verbose=T,
                 .combine='rbind',
                 .packages=c('sf','raster'),
                 .inorder=T) %dopar% {

                   if (landuses[i,'shannon_index']>=median_shannon) {

                     buffer <- st_buffer(dataset_sf[i,], dist = 2000)

                     if (dataset_sf$Ano_CS[i]==2016) {
                       stack2016 <- stack(raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_Mean_2016.asc'),
                                          raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_StdDev_2016.asc'))
                       values <- na.omit(as.data.frame(raster::extract(stack2016, buffer)))
                       rm(stack2016)
                     }
                     else if (dataset_sf$Ano_CS[i]==2017) {
                       stack2017 <- stack(list(raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_Mean_2017.asc'),
                                               raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_StdDev_2017.asc')))
                       values <- na.omit(as.data.frame(raster::extract(stack2017, buffer)))
                       rm(stack2017)
                     }
                     else if (dataset_sf$Ano_CS[i]==2018) {
                       stack2018 <- stack(list(raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_Mean_2018.asc'),
                                               raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_StdDev_2018.asc')))
                       values <- na.omit(as.data.frame(raster::extract(stack2018, buffer)))
                       rm(stack2018)
                     }

                     mean_cont <- colMeans(values)

                     return(mean_cont)
                   }

                   else {

                     buffer <- st_buffer(dataset_sf[i,], dist = 4500)

                     if (dataset_sf$Ano_CS[i]==2016) {
                       stack2016 <- stack(list(raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_Mean_2016.asc'),
                                               raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_StdDev_2016.asc')))
                       values <- na.omit(as.data.frame(raster::extract(stack2016, buffer)))
                       rm(stack2016)
                     }
                     else if (dataset_sf$Ano_CS[i]==2017) {
                       stack2017 <- stack(list(raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_Mean_2017.asc'),
                                               raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_StdDev_2017.asc')))
                       values <- na.omit(as.data.frame(raster::extract(stack2017, buffer)))
                       rm(stack2017)
                     }
                     else if (dataset_sf$Ano_CS[i]==2018) {
                       stack2018 <- stack(list(raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_Mean_2018.asc'),
                                               raster('/mnt/lustre/scratch/nlsas/home/csic/byc/abl/Falco_naumanni/spatial_data/bndvi/bndvi_StdDev_2018.asc')))
                       values <- na.omit(as.data.frame(raster::extract(stack2018, buffer)))
                       rm(stack2018)
                     }

                     mean_cont <- colMeans(values)

                     return(mean_cont)
                   }
                 }

bndvi <- as.data.frame(bndvi)
colnames(bndvi) <- c('Mean_BNDVI', 'StdDev_BNDVI')
mod_data <- cbind(mod_data, bndvi)

write.csv(mod_data, 'data/modeling_data_v4.csv', row.names=F)
