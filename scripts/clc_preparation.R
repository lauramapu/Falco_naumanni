
# CLC 10m for land uses to simplify our models because these bichos are not specialists

rm(list=ls())
source('scripts/utils.R')

# load dataset
dataset <- read.csv2('data/modeling_data.csv')[,1:5]
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

spain <- getpeninsularspain()

mask_raster <- function(x) {
  reproj <- st_transform(spain, crs=crs(x))
  cropp <- crop(x, extent(reproj))
  maskk <- mask(cropp, reproj)
  return(maskk)
}

clc <- raster('spatial_data/clc_2018/U2018_CLC2018_V2020_20u1.tif') %>%
  mask_raster()
saveRDS(clc, 'objects/clc_masked.rds')
# clc <- readRDS('objects/clc_masked.rds')

legend <- data.frame(
  ID = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
         31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 48, NA),
  Description = c(
    "Continuous urban fabric", "Discontinuous urban fabric", "Industrial or commercial units", "Road and rail networks and associated land",
    "Port areas", "Airports", "Mineral extraction sites", "Dump sites", "Construction sites", "Green urban areas", 
    "Sport and leisure facilities", "Non-irrigated arable land", "Permanently irrigated land", "Rice fields", 
    "Vineyards", "Fruit trees and berry plantations", "Olive groves", "Pastures", "Annual crops associated with permanent crops", 
    "Complex cultivation patterns", "Land principally occupied by agriculture with significant areas of natural vegetation", 
    "Agro-forestry areas", "Broad-leaved forest", "Coniferous forest", "Mixed forest", "Natural grasslands", 
    "Moors and heathland", "Sclerophyllous vegetation", "Transitional woodland-shrub", "Beaches - dunes - sands", 
    "Bare rocks", "Sparsely vegetated areas", "Burnt areas", "Glaciers and perpetual snow", "Inland marshes", 
    "Peat bogs", "Salt marshes", "Salines", "Intertidal flats", "Water courses", "Water bodies", "Coastal lagoons", 
    "Estuaries", "Sea and ocean", "NODATA", NA
  )
)
legend

legend$new <- c('urban', 'urban', 'urban', 'urban', 'urban', 'urban', 'urban', 'urban', 'urban', 'urban',
                'urban', 'dry_cropland', 'wet_cropland', 'wet_cropland', 'perm_cropland', 'perm_cropland',
                'perm_cropland', 'grassland', 'mixed_cropland', 'mixed_cropland', 'cropland_natural', 
                'agro_forestry', 'forest', 'forest', 'forest', 'grassland', 'shrubland', 'shrubland', 
                'shrubland', 'bare_soil', 'bare_soil', 'sparse_vegetation', 'burnt', 'frozen', 'wetland',
                'wetland', 'wetland', 'salines', 'salines', 'water', 'water', 'water', 'water', 'water', 'ND', NA
                )

# replace categories with numbers so that we can assign to raster
mapping <- setNames(seq_along(unique(legend$new)), unique(legend$new))
mapping[19] <- NA
legend$new_num <- as.numeric(mapping[legend$new])

write.csv(legend, 'spatial_data/clc_2018/legend.csv', row.names=F)
# legend <- read.csv('spatial_data/clc_2018/legend.csv')

# reclass rasterlayer with new categories
clc_reclass <- reclassify(clc, as.matrix(legend[,c('ID','new_num')]))
writeRaster(clc_reclass, 'spatial_data/clc_2018/clc_reclass.asc', format='ascii', overwrite=T)
# clc_reclass <- raster('spatial_data/clc_2018/clc_reclass.asc')

# detect cores
if (Sys.getenv("SLURM_JOB_ID") != "") {
  n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK"))  # SLURM
} else {
  n.cores <- parallel::detectCores() - 1 # local
}
# create cluster
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

landuses <- foreach(i = 1:nrow(dataset),
                    .verbose=T,
                    .combine='rbind',
                    .packages=c('sf','raster'),
                    .inorder=T) %dopar% {
  
  buffer <- st_buffer(dataset_sf[i,], dist=20000) %>%
    st_transform(crs=crs(clc_reclass))
  
  # extract values and calculate proportions
  values <- na.omit(unlist(raster::extract(clc_reclass, buffer)))
  proportions <- table(factor(values, levels = unique(legend$new_num))) / sum(values)
  
  # shannon index
  shannon <- -sum(proportions[proportions > 0] * log(proportions[proportions > 0]))
  
  # return result
  return(c(proportions, shannon))
}

colnames(landuses) <- c(unique(legend$new)[-19], 'shannon_index')
write.csv(landuses, 'data/landuses_12km.csv', row.names=F)
# landuses <- read.csv('data/landuses_12km.csv')

median_shannon <- median(landuses[,'shannon_index']) # 0.6836354

 # load env stacks

env2016 <- readRDS("objects/env2016.rds")
env2017 <- readRDS("objects/env2017.rds")
env2018 <- readRDS("objects/env2018.rds")

columnnames <- c('tmean', 'tmax', 'tmin', 'tseas', 'pr', 'prseas', 'elev', 'slope', 'linear_infrastr', 'streams', 'preys',
                 'forest', 'water', 'shrub', 'bare', 'for_shr_bare', 'settle_med', 'settle_low', 'for_shr_grass',
                 'for_shr_crop', 'for_shr_agric', 'mosaic_low', 'settle_high', 'crop_med', 'grass_low', 
                 'crop_high', 'grass_med', 'grass_high', 'mosaic_med', 'ext_perm_crop', 'crop_low', 'int_perm_crop', 
                 'mosaic_high')
names(env2016) <- columnnames
names(env2017) <- columnnames
names(env2018) <- columnnames

# convert to stack
stack2016 <- stack(env2016[1:11])
stack2017 <- stack(env2017[1:11])
stack2018 <- stack(env2018[1:11])

# calculate all the other values based on shannon index
# high index = small buffer, low index = large buffer
envshannon <- foreach(i=1:nrow(dataset_sf),
                      .verbose=T,
                      .combine='rbind',
                      .packages=c('sf','raster'),
                      .inorder=T) %dopar% {
  
  if (landuses[i,'shannon_index']>=median_shannon) {
    
    buffer <- st_buffer(dataset_sf[i,], dist = 2000)
    
    if (dataset_sf$Ano_CS[i]==2016) {
      values <- na.omit(as.data.frame(raster::extract(stack2016, buffer)))
    }
    else if (dataset_sf$Ano_CS[i]==2017) {
      values <- na.omit(as.data.frame(raster::extract(stack2017, buffer)))
    }
    else if (dataset_sf$Ano_CS[i]==2018) {
      values <- na.omit(as.data.frame(raster::extract(stack2018, buffer)))
    }
    
    mean_cont <- colMeans(values) 
    
    return(as.numeric(c(mean_cont, landuses[i,])))
  }
  
  else {
  
  buffer <- st_buffer(dataset_sf[i,], dist = 4500)
  
  if (dataset_sf$Ano_CS[i]==2016) {
    values <- na.omit(as.data.frame(raster::extract(stack2016, buffer)))
  }
  else if (dataset_sf$Ano_CS[i]==2017) {
    values <- na.omit(as.data.frame(raster::extract(stack2017, buffer)))
  }
  else if (dataset_sf$Ano_CS[i]==2018) {
    values <- na.omit(as.data.frame(raster::extract(stack2018, buffer)))
  }
  
  mean_cont <- colMeans(values) 
  
  # recalculate landuses also with new buffer
  values <- na.omit(unlist(raster::extract(clc_reclass, buffer)))
  proportions <- table(factor(values, levels = unique(legend$new_num))) / sum(values)
  
  return(as.numeric(c(mean_cont, proportions, landuses[i,'shannon_index'])))
  }
}

colnames(envshannon) <- c(columnnames[1:11], unique(legend$new)[-19], 'shannon_index')

finaldata <- cbind(dataset, envshannon) %>%
  dplyr::select(-ND)

# select only env cols
corrdata <- finaldata[,6:34] # same variables but with different buffer size
zero_cols <- sapply(corrdata, function(col) all(col == 0))
corrdata <- corrdata[, !zero_cols]
# Calculate the correlation matrix
cor_matrix <- cor(corrdata, use = "complete.obs")
# to df
cor_df <- as.data.frame(as.table(cor_matrix))
# filter pairs with cor higher than 0.80
high_corr <- subset(cor_df, abs(Freq) > 0.75 & abs(Freq) < 1)
# order per the magnitude of the correlation
high_corr <- high_corr[order(-abs(high_corr$Freq)), ]
# print
print(high_corr)

# we erase frozen because all zeros
# tmin, tmax and elev because of high cor

finaldata <- finaldata %>%
  dplyr::select(-frozen, -tmin, -tmax, -elev)

write.csv(finaldata, 'data/modeling_data_v3.csv', row.names=F)
