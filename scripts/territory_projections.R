# we need to generate predictions to all the territory
# to do this, we need env layers for all peninsular spain
# we're doing it for all raster centroids the same way we did for observations

rm(list=ls())
source('scripts/utils.R')

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

env.list <- readRDS("objects/env2018.rds")

columnnames <- c('tmean', 'tmax', 'tmin', 'tseas', 'pr', 'prseas', 'elev', 'slope', 'linear_infrastr', 'streams', 'preys',
                 'forest', 'water', 'shrub', 'bare', 'for_shr_bare', 'settle_med', 'settle_low', 'for_shr_grass',
                 'for_shr_crop', 'for_shr_agric', 'mosaic_low', 'settle_high', 'crop_med', 'grass_low', 
                 'crop_high', 'grass_med', 'grass_high', 'mosaic_med', 'ext_perm_crop', 'crop_low', 'int_perm_crop', 
                 'mosaic_high')

# we need to calculate same buffers we made for the modeling data for every pixel centroid in the territory
# then, basing on the smallest buffer, we calculate shannon diversity index
# with that value, we can choose one buffer or another for each point

# we need 9 dfs in the process:
# 3 distances (zero, small buffer, large buffer)
# for each period/scenario (present, future one, future two)

# extract env centroids with data as df
env.stack <- stack(env.list)
env.df <- na.omit(as.data.frame(env.stack, xy=T))
colnames(env.df)[3:35] <- columnnames
# # erase correlated and zero sdtdev variables
# env.df <- dplyr::select(env.df, -tmax, -tmin, -elev, -grass_high, -mosaic_high)
# # erase as well from stack
# env.stack <- dropLayer(env.stack, c(2, 3, 7, 27, 32))

# generate points from previous df
centroids <- env.df %>%
  st_as_sf(coords=c('x','y'), crs=crs(env.list[[1]])) %>%
  dplyr::select(geometry) # only geometry

# calculate values of the small buffer and shannon index
small.buffer <- foreach(i = 1:nrow(centroids),
                        .packages=c("sf", "dplyr"),
                        .combine='rbind',
                        .verbose=F) %dopar% {
  
                      # generate buffer
                      buffer <- st_buffer(centroids[i,], dist = 2000)
                      # extract values in buffer
                      values <- na.omit(as.data.frame(raster::extract(env.stack, buffer)))
                      
                      # calculate mean for continuous and proportion for categorical
                      mean_cont <- colMeans(values[,1:11]) # continuous variables
                      
                      landuses <- colSums(values[,12:33]) # presence of each land use
                      total <- sum(landuses)
                      mean_cat <- landuses / total # proportion of ones vs total
                      # calculate shannon diversity index (H') for land uses
                      # max H' is log(22) = 3.091042
                      shannon <- -sum(mean_cat[mean_cat > 0] * log(mean_cat[mean_cat > 0]))
                      # bind
                      mean_values <- c(mean_cont, mean_cat, shannon)
                      
                      return(mean_values)
}

colnames(small.buffer) <- c(columnnames, 'shannon_index')

# extract shannon median value from the modeling points
median.shannon <- median(read.csv2('data/fnaumanni_envdata.csv')[,'shannon_index'])

# save
write.csv2(small.buffer, 'data/env2018_smallbuffer.csv', row.names=F)
# small.buffer <- read.csv2('data/env2018_smallbuffer.csv')

# assign small buffer to an obs if H'>median, else calculate large buffer and assign

shannon.df <- foreach(i = 1:nrow(centroids),
                      .packages=c("sf", "dplyr"),
                      .combine='rbind',
                      .verbose=F) %dopar% {
  
                if (small.buffer[i, 'shannon_index'] >= median.shannon) { # assign small buffer
                  return(small.buffer[i,])
                }
                else { # calculate large buffer and assign
                  
                  # generate buffer
                  buffer <- st_buffer(centroids[i,], dist = 4500)
                  # extract values in buffer
                  values <- na.omit(as.data.frame(raster::extract(env.stack, buffer)))
                  
                  # calculate mean for continuous and proportion for categorical
                  mean_cont <- colMeans(values[,1:11]) # continuous variables
                  
                  landuses <- colSums(values[,12:33]) # presence of each land sue
                  total <- sum(landuses)
                  mean_cat <- landuses / total # proportion of ones vs total
                  
                  # bind
                  mean_values <- c(mean_cont, mean_cat, small.buffer[i, 'shannon_index'])
                  
                  return(mean_values)
                }
}

colnames(shannon.df) <- c(columnnames, 'shannon_index')

# save
write.csv2(shannon.df, 'data/env2018_buffered.csv', row.names=F)

# same with both futures

env.list <- readRDS("objects/envf1.rds")
columnnames <- c('tmean', 'tmax', 'tmin', 'tseas', 'pr', 'prseas', 'elev', 'slope', 'linear_infrastr', 'streams',
                 'forest', 'water', 'shrub', 'bare', 'for_shr_bare', 'settle_med', 'settle_low', 'for_shr_grass',
                 'for_shr_crop', 'for_shr_agric', 'mosaic_low', 'settle_high', 'crop_med', 'grass_low', 
                 'crop_high', 'grass_med', 'grass_high', 'mosaic_med', 'ext_perm_crop', 'crop_low', 'int_perm_crop', 
                 'mosaic_high')
env.stack <- stack(env.list)

shannon.fut1 <- foreach(i = 1:nrow(centroids),
                        .packages=c("sf", "dplyr"),
                        .combine='rbind',
                        .verbose=F) %dopar% {
                          
  if (small.buffer[i, 'shannon_index'] >= median.shannon) { # assign small buffer
    
    # generate buffer
    buffer <- st_buffer(centroids[i,], dist = 2000)
    # extract values in buffer
    values <- na.omit(as.data.frame(raster::extract(env.stack, buffer)))
    
    # calculate mean for continuous
    mean_cont <- colMeans(values[,1:6])
    
    return(mean_cont)
    
  }
  else { # calculate large buffer and assign
    
    # generate buffer
    buffer <- st_buffer(centroids[i,], dist = 4500)
    # extract values in buffer
    values <- na.omit(as.data.frame(raster::extract(env.stack, buffer)))
    
    # calculate mean for continuous 
    mean_cont <- colMeans(values[,1:6])
    
    return(mean_cont)
  }                        
}

colnames(shannon.fut1) <- columnnames[1:6]
write.csv2(shannon.fut1, 'data/envf1_buffered.csv', row.names=F)

rm(shannon.fut1)

env.list <- readRDS("objects/envf5.rds")
columnnames <- c('tmean', 'tmax', 'tmin', 'tseas', 'pr', 'prseas', 'elev', 'slope', 'linear_infrastr', 'streams',
                 'forest', 'water', 'shrub', 'bare', 'for_shr_bare', 'settle_med', 'settle_low', 'for_shr_grass',
                 'for_shr_crop', 'for_shr_agric', 'mosaic_low', 'settle_high', 'crop_med', 'grass_low', 
                 'crop_high', 'grass_med', 'grass_high', 'mosaic_med', 'ext_perm_crop', 'crop_low', 'int_perm_crop', 
                 'mosaic_high')
env.stack <- stack(env.list)
                 
shannon.fut5 <- foreach(i = 1:nrow(centroids),
                        .packages=c("sf", "dplyr"),
                        .combine='rbind',
                        .verbose=F) %dopar% {
 
   if (small.buffer[i, 'shannon_index'] >= median.shannon) { # assign small buffer
     
     # generate buffer
     buffer <- st_buffer(centroids[i,], dist = 2000)
     # extract values in buffer
     values <- na.omit(as.data.frame(raster::extract(env.stack, buffer)))
     
     # calculate mean for continuous
     mean_cont <- colMeans(values[,1:6])
     
     return(mean_cont)
     
   }
   else { # calculate large buffer and assign
     
     # generate buffer
     buffer <- st_buffer(centroids[i,], dist = 4500)
     # extract values in buffer
     values <- na.omit(as.data.frame(raster::extract(env.stack, buffer)))
     
     # calculate mean for continuous 
     mean_cont <- colMeans(values[,1:6])
     
     return(mean_cont)
   }                        
}

colnames(shannon.fut5) <- columnnames[1:6]
write.csv2(shannon.fut5, 'data/envf5_buffered.csv', row.names=F)

# visualize variables in env2018 to check
raster_list <- list()

xy <- readRDS("objects/env2018.rds") %>% stack() %>% as.data.frame(xy=T) %>%
  na.omit() %>% dplyr::select(x,y) %>% rename(X=x, Y=y)

env2018 <- read.csv2('data/env2018_buffered.csv')

# iterate through all cols, generate rasterlayers and store
for (i in 1:ncol(env2018[, 1:29])) {
  column_name <- colnames(env2018)[i]
  raster_layer <- rasterFromXYZ(data.frame(xy, env2018[[i]]),
                                res = c(0.008333333, 0.008333333), crs = 4326)
  raster_list[[column_name]] <- raster_layer
}

# generate a plot for each raster
plot_list <- lapply(names(raster_list), function(layer_name) {
  
  raster_df <- as.data.frame(raster_list[[layer_name]], xy = TRUE, na.rm = TRUE)
  colnames(raster_df)[3] <- "value"
  
  ggplot(raster_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_viridis_c() +  
    theme_minimal() +
    labs(title = paste(layer_name), x = "Longitude", y = "Latitude")
})

# combine
env2018_vars <- grid.arrange(grobs = plot_list, ncol = 7)
ggsave('data/env2018_variables.jpg', env2018_vars, width = 35, height = 15)

# # now we need to calculate spatial predictors with the same method we used for modeling data
# 
# env.list <- readRDS("objects/env2018.rds")
# env.stack <- stack(env.list)
# env.df <- na.omit(as.data.frame(env.stack, xy=T))
# 
# # coordinates of the cases
# xy <- env.df[, c("x", "y")]
# # distance matrix
# distance.matrix <- as.matrix(dist(xy))
# # distance thresholds (same units as distance_matrix)
# # same units as projection!!!!! (degrees)
# distance.thresholds <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
# S
# pca <- pca_multithreshold(
#   distance.matrix = distance.matrix,
#   distance.thresholds = distance.thresholds,
#   max.spatial.predictors = 10)
