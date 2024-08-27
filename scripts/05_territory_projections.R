# we need to generate predictions to all the territory
# to do this, we need env layers for all peninsular spain
# we're doing it for all raster centroids the same way we did for observations

library(sf)
library(raster)
library(mapview)
library(dplyr)

env2018 <- readRDS("objects/env2018.rds")
envf1 <- readRDS("objects/envf1.rds")
envf5 <- readRDS("objects/envf5.rds")

columnnames <- c('tmean', 'tmax', 'tmin', 'tseas', 'pr', 'prseas', 'elev', 'slope', 'linear_infrastr',
                 'forest', 'water', 'shrub', 'bare', 'for_shr_bare', 'settle_med', 'settle_low', 'for_shr_grass',
                 'for_shr_crop', 'for_shr_agric', 'mosaic_low', 'settle_high', 'crop_med', 'grass_low', 
                 'crop_high', 'grass_med', 'grass_high', 'mosaic_med', 'ext_perm_crop', 'crop_low', 'int_perm_crop', 
                 'mosaic_high')

# outputs:
# - env df 2018: no buffer, 2km buffer, 4.5km buffer
# - env df future 1: same
# - env df future 5: same
# 9 extractions in total

# convert to stack
stack2018 <- stack(env2018)
stackf1 <- stack(envf1)
stackf5 <- stack(envf5)

# local variables can be extracted to centroids just by applying asdataframe to the stacks
pres_0 <- na.omit(as.data.frame(stack2018, xy=T))
colnames(pres_0)[3:33] <- columnnames

fut1_0 <- na.omit(as.data.frame(stackf1, xy=T))
colnames(fut1_0)[3:33] <- columnnames

fut5_0 <- na.omit(as.data.frame(stackf5, xy=T))
colnames(fut5_0)[3:33] <- columnnames

# extract centroids (sf) from previous df
centroids <- pres_0 %>%
  st_as_sf(coords=c('x','y'), crs=crs(stack2018[[1]])) %>%
  dplyr::select(geometry)

# empty df
env1 <- data.frame(matrix(ncol=length(env2018), nrow=0))
env2 <- data.frame(matrix(ncol=6, nrow=0)) # only climatic 
env3 <- data.frame(matrix(ncol=6, nrow=0))
env4 <- data.frame(matrix(ncol=length(env2018), nrow=0))
env5 <- data.frame(matrix(ncol=6, nrow=0)) # only climatic 
env6 <- data.frame(matrix(ncol=6, nrow=0))

# here we also calculate shannon diversity index
for (i in 1:nrow(centroids)) {
  
  # generate buffer
  buffer <- st_buffer(centroids[i,], dist = 2000)
  # extract values in buffer
  values1 <- na.omit(as.data.frame(extract(stack2018, buffer)))
  values2 <- na.omit(as.data.frame(extract(stackf1[[1:6]], buffer)))
  values3 <- na.omit(as.data.frame(extract(stackf5[[1:6]], buffer)))
  
  # calculate mean for continuous and proportion for categorical
  mean_cont1 <- colMeans(values1[,1:9]) # continuous variables
  mean_cont2 <- colMeans(values2) 
  mean_cont3 <- colMeans(values3) 
  
  landuses <- colSums(values1[,10:31]) # presence of each land sue
  total <- sum(landuses)
  mean_cat <- landuses / total # proportion of ones vs total
  # calculate shannon diversity index (H') for land uses
  # max H' is log(22) = 3.091042
  shannon <- -sum(mean_cat[mean_cat > 0] * log(mean_cat[mean_cat > 0]))
  # bind
  mean_values <- c(mean_cont1, mean_cat, shannon)
  env1 <- rbind(env1, mean_values)
  env2 <- rbind(env2, mean_cont2)
  env3 <- rbind(env3, mean_cont3)
  
  # generate buffer
  buffer <- st_buffer(centroids[i,], dist = 4500)
  # extract values in buffer
  values4 <- na.omit(as.data.frame(extract(stack2018, buffer)))
  values5 <- na.omit(as.data.frame(extract(stackf1[[1:6]], buffer)))
  values6 <- na.omit(as.data.frame(extract(stackf5[[1:6]], buffer)))
  
  # calculate mean for continuous and proportion for categorical
  mean_cont4 <- colMeans(values4[,1:9]) # continuous variables
  mean_cont5 <- colMeans(values5) 
  mean_cont6 <- colMeans(values6) 
  
  landuses <- colSums(values4[,10:31]) # presence of each land sue
  total <- sum(landuses)
  mean_cat <- landuses / total # proportion of ones vs total

  # bind
  mean_values <- c(mean_cont4, mean_cat)
  env4 <- rbind(env4, mean_values)
  env5 <- rbind(env5, mean_cont5)
  env6 <- rbind(env6, mean_cont6)
}

pres <- cbind(pres_0, env1, env4)
colnames(pres)[3:33] <- columnnames
colnames(env_df_15)[36] <- 'shannon_index'
# scale buffer results
env_df_15[,14:35] <- scale(st_drop_geometry(env_df_15)[,14:35]) # zero variance variables get NA
env_df_15[is.na(env_df_15)] <- 0
# save
write.csv2(st_drop_geometry(env_df_15), "data/fnaumanni_15km.csv", row.names=F)

# 4.5km buffer (~64km2)
# empty df
env <- data.frame(matrix(ncol=length(env2016), nrow=0))

for (i in 1:nrow(centroids)) {
  buffer <- st_buffer(centroids[i,], dist = 4500)
  if (centroids$Ano_CS[i]==2016) {
    values <- na.omit(as.data.frame(extract(stack2016, buffer)))
  }
  else if (centroids$Ano_CS[i]==2017) {
    values <- na.omit(as.data.frame(extract(stack2017, buffer)))
  }
  else if (centroids$Ano_CS[i]==2018) {
    values <- na.omit(as.data.frame(extract(stack2018, buffer)))
  }
  mean_cont <- colMeans(values[,1:9]) 
  mean_cat <- colSums(values[,10:31]) / nrow(values) 
  mean_values <- c(mean_cont, mean_cat)
  env <- rbind(env, mean_values)
}

env_df_60 <- cbind(centroids, env)
colnames(env_df_60)[5:35] <- columnnames
# scale buffer results
env_df_60[,14:35] <- scale(st_drop_geometry(env_df_60)[,14:35])
env_df_60[is.na(env_df_60)] <- 0
# save
write.csv2(st_drop_geometry(env_df_60), "data/fnaumanni_60km.csv", row.names=F)

# add suffix to col names in each dataset and join
colnames(env_df_0)[5:35] <- paste(colnames(env_df_0)[5:35], "0", sep = "_")
colnames(env_df_15)[5:35] <- paste(colnames(env_df_15)[5:35], "15", sep = "_")
colnames(env_df_60)[5:35] <- paste(colnames(env_df_60)[5:35], "60", sep = "_")
env_df <- cbind(env_df_0, env_df_15[,5:35], env_df_60[,5:35], env_df_15[,'shannon_index'])
env_df <- env_df %>% dplyr::select(-geometry.1, -geometry.2, -geometry.3) %>% st_drop_geometry() # erase geometry

# add log column of response because hightly imbalanced
env_df$logTotal <- log(env_df$Total)

# save complete database
write.csv2(env_df, "data/fnaumanni_envdata.csv", row.names=F)

# modelling observations versus shannon index
# just to preeliminairly check if there's any relation
plot(env_df$logTotal, env_df$shannon_index)
model <- lm(logTotal ~ shannon_index, data = env_df)
summary(model)
# shannon index has no significant relation with abundance (p-value>0.1)

# we're choosing 15km2 buffer when H' > mean(range(H')) and 60km2 when lower
# (lower habitat area when diversity is high and viceversa)
medium_shannon <- mean(range(env_df$shannon_index))
mean_shannon <- mean(env_df$shannon_index)
median_shannon <- median(env_df$shannon_index)
# we're choosing median to cut obs in half

shannon_env <- matrix(ncol=31,nrow=0)
for (i in 1:nrow(env_df)) {
  if (env_df[i, 'shannon_index'] >= median_shannon) {
    shannon_env <- rbind(shannon_env, as.numeric(st_drop_geometry(env_df_15[i,5:35])))
  }
  else if (env_df[i, 'shannon_index'] < median_shannon) {
    shannon_env <- rbind(shannon_env, as.numeric(st_drop_geometry(env_df_60[i,5:35])))
  }
}
shannon_env_full <- cbind(env_df[,1:4], logTotal = env_df[,'logTotal'],
                          st_drop_geometry(env_df_0[,5:35]),
                          st_drop_geometry(env_df_15[,'shannon_index']),
                          shannon_env)
colnames(shannon_env_full)[38:68] <- columnnames
write.csv2(shannon_env_full, 'data/fnaumanni_byshannon.csv', row.names=F)
