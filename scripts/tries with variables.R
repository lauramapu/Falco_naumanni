
columnnames <- c('tmean', 'tmax', 'tmin', 'tseas', 'pr', 'prseas', 'elev', 'slope', 'linear_infrastr', 'streams',
                 'forest', 'water', 'shrub', 'bare', 'for_shr_bare', 'settle_med', 'settle_low', 'for_shr_grass',
                 'for_shr_crop', 'for_shr_agric', 'mosaic_low', 'settle_high', 'crop_med', 'grass_low', 
                 'crop_high', 'grass_med', 'grass_high', 'mosaic_med', 'ext_perm_crop', 'crop_low', 'int_perm_crop', 
                 'mosaic_high')

####################################################################
#############################################################################
######################################################################################
# primero todas las variables con el buffer que corresponda y en local todas menos las climaticas y la elevacion
######################################################################################
#############################################################################
####################################################################

env_df_0 <- read.csv2("data/fnaumanni_0km.csv")
env_df_15 <- read.csv2("data/fnaumanni_15km.csv")
env_df_60 <- read.csv2("data/fnaumanni_60km.csv")

env_df <- read.csv2("data/fnaumanni_envdata.csv")

# we're choosing 15km2 buffer when H' > mean(range(H')) and 60km2 when lower
# (lower habitat area when diversity is high and viceversa)
medium_shannon <- mean(range(env_df$shannon_index))
mean_shannon <- mean(env_df$shannon_index)
median_shannon <- median(env_df$shannon_index)
# we're choosing median to cut obs in half
# mean and median are very similar

# add these variables to the final dataset
shannon_env <- matrix(ncol=32, nrow=0)

# select buffer size per observation by shannon index
for (i in 1:nrow(env_df)) {
  if (env_df[i, 'shannon_index'] >= median_shannon) {
    shannon_env <- rbind(shannon_env, as.numeric(st_drop_geometry(env_df_15[i,5:36])))
  }
  else if (env_df[i, 'shannon_index'] < median_shannon) {
    shannon_env <- rbind(shannon_env, as.numeric(st_drop_geometry(env_df_60[i,5:36])))
  }
}

#shannon_df <- as.data.frame(shannon_env)
colnames(shannon_env) <- columnnames

env_final <- cbind(env_df[,c(1:4, 101:102, 12:36)],
                   shannon_env) 
# calculate std dev of each column
std_devs <- sapply(env_final, sd, na.rm = TRUE)
# find the ones with zero std dev (no change)
cols_to_remove <- names(std_devs[std_devs == 0])
env_final <- env_final %>% dplyr::select(-all_of(cols_to_remove))

# CORRELATION TEST

# all cols as numeric
env_final[ , sapply(env_final, is.numeric)] <- lapply(env_final[ , sapply(env_final, is.numeric)], as.numeric)

# select only env cols
corrdata <- env_final[,7:57]

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
# we erase tmin, tmax and elev
env_final <- env_final %>% dplyr::select(-tmin, -tmax, -elev)

# models
dataset <- env_final %>%
  rename(x=X, y=Y)
dataset[ , sapply(dataset, is.numeric)] <- lapply(dataset[ , sapply(dataset, is.numeric)], as.numeric)

# names of the response variable and the predictors
dependent.variable.name <- "logTotal"
predictor.variable.names <- colnames(dplyr::select(dataset, -x, -y, -Total, -logTotal, -Ano_CS))

# coordinates of the cases
xy <- dataset[, c("x", "y")]

# distance matrix
distance.matrix <- as.matrix(dist(xy))

# distance thresholds (same units as distance_matrix)
# same units as projection!!!!! (degrees)
distance.thresholds <- c(0.05, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4, 0.5, 1)

# random seed for reproducibility
random.seed <- 21

# moran's index per distance thresholds
spatialRF::plot_training_df_moran(
  data = dataset,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  fill.color = viridis::viridis(
    100,
    option = "F",
    direction = -1
  ),
  point.color = "gray40"
)

# FITTING A SPATIAL MODEL

# first non-spatial to check moran index 
model.non.spatial <- spatialRF::rf(
  data = dataset,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = FALSE
)

spatialRF::plot_moran(
  model.non.spatial, 
  verbose = TRUE
)
# max dist and min p-value at ~15km (0.10 degrees)

model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = T,
  seed = random.seed
)
# The model residuals are not spatially correlated, there is no need to fit a spatial model

spatialRF::plot_moran(
  model.spatial, 
  verbose = FALSE
)

plot_response_curves(
  model = model.spatial,
  variables = NULL,
  quantiles = c(0.1, 0.5, 0.9),
  grid.resolution = 200,
  line.color = viridis::viridis(length(quantiles), option = "F", end = 0.9),
  ncol = 5,
  show.data = FALSE,
  verbose = TRUE
)

####################################################################
#############################################################################
######################################################################################
# ahora solo con las variables de buffer que correspondan
######################################################################################
#############################################################################
####################################################################

rm(list=ls())

env_df_0 <- read.csv2("data/fnaumanni_0km.csv")
env_df_15 <- read.csv2("data/fnaumanni_15km.csv")
env_df_60 <- read.csv2("data/fnaumanni_60km.csv")

env_df <- read.csv2("data/fnaumanni_envdata.csv")

# we're choosing 15km2 buffer when H' > mean(range(H')) and 60km2 when lower
# (lower habitat area when diversity is high and viceversa)
medium_shannon <- mean(range(env_df$shannon_index))
mean_shannon <- mean(env_df$shannon_index)
median_shannon <- median(env_df$shannon_index)
# we're choosing median to cut obs in half
# mean and median are very similar

# add these variables to the final dataset
shannon_env <- matrix(ncol=32, nrow=0)

# select buffer size per observation by shannon index
for (i in 1:nrow(env_df)) {
  if (env_df[i, 'shannon_index'] >= median_shannon) {
    shannon_env <- rbind(shannon_env, as.numeric(st_drop_geometry(env_df_15[i,5:36])))
  }
  else if (env_df[i, 'shannon_index'] < median_shannon) {
    shannon_env <- rbind(shannon_env, as.numeric(st_drop_geometry(env_df_60[i,5:36])))
  }
}

#shannon_df <- as.data.frame(shannon_env)
colnames(shannon_env) <- columnnames

env_final <- cbind(env_df[,c(1:4, 101:102)],
                   shannon_env) 
# calculate std dev of each column
std_devs <- sapply(env_final, sd, na.rm = TRUE)
# find the ones with zero std dev (no change)
cols_to_remove <- names(std_devs[std_devs == 0])
env_final <- env_final %>% dplyr::select(-all_of(cols_to_remove))

env_final$logTotal <- log10(env_final$Total)

# CORRELATION TEST

# all cols as numeric
env_final[ , sapply(env_final, is.numeric)] <- lapply(env_final[ , sapply(env_final, is.numeric)], as.numeric)

# select only env cols
corrdata <- env_final[,7:36] 

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
# we erase tmin, tmax and elev
env_final <- env_final %>% dplyr::select(-tmin, -tmax, -elev)

# models
dataset <- env_final %>%
  rename(x=X, y=Y)

# names of the response variable and the predictors
dependent.variable.name <- "logTotal"
predictor.variable.names <- colnames(dplyr::select(dataset, -x, -y, -Total, -logTotal, -Ano_CS))

# coordinates of the cases
xy <- dataset[, c("x", "y")]

# distance matrix
distance.matrix <- as.matrix(dist(xy))

# distance thresholds (same units as distance_matrix)
# same units as projection!!!!! (degrees)
distance.thresholds <- c(0.05, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4, 0.5, 1)

# random seed for reproducibility
random.seed <- 21

# FITTING A SPATIAL MODEL

# first non-spatial to check moran index 
model.non.spatial <- spatialRF::rf(
  data = dataset,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy, #not needed by rf, but other functions read it from the model
  seed = random.seed,
  verbose = FALSE
)

spatialRF::plot_moran(
  model.non.spatial, 
  verbose = TRUE
)
# max dist and min p-value at ~15km (0.10 degrees)

model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = T,
  seed = random.seed
)
# The model residuals are not spatially correlated, there is no need to fit a spatial model

# spatialRF::plot_moran(
#   model.spatial, 
#   verbose = FALSE
# )

print_performance(model.non.spatial)
print_importance(model.non.spatial)

plot_response_curves(
  model = model.non.spatial,
  variables = NULL,
  quantiles = c(0.1, 0.5, 0.9),
  grid.resolution = 200,
  line.color = viridis::viridis(length(c(0.1, 0.5, 0.9)), option = "F", end = 0.9),
  ncol = 5,
  show.data = FALSE,
  verbose = TRUE
)

# another try with random cross validation

model.repeated <- rf_repeat(
  model = model.non.spatial,
  data = dataset,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy,
  repetitions = 5,
  keep.models = TRUE,
  seed = 21,
  verbose = TRUE,
  n.cores = parallel::detectCores() - 1,
  cluster = NULL
)
print_performance(model.repeated) # very similar

# save this env_final as modeling data
write.csv2(env_final, "data/modeling_data_v2.csv", row.names=F)
