
# model with final hyperparameters and cv method 
# importance calculation and pdps

rm(list=ls())
source('scripts/utils.R')

# load dataset
dataset <- read.csv('data/modeling_data_v4.csv')
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

# load territory data and add missing variables to predict (x, y, total, logtotal, ano_cs)

# xy <- readRDS('objects/env2018.rds') %>% stack() %>% as.data.frame(xy=T) %>%
#   na.omit() %>% dplyr::select(x,y) %>% rename(X=x, Y=y)
# 
# missing <- data.frame(logTotal = rep(0, nrow(xy)),
#                       Total = rep(0, nrow(xy)),
#                       Ano_CS = rep(0, nrow(xy)))
# 
# env2018 <- read.csv2('data/env2018_buffered.csv') %>%
#   cbind(xy, missing) %>%
#   dplyr::select(-tmax, -tmin, -elev, -grass_high, -mosaic_high)
# 
# envf1 <- read.csv2('data/envf1_buffered.csv') %>%
#   dplyr::select(-tmax, -tmin) %>%
#   cbind(env2018[,5:33])
# 
# envf5 <- read.csv2('data/envf5_buffered.csv') %>%
#   dplyr::select(-tmax, -tmin) %>%
#   cbind(env2018[,5:33])

# create results dir
if (!dir.exists('results/random_cv_v2')) {
  dir.create('results/random_cv_v2')
  cat('results/random_cv_v2 created\n')
} else {
  cat('dir results/random_cv_v2 already exists\n')
}

# load best hyperparameters
hyperpar <- read.csv2('results/hyperparameters_mean.csv')
trees <- read.csv2('results/hyperparameters_ntrees_mean.csv')

mtry <- hyperpar$mtry[which.max(hyperpar$Rsquared_test)]
nodesize <- hyperpar$nodesize[which.max(hyperpar$Rsquared_test)]
ntree <- trees$ntree[which.max(trees$Rsquared_test)]

# can't parallelize because territory predictions takes almost all ram in each iteration

tic()
for (i in 1:10) {
          
    # set unique seed for each fold
    set.seed(i)
    
    # use caret to randomly split train/test data (80-20)
    index <- createDataPartition(dataset$logTotal, p = 0.8, list = FALSE)
    train <- dataset[index, ]
    test <- dataset[-index, ]
    
    # train the model
    model <- ranger(logTotal ~ . -Total -X -Y -Ano_CS -linear_infrastr,
                    data = train,
                    num.tree = ntree,
                    mtry = mtry,
                    min.node.size = nodesize,
                    importance = 'permutation')
    
    # store model
    filepath <- paste0('results/random_cv_v2/model_', i, '.rds')
    saveRDS(model, filepath)
    
    # generate training predictions 
    preds <- predict(model, data=train)
    train$pred <- preds[[1]]
    
    # save training preds
    filepath <- paste0('results/random_cv_v2/train_', i, '.csv')
    write.csv2(preds, filepath, row.names=T)
    
    # calculate metrics 
    mae_train <- MAE(preds[[1]], train$logTotal)
    mse_train <- mean((preds[[1]] - train$logTotal)^2)
    rmse_train <- RMSE(pred = preds[[1]], obs = train$logTotal)
    Rsq_train <- 1 - (sum((train$logTotal - preds[[1]])^2) / sum((train$logTotal - mean(train$logTotal))^2))
    
    # generate testing predictions
    preds <- predict(model, data=test)
    test$pred <- preds[[1]]
    
    # save testing preds
    filepath <- paste0('results/random_cv_v2/test_', i, '.csv')
    write.csv2(preds, filepath, row.names=T)
    
    # calculate testing metrics
    mae_test <- MAE(preds[[1]], test$logTotal)
    mse_test <- mean((preds[[1]] - test$logTotal)^2)
    rmse_test <- RMSE(pred = preds[[1]], obs = test$logTotal)
    Rsq_test <- 1 - (sum((test$logTotal - preds[[1]])^2) / sum((test$logTotal - mean(test$logTotal))^2))
    
    # store all these values in the metrics dataframe
    metrics <- c(i, mae_train, mse_train, rmse_train, Rsq_train,
                 mae_test, mse_test, rmse_test, Rsq_test)
    filepath <- paste0('results/random_cv_v2/metrics_', i, '.csv')
    write.csv2(metrics, filepath, row.names=F)
    
    # store variable importance
    varimp <- as.data.frame(t(model$variable.importance))
    filepath <- paste0('results/random_cv_v2/varimp_', i, '.csv')
    write.csv2(varimp, filepath, row.names=F)
    
    # # predict over whole territory
    # preds <- predict(model, data=env2018)[[1]]
    # filepath <- paste0('results/random_cv_v2/prespreds_', i, '.csv')
    # write.csv2(preds, filepath, row.names=F)
    # 
    # preds <- predict(model, data=envf1)[[1]]
    # filepath <- paste0('results/random_cv_v2/sustpreds_', i, '.csv')
    # write.csv2(preds, filepath, row.names=F)
    # 
    # preds <- predict(model, data=envf5)[[1]]
    # filepath <- paste0('results/random_cv_v2/fosspreds_', i, '.csv')
    # write.csv2(preds, filepath, row.names=F)
    # 
    # rm(preds)
    
    # generate PDPs for each predictor
    predictors <- names(train)[!(names(train) %in%
                                   c('Total', 'X', 'Y', 'Ano_CS', 'logTotal', 'pred', 'linear_infrastr'))]
    
    for (var in predictors) {
      
      pdp <- as.data.frame(partial(model, pred.var = var, train = train))
      filepath <- paste0('results/random_cv_v2/pdp_', var, '_fold_', i, '.csv')
      write.csv2(pdp, filepath, row.names = F)
      
    }
}
toc() # 1802.56 sec

# join all metrics rows into one dataframe
column_names <- c('fold','MAE_train','MSE_train','RMSE_train', 'Rsq_train',
                  'MAE_test','MSE_test','RMSE_test', 'Rsq_test')
metrics <- data.frame(matrix(ncol = length(column_names), nrow = 0))
file_list <- list.files('results/random_cv_v2', pattern = 'metrics_', full.names = T)
for (file in file_list) {
  metrics_row <- t(read.csv2(file))
  metrics <- rbind(metrics, metrics_row)
}
colnames(metrics) <- column_names
write.csv2(metrics, 'results/random_cv_v2/metrics.csv', row.names=F)

summary(metrics$Rsq_train)
summary(metrics$Rsq_test)
summary(metrics$RMSE_train)
summary(metrics$RMSE_test)

# join all variable importance
# load models filenames
file_list <- list.files('results/random_cv_v2', pattern = 'varimp', full.names = T)

# empty df to store accuracy values for each model
accu <- data.frame(matrix(nrow = 0, ncol = ncol(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y))))

# iterate through the model list 
for (i in 1:length(file_list)) {
  
  file <- read.csv2(file_list[[i]])
  accu <- rbind(accu, file)
}

# define variable categories for color legend
categories <- list(
  climatic = c('tmean', 'tseas', 'pr', 'prseas'),
  topographic = c('slope'),
  disturbance = c('linear_infrastr'),
  resources = c('streams', 'preys', 'Mean_BNDVI','StdDev_BNDVI'),
  diversity = c('shannon_index')
)

# vector with all variable names
all_vars <- colnames(accu)

# define land uses category by exclusion
land_uses <- setdiff(all_vars, unlist(categories))

# calculate absolute means of each variable
accu_mean_abs <- accu %>%
  summarise(across(everything(), ~ mean(abs(.)))) %>%
  pivot_longer(everything(), names_to = 'Variable', values_to = 'Mean_Abs') %>%
  arrange(Mean_Abs) # order ascending

# assign category to each variable
accu_mean_abs$category <- case_when(
  accu_mean_abs$Variable %in% categories$climatic ~ 'Climate',
  accu_mean_abs$Variable %in% categories$topographic ~ 'Topography',
  accu_mean_abs$Variable %in% categories$disturbance ~ 'Anthropic',
  accu_mean_abs$Variable %in% categories$resources ~ 'Resources',
  accu_mean_abs$Variable %in% categories$diversity ~ 'Landscape Diversity',
  accu_mean_abs$Variable %in% land_uses ~ 'Land Uses'
)

# vector of variables in the defined order
ordered_vars <- accu_mean_abs$Variable

# transform to long format and assign category
accu_long <- accu %>%
  pivot_longer(everything(), names_to = 'Variable', values_to = 'Value') %>%
  mutate(Variable = factor(Variable, levels = ordered_vars)) %>%
  mutate(category = case_when(
    Variable %in% categories$climatic ~ 'Climate',
    Variable %in% categories$topographic ~ 'Topography',
    Variable %in% categories$disturbance ~ 'Anthropic',
    Variable %in% categories$resources ~ 'Resources',
    Variable %in% categories$diversity ~ 'Landscape Diversity',
    TRUE ~ 'Land Uses'
  ))

# define colors for each category
colors <- c('Climate' = '#95b8f6', 
            'Topography' = '#f9d99a', 
            'Anthropic' = '#dcd9f8',
            'Land Uses' = '#fa5f49',
            'Resources' = 'mediumpurple1',
            'Landscape Diversity' = 'palegreen')

# define human-readable labels for variables
legend <- read.csv('spatial_data/clc_2018/legend.csv')

variable_labels <- c(
  'tmean' = 'Mean Temperature',
  'tseas' = 'Temperature Seasonality',
  'pr' = 'Precipitation',
  'prseas' = 'Precipitation Seasonality',
  'slope' = 'Slope',
  'linear_infrastr' = 'Linear Infrastructures',
  'streams' = 'Water Streams',
  'preys' = 'Preys Richness',
  'shannon_index' = 'Landscape Shannon Index',
  'forest' = 'Forest',
  'water' = 'Water Body',
  'shrub' = 'Shrub',
  'bare' = 'Bare Soil',
  'for_shr_bare' = 'Forest, Shrub & Bare',
  'for_shr_agric' = 'Forest, Shrub & Agric.',
  'settle_med' = 'Medium Int. Setllements',
  'settle_low' = 'Low Int. Settlements',
  'for_shr_grass' = 'Forest, Shrub & Grasslands',
  'for_shr_crop' = 'Forest, Shrub & Rain. Crops',
  'mosaic_low' = 'Low Int. Agriculture Mosaic',
  'settle_high' = 'High Int. Settelements',
  'crop_med' = 'Medium Int. Rainfed Crops',
  'grass_low' = 'Low Int. Grasslands',
  'crop_high' = 'High Int. Rainfed Crops',
  'grass_med' = 'Medium Int. Grasslands',
  'mosaic_med' = 'Medium Int. Agriculture Mosaic',
  'ext_perm_crop' = 'Extensive Permanent Crops',
  'crop_low' = 'Low Int. Rainfed Crops',
  'int_perm_crop' = 'Intensive Permanent Crops'
)

# Create and save boxplots
varimp_plot <- ggplot(accu_long, aes(x = Variable, y = abs(Value), fill = category)) +
  geom_boxplot() +
  coord_flip() +
  theme_minimal() +
  labs(title = 'Abs values of Mean Decrease in Accuracy',
       x = 'Variable', y = 'Mean Decrease in Accuracy', fill = 'Category') +
  scale_fill_manual(values = colors)+
  scale_x_discrete(labels = variable_labels)  # Replace x-axis labels with descriptive names

# Save the plot
ggsave('results/random_cv_v2/varimp.jpg', varimp_plot, width = 9, height = 6)

# # join all preds
# 
# present.df <- data.frame(matrix(ncol = 0, nrow = nrow(env2018)))
# file_list <- list.files('results/random_cv_v2', pattern = 'prespreds_', full.names = T)
# for (file in file_list) {
#   pred <- read.csv2(file)
#   present.df <- cbind(present.df, pred)
# }
# colnames(present.df) <- c('x1','x2','x3','x4','x5','x6','x7','x8','x9','x10')
# present.df$mean <- rowMeans(present.df)
# write.csv2(present.df, 'results/random_cv_v2/present_preds.csv', row.names=F)
# # present.df <- read.csv2('results/random_cv_v2/present_preds.csv')
# 
# present <- rasterFromXYZ(data.frame(xy, present.df$mean),
#                          res = c(0.008333333, 0.008333333), crs = 4326)
# writeRaster(present, 'results/random_cv_v2/present_prediction.asc', format='ascii', overwrite=T)
# 
# sustainable.df <- data.frame(matrix(ncol = 0, nrow = nrow(env2018)))
# file_list <- list.files('results/random_cv_v2', pattern = 'sustpreds_', full.names = T)
# for (file in file_list) {
#   pred <- read.csv2(file)
#   sustainable.df <- cbind(sustainable.df, pred)
# }
# colnames(sustainable.df) <- c('x1','x2','x3','x4','x5','x6','x7','x8','x9','x10')
# sustainable.df$mean <- rowMeans(sustainable.df)
# write.csv2(sustainable.df, 'results/random_cv_v2/sustainable_preds.csv', row.names=F)
# # sustainable.df <- read.csv2('results/random_cv_v2/sustainable_preds.csv')
# 
# sustainable <- rasterFromXYZ(data.frame(xy, sustainable.df$mean),
#                          res = c(0.008333333, 0.008333333), crs = 4326)
# writeRaster(sustainable, 'results/random_cv_v2/sustainable_prediction.asc', format='ascii', overwrite=T)
# 
# fossil.df <- data.frame(matrix(ncol = 0, nrow = nrow(env2018)))
# file_list <- list.files('results/random_cv_v2', pattern = 'fosspreds_', full.names = T)
# for (file in file_list) {
#   pred <- read.csv2(file)
#   fossil.df <- cbind(fossil.df, pred)
# }
# colnames(fossil.df) <- c('x1','x2','x3','x4','x5','x6','x7','x8','x9','x10')
# fossil.df$mean <- rowMeans(fossil.df)
# write.csv2(fossil.df, 'results/random_cv_v2/fossil_preds.csv', row.names=F)
# # fossil.df <- read.csv2('results/random_cv_v2/fossil_preds.csv')
# 
# fossil <- rasterFromXYZ(data.frame(xy, fossil.df$mean),
#                              res = c(0.008333333, 0.008333333), crs = 4326)
# writeRaster(fossil, 'results/random_cv_v2/fossil_prediction.asc', format='ascii', overwrite=T)
# 
# present.plot.df <- as.data.frame(present, xy=T)
# present.plot <- ggplot(present.plot.df, aes(x = x, y = y, fill = present.df.mean)) +
#   geom_raster() +
#   scale_fill_viridis_c(na.value = 'transparent', name = "Prediction") +
#   coord_equal() +  # Keep this to maintain aspect ratio
#   theme_bw() +  # This theme includes a light background grid
#   theme(
#     panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
#     panel.grid.minor = element_line(color = "gray90", linetype = "dotted"),
#     legend.position = c(0.95, 0.05),  # This positions the legend in the bottom right
#     legend.justification = c(1, 0),   # This aligns the legend to the bottom right
#     legend.box.background = element_rect(color = "black", size = 0.5),
#     legend.box.margin = margin(6, 6, 6, 6)
#   ) +
#   labs(title = "Present Predictions",
#        x = "Longitude",
#        y = "Latitude")
# 
# sustainable.plot.df <- as.data.frame(sustainable, xy=T)
# sustainable.plot <- ggplot(sustainable.plot.df, aes(x = x, y = y, fill = sustainable.df.mean)) +
#   geom_raster() +
#   scale_fill_viridis_c(na.value = 'transparent', name = "Prediction") +
#   coord_equal() +  # Keep this to maintain aspect ratio
#   theme_bw() +  # This theme includes a light background grid
#   theme(
#     panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
#     panel.grid.minor = element_line(color = "gray90", linetype = "dotted"),
#     legend.position = c(0.95, 0.05),  # This positions the legend in the bottom right
#     legend.justification = c(1, 0),   # This aligns the legend to the bottom right
#     legend.box.background = element_rect(color = "black", size = 0.5),
#     legend.box.margin = margin(6, 6, 6, 6)
#   ) +
#   labs(title = "Sustainable Predictions",
#        x = "Longitude",
#        y = "Latitude")
# 
# fossil.plot.df <- as.data.frame(fossil, xy=T)
# fossil.plot <- ggplot(fossil.plot.df, aes(x = x, y = y, fill = fossil.df.mean)) +
#   geom_raster() +
#   scale_fill_viridis_c(na.value = 'transparent', name = "Prediction") +
#   coord_equal() +  # Keep this to maintain aspect ratio
#   theme_bw() +  # This theme includes a light background grid
#   theme(
#     panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
#     panel.grid.minor = element_line(color = "gray90", linetype = "dotted"),
#     legend.position = c(0.95, 0.05),  # This positions the legend in the bottom right
#     legend.justification = c(1, 0),   # This aligns the legend to the bottom right
#     legend.box.background = element_rect(color = "black", size = 0.5),
#     legend.box.margin = margin(6, 6, 6, 6)
#   ) +
#   labs(title = "Fossil-Fueled Predictions",
#        x = "Longitude",
#        y = "Latitude")
# 
# predictions <- grid.arrange(present.plot, sustainable.plot, fossil.plot, ncol=3)
# ggsave('results/random_cv_v2/pres_sust_foss_preds.jpg', predictions, width=18, height=10)
# 
# # visualize all preds per model
# # iterate through all cols, generate rasterlayers and store
# raster_list <- list()
# 
# for (i in 1:ncol(present.df[,1:10])) {
#   column_name <- colnames(present.df)[i]
#   raster_layer <- rasterFromXYZ(data.frame(xy, present.df[[i]]),
#                                 res = c(0.008333333, 0.008333333), crs = 4326)
#   raster_list[[column_name]] <- raster_layer
# }
# 
# # generate a plot for each raster
# plot_list <- lapply(names(raster_list), function(layer_name) {
#   
#   raster_df <- as.data.frame(raster_list[[layer_name]], xy = TRUE, na.rm = TRUE)
#   colnames(raster_df)[3] <- 'value'
#   
#   ggplot(raster_df, aes(x = x, y = y, fill = value)) +
#     geom_raster() +
#     scale_fill_viridis_c() +  
#     theme_minimal() +
#     labs(title = paste(layer_name), x = 'Longitude', y = 'Latitude')
# })
# 
# # combine
# all_pres_preds <- grid.arrange(grobs = plot_list, ncol = 5)
# ggsave('results/random_cv_v2/all_present_predictions.jpg', all_pres_preds, width = 30, height = 10)

# # see shannon landscape diversity in modeling points
# dataset_sf$shannon <- ifelse(dataset$shannon_index>= median(dataset$shannon_index), 1, 0)
# mapview(dataset_sf, zcol='shannon')

# generate plots for all the pdps 

predictors <- names(dataset)[!(names(dataset) %in% c('Total', 'X', 'Y', 'Ano_CS', 'logTotal', 'linear_infrastr'))]

pdp_plot <- function(predictor) {
  predictions_df <- data.frame(matrix(ncol = 2, nrow = 0))
  
  file_list <- list.files('results/random_cv_v2', pattern = paste0('pdp_', predictor,'_'), full.names = T)
  
  # read and concatenate
  for (file in file_list) {
    prec_file <- read.csv2(file)
    predictions_df <- rbind(predictions_df, prec_file)
  }
  
  # calculate mean and stddev
  summary_stats <- predictions_df %>%
    group_by(!!sym(predictor)) %>%
    summarise(mean_yhat = mean(yhat), sd_yhat = sd(yhat))
  
  # plot
  p <- ggplot(summary_stats, aes_string(x = predictor, y = 'mean_yhat')) +
    geom_line() + 
    geom_ribbon(aes(ymin = mean_yhat - sd_yhat, ymax = mean_yhat + sd_yhat), alpha = 0.3) +
    theme_bw() +
    labs(x = predictor, y = 'log10 Total') +
    ggtitle(predictor) # name of the variable
  
  return(p) 
}

pdp_plots <- list()
# generate and save plots in list
for (predictor in predictors) {
  pdp_plots[[predictor]] <- pdp_plot(predictor)
}

# combine
all_pdps <- grid.arrange(grobs = pdp_plots, ncol = 8)
ggsave('results/random_cv_v2/pdps.png', plot = all_pdps, width = 20, height = 10)
