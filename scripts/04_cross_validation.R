
library(ranger)
library(caret)
library(dplyr)
library(tictoc)
library(doParallel)
library(plotly)
library(ggplot2)
library(MazamaSpatialUtils)
library(mapview)

# spatial block cross validation with iterations

rm(list=ls())

# load dataset
dataset <- read.csv2("data/modeling_data.csv")
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

y.max <- max(dataset$Y)
y.min <- min(dataset$Y)
folds <- seq(y.min, y.max, by = (y.max-y.min)/5)

dataset$fold <- cut(dataset$Y, 
                    breaks = folds, 
                    labels = FALSE, 
                    include.lowest = TRUE, 
                    right = FALSE)
# check
dataset_sf$fold <- dataset$fold
mapview(dataset_sf, zcol='fold')

# now we have 10 spatially 'independent' blocks to cross validate
# folds ordered by id so there's a latitudinal gradient

# new dir for latitudinal cross validation
dir.create('results/lat_cv')

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

ntree = 1500
nodesize = 10
mtry = 22

foreach(i = 1:length(unique(dataset$fold)),
       .packages=c("ranger", "sf", "dplyr", "caret"),
       .combine='rbind',
       verbose=T) %dopar% {
         
         # split train/test per fold
         indexes <- which(dataset$fold == i)
         test <- dataset[indexes,]
         train <- dataset[-indexes,]
         
         # set unique seed for each fold
         set.seed(i**2)
         
         # train the model
         model <- ranger(Total ~ . -X -Y -Ano_CS -fold,
                         data = train,
                         num.tree = ntree,
                         mtry = mtry,
                         min.node.size = nodesize)
         
         # generate training predictions (only select probs of presence)
         preds <- predict(model, data=train)
         train$pred <- preds[[1]]
         
         # calculate metrics (caret)
         mae_train <- MAE(preds[[1]], train$Total)
         mse_train <- mean((preds[[1]] - train$Total)^2)
         rmse_train <- RMSE(preds[[1]], train$Total)
         Rsq_train <- 1 - (sum((train$Total - preds[[1]])^2) / sum((train$Total - mean(train$Total))^2))
         
         # generate testing predictions and binaries from maxSSS
         preds <- predict(model, data=test)
         test$pred <- preds[[1]]
         
         # calculate metrics (caret)
         mae_test <- MAE(preds[[1]], test$Total)
         mse_test <- mean((preds[[1]] - test$Total)^2)
         rmse_test <- RMSE(preds[[1]], test$Total)
         Rsq_test <- 1 - (sum((test$Total - preds[[1]])^2) / sum((test$Total - mean(test$Total))^2))
         
         # store all the relevant values in the metrics dataframe
         metrics <- c(i, mae_train, mse_train, rmse_train, Rsq_train,
                      mae_test, mse_test, rmse_test, Rsq_test)
         filepath <- paste0('results/lat_cv/metrics_', i, '.csv')
         write.csv2(metrics, filepath, row.names=F)
       }

# join all metrics rows into one dataframe
column_names <- c('lat_S/N','MAE_train','MSE_train','RMSE_train', 'Rsq_train',
                  'MAE_test','MSE_test','RMSE_test','Rsq_test')
metrics <- data.frame(matrix(ncol = length(column_names), nrow = 0))
file_list <- list.files("results/lat_cv", pattern = 'metrics_', full.names = T)
for (file in file_list) {
  metrics_row <- t(read.csv2(file))
  metrics <- rbind(metrics, metrics_row)
}
colnames(metrics) <- column_names
write.csv2(metrics, "results/lat_cv/metrics.csv", row.names=F)

# latitudinal cross validation performs poorly in all folds except the southest 
# we're trying clustering points per environmental conditions

rm(list=ls())

# load dataset
dataset <- read.csv2("data/modeling_data.csv")
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

# variables to cluster
cols <- colnames(dataset[,6:49])
datatocluster <- dataset[, c('X', 'Y', cols)]

set.seed(123) 
# elbow method to choose cluster number
wss <- (nrow(datatocluster)-1)*sum(apply(datatocluster, 2, var))
for (i in 2:15) wss[i] <- sum(kmeans(datatocluster, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Número de clusters", ylab="Suma de los cuadrados dentro del cluster (WSS)")
# 5 seems fine

k <- 5 # folds
kmeans_result <- kmeans(datatocluster, centers = k, nstart = 25)
# assign clusters
dataset$cluster <- kmeans_result$cluster
table(dataset$cluster) # still imbalanced

# visualize
dataset_sf <- cbind(dataset_sf, cluster = dataset$cluster)
mapview(dataset_sf, zcol='cluster')
# clusters have some spatial pattern so might be coherent

# new dir for kmeans cross validation
dir.create('results/kmeans_cv')

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

ntree = 1500
nodesize = 10
mtry = 22

metrics <- foreach(i = 1:length(unique(dataset$cluster)),
                  .packages=c("ranger", "sf", "dplyr", "caret"),
                  .combine='rbind') %dopar% {
                    
                    # split train/test per cluster
                    indexes <- which(dataset$cluster == i)
                    test <- dataset[indexes,]
                    train <- dataset[-indexes,]
                    
                    # set unique seed for each cluster
                    set.seed(i**2)
                    
                    # train the model
                    model <- ranger(Total ~ . -X -Y -Ano_CS -cluster,
                                    data = train,
                                    num.tree = ntree,
                                    mtry = mtry,
                                    min.node.size = nodesize)
                    
                    # generate training predictions (only select probs of presence)
                    preds <- predict(model, data=train)
                    train$pred <- preds[[1]]
                    
                    # calculate metrics (caret)
                    mae_train <- MAE(preds[[1]], train$Total)
                    mse_train <- mean((preds[[1]] - train$Total)^2)
                    rmse_train <- RMSE(preds[[1]], train$Total)
                    Rsq_train <- 1 - (sum((train$Total - preds[[1]])^2) / sum((train$Total - mean(train$Total))^2))
                    
                    # generate testing predictions and binaries from maxSSS
                    preds <- predict(model, data=test)
                    test$pred <- preds[[1]]
                    
                    # calculate metrics (caret)
                    mae_test <- MAE(preds[[1]], test$Total)
                    mse_test <- mean((preds[[1]] - test$Total)^2)
                    rmse_test <- RMSE(preds[[1]], test$Total)
                    Rsq_test <- 1 - (sum((test$Total - preds[[1]])^2) / sum((test$Total - mean(test$Total))^2))
                    
                    # store all the relevant values in the metrics dataframe
                    metrics <- c(i, mae_train, mse_train, rmse_train, Rsq_train,
                                 mae_test, mse_test, rmse_test, Rsq_test)
                    return(metrics)
                  }

# join all metrics rows into one dataframe
metrics <- as.data.frame(metrics)
column_names <- c('cluster','MAE_train','MSE_train','RMSE_train', 'Rsq_train',
                  'MAE_test','MSE_test','RMSE_test','Rsq_test')
colnames(metrics) <- column_names
write.csv2(metrics, "results/metrics_kmeans.csv", row.names = F)

# predictions for all clusters look just fine