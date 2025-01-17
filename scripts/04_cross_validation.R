
# spatial block cross validation with iterations

rm(list=ls())
source('scripts/utils.R')

# load dataset
dataset <- read.csv2("data/modeling_data_v3.csv")
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

y.max <- max(dataset$Y)
y.min <- min(dataset$Y)
folds <- seq(y.min, y.max, by = (y.max-y.min)/10)

dataset$fold <- cut(dataset$Y, 
                    breaks = folds, 
                    labels = FALSE, 
                    include.lowest = TRUE, 
                    right = FALSE)
# check
dataset_sf$fold <- dataset$fold
mapview(dataset_sf, zcol='fold')

# now we have 5 spatially 'independent' blocks to cross validate
# folds ordered by id so there's a latitudinal gradient

# new dir for latitudinal cross validation
dir.create('results/lat_cv')

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

ntree = 1000
nodesize = 50
mtry = 5

metrics <- foreach(i = 1:length(unique(dataset$fold)),
                   .packages=c("ranger", "sf", "dplyr", "caret"),
                   .combine='rbind',
                   .verbose=T) %dopar% {
                     
                     # split train/test per fold
                     indexes <- which(dataset$fold == i)
                     test <- dataset[indexes,]
                     train <- dataset[-indexes,]
                     
                     # set unique seed for each fold
                     set.seed(i**2)
                     
                     # train the model
                     model <- ranger(logTotal ~ . -Total -X -Y -Ano_CS -fold,
                                     data = train,
                                     num.tree = ntree,
                                     mtry = mtry,
                                     min.node.size = nodesize,
                                     importance='none')
                     
                     # generate training predictions 
                     preds <- predict(model, data=train)
                     train$pred <- preds[[1]]
                     
                     # calculate metrics 
                     mae_train <- MAE(preds[[1]], train$logTotal)
                     mse_train <- mean((preds[[1]] - train$logTotal)^2)
                     rmse_train <- RMSE(pred = preds[[1]], obs = train$logTotal)
                     Rsq_train <- 1 - (sum((train$logTotal - preds[[1]])^2) / sum((train$logTotal - mean(train$logTotal))^2))
                     
                     # generate testing predictions
                     preds <- predict(model, data=test)
                     test$pred <- preds[[1]]
                     
                     # calculate testing metrics 
                     mae_test <- MAE(preds[[1]], test$logTotal)
                     mse_test <- mean((preds[[1]] - test$logTotal)^2)
                     rmse_test <- RMSE(pred = preds[[1]], obs = test$logTotal)
                     Rsq_test <- 1 - (sum((test$logTotal - preds[[1]])^2) / sum((test$logTotal - mean(test$logTotal))^2))
                     
                     # store all these values in the metrics dataframe
                     metrics <- c(i, mae_train, mse_train, rmse_train, Rsq_train,
                                  mae_test, mse_test, rmse_test, Rsq_test)
                     return(metrics)
                   }

# join all metrics rows into one dataframe
column_names <- c('lat_S/N','MAE_train','MSE_train','RMSE_train', 'Rsq_train',
                  'MAE_test','MSE_test','RMSE_test','Rsq_test')
colnames(metrics) <- column_names
mean(metrics[,'Rsq_test'])
write.csv2(metrics, "results/lat_cv/metrics.csv", row.names = F)

# latitude cross validation performs poorly 

# let's try latitude cv with balanced folds

rm(list=ls())
source('scripts/utils.R')

# load dataset
dataset <- read.csv2("data/modeling_data_v3.csv")
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

dataset$fold <- cut(dataset$Y,
                    breaks = quantile(dataset$Y, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
                    include.lowest = TRUE,
                    labels = 1:10)

# check
dataset_sf$fold <- dataset$fold
mapview(dataset_sf, zcol='fold')

# new dir for latitudinal cross validation
dir.create('results/lat_cv')

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

ntree = 1000
nodesize = 50
mtry = 5

metrics <- foreach(i = 1:length(unique(dataset$fold)),
                   .packages=c("ranger", "sf", "dplyr", "caret"),
                   .combine='rbind',
                   .verbose=T) %dopar% {
                     
                     # split train/test per fold
                     indexes <- which(dataset$fold == i)
                     test <- dataset[indexes,]
                     train <- dataset[-indexes,]
                     
                     # set unique seed for each fold
                     set.seed(i**2)
                     
                     # train the model
                     model <- ranger(logTotal ~ . -Total -X -Y -Ano_CS -fold,
                                     data = train,
                                     num.tree = ntree,
                                     mtry = mtry,
                                     min.node.size = nodesize,
                                     importance='none')
                     
                     # generate training predictions 
                     preds <- predict(model, data=train)
                     train$pred <- preds[[1]]
                     
                     # calculate metrics 
                     mae_train <- MAE(preds[[1]], train$logTotal)
                     mse_train <- mean((preds[[1]] - train$logTotal)^2)
                     rmse_train <- RMSE(pred = preds[[1]], obs = train$logTotal)
                     Rsq_train <- 1 - (sum((train$logTotal - preds[[1]])^2) / sum((train$logTotal - mean(train$logTotal))^2))
                     
                     # generate testing predictions
                     preds <- predict(model, data=test)
                     test$pred <- preds[[1]]
                     
                     # calculate testing metrics 
                     mae_test <- MAE(preds[[1]], test$logTotal)
                     mse_test <- mean((preds[[1]] - test$logTotal)^2)
                     rmse_test <- RMSE(pred = preds[[1]], obs = test$logTotal)
                     Rsq_test <- 1 - (sum((test$logTotal - preds[[1]])^2) / sum((test$logTotal - mean(test$logTotal))^2))
                     
                     # store all these values in the metrics dataframe
                     metrics <- c(i, mae_train, mse_train, rmse_train, Rsq_train,
                                  mae_test, mse_test, rmse_test, Rsq_test)
                     return(metrics)
                   }

# join all metrics rows into one dataframe
column_names <- c('lat_S/N','MAE_train','MSE_train','RMSE_train', 'Rsq_train',
                  'MAE_test','MSE_test','RMSE_test','Rsq_test')
colnames(metrics) <- column_names
mean(metrics[,'Rsq_test'])
write.csv2(metrics, "results/lat_cv/metrics_q.csv", row.names = F)

# let's try longitude cross validation

rm(list=ls())
source('scripts/utils.R')

# load dataset
dataset <- read.csv2("data/modeling_data_v3.csv")
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

dataset$fold <- cut(dataset$X,
                    breaks = quantile(dataset$X, probs = seq(0, 1, by = 0.1), na.rm = TRUE),
                    include.lowest = TRUE,
                    labels = 1:10)
# check
dataset_sf$fold <- dataset$fold
mapview(dataset_sf, zcol='fold')

# new dir for latitudinal cross validation
dir.create('results/lon_cv')

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

ntree = 1000
nodesize = 50
mtry = 5

metrics <- foreach(i = 1:length(unique(dataset$fold)),
                   .packages=c("ranger", "sf", "dplyr", "caret"),
                   .combine='rbind',
                   .verbose=T) %dopar% {
                     
                     # split train/test per fold
                     indexes <- which(dataset$fold == i)
                     test <- dataset[indexes,]
                     train <- dataset[-indexes,]
                     
                     # set unique seed for each fold
                     set.seed(i**2)
                     
                     # train the model
                     model <- ranger(logTotal ~ . -Total -X -Y -Ano_CS -fold,
                                     data = train,
                                     num.tree = ntree,
                                     mtry = mtry,
                                     min.node.size = nodesize,
                                     importance='none')
                     
                     # generate training predictions 
                     preds <- predict(model, data=train)
                     train$pred <- preds[[1]]
                     
                     # calculate metrics 
                     mae_train <- MAE(preds[[1]], train$logTotal)
                     mse_train <- mean((preds[[1]] - train$logTotal)^2)
                     rmse_train <- RMSE(pred = preds[[1]], obs = train$logTotal)
                     Rsq_train <- 1 - (sum((train$logTotal - preds[[1]])^2) / sum((train$logTotal - mean(train$logTotal))^2))
                     
                     # generate testing predictions
                     preds <- predict(model, data=test)
                     test$pred <- preds[[1]]
                     
                     # calculate testing metrics 
                     mae_test <- MAE(preds[[1]], test$logTotal)
                     mse_test <- mean((preds[[1]] - test$logTotal)^2)
                     rmse_test <- RMSE(pred = preds[[1]], obs = test$logTotal)
                     Rsq_test <- 1 - (sum((test$logTotal - preds[[1]])^2) / sum((test$logTotal - mean(test$logTotal))^2))
                     
                     # store all these values in the metrics dataframe
                     metrics <- c(i, mae_train, mse_train, rmse_train, Rsq_train,
                                  mae_test, mse_test, rmse_test, Rsq_test)
                     return(metrics)
                   }

# join all metrics rows into one dataframe
column_names <- c('lon_W/E','MAE_train','MSE_train','RMSE_train', 'Rsq_train',
                  'MAE_test','MSE_test','RMSE_test','Rsq_test')
colnames(metrics) <- column_names
mean(metrics[,'Rsq_test'])
write.csv2(metrics, "results/lon_cv/metrics_q.csv", row.names = F)

# better performance than latitude but still worse than random cv

# we're trying clustering points per environmental conditions

rm(list=ls())
source('scripts/utils.R')

# load dataset
dataset <- read.csv2("data/modeling_data.csv")
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

# variables to cluster
cols <- colnames(dplyr::select(dataset, -logTotal, -Total, -Ano_CS))
datatocluster <- dataset[, cols]

set.seed(123) 
# elbow method to choose cluster number
wss <- (nrow(datatocluster)-1)*sum(apply(datatocluster, 2, var))
for (i in 2:15) wss[i] <- sum(kmeans(datatocluster, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of clusters", ylab="Sum of squares in the cluster (WSS)")

k <- 7 # folds
kmeans_result <- kmeans(datatocluster, centers = k, nstart = 25)
# assign clusters
dataset$cluster <- kmeans_result$cluster
table(dataset$cluster) # pretty balanced

# visualize
dataset_sf <- cbind(dataset_sf, cluster = dataset$cluster)
mapview(dataset_sf, zcol='cluster')

# create results dir
if (!dir.exists('results/kmeans_cv')) {
  dir.create('results/kmeans_cv')
  cat('results/kmeans_cv created\n')
} else {
  cat('dir results/kmeans_cv already exists\n')
}

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

# load best hyperparameters
hyperpar <- read.csv2('results/hyperparameters_mean.csv')
trees <- read.csv2('results/hyperparameters_ntrees_mean.csv')

mtry <- hyperpar$mtry[which.max(hyperpar$Rsquared_test)]
nodesize <- hyperpar$nodesize[which.max(hyperpar$Rsquared_test)]
ntree <- trees$ntree[which.max(trees$Rsquared_test)]

metrics <- foreach(i = 1:length(unique(dataset$cluster)),
                  .packages=c("ranger", "sf", "dplyr", "caret"),
                  .combine='rbind',
                  .verbose=T) %dopar% {
                    
                    # split train/test per cluster
                    indexes <- which(dataset$cluster == i)
                    test <- dataset[indexes,]
                    train <- dataset[-indexes,]
                    
                    # set unique seed for each cluster
                    set.seed(i)
                    
                    # train the model
                    model <- ranger(logTotal ~ . -Total -X -Y -Ano_CS -cluster,
                                    data = train,
                                    num.tree = ntree,
                                    mtry = mtry,
                                    min.node.size = nodesize)
                    
                    # generate training predictions (only select probs of presence)
                    preds <- predict(model, data=train)
                    train$pred <- preds[[1]]
                    
                    # calculate metrics (caret)
                    mae_train <- MAE(preds[[1]], train$logTotal)
                    mse_train <- mean((preds[[1]] - train$logTotal)^2)
                    rmse_train <- RMSE(preds[[1]], train$logTotal)
                    Rsq_train <- 1 - (sum((train$logTotal - preds[[1]])^2) / sum((train$logTotal - mean(train$logTotal))^2))
                    
                    # generate testing predictions and binaries from maxSSS
                    preds <- predict(model, data=test)
                    test$pred <- preds[[1]]
                    
                    # calculate metrics (caret)
                    mae_test <- MAE(preds[[1]], test$logTotal)
                    mse_test <- mean((preds[[1]] - test$logTotal)^2)
                    rmse_test <- RMSE(preds[[1]], test$logTotal)
                    Rsq_test <- 1 - (sum((test$logTotal - preds[[1]])^2) / sum((test$logTotal - mean(test$logTotal))^2))
                    
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
mean(metrics[,'Rsq_test'])
write.csv2(metrics, "results/kmeans_cv/metrics.csv", row.names = F)

# predictions per env cluster show poor performance as well, so we'll try deleting as few obs per fold as possible with SLOO

rm(list=ls())
source('scripts/utils.R')

# load dataset
dataset <- read.csv2("data/modeling_data.csv")
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

# create results dir
if (!dir.exists('results/SLOO')) {
  dir.create('results/SLOO')
  cat('results/SLOO created\n')
} else {
  cat('dir results/SLOO already exists\n')
}

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

# load best hyperparameters
hyperpar <- read.csv2('results/hyperparameters_mean.csv')
trees <- read.csv2('results/hyperparameters_ntrees_mean.csv')

mtry <- hyperpar$mtry[which.max(hyperpar$Rsquared_test)]
nodesize <- hyperpar$nodesize[which.max(hyperpar$Rsquared_test)]
ntree <- trees$ntree[which.max(trees$Rsquared_test)]

tic('SLOO')
foreach(i = 1:nrow(dataset),
        .packages=c("ranger", "sf", "dplyr", "caret"),
        .combine='rbind',
        .verbose=T) %dopar% {
          
          # generate 20km buffer around test point
          buffer <- st_buffer(dataset_sf[i,], dist = 20000) %>% dplyr::select(geometry)
          
          # store all points inside the buffer
          to_erase <- st_intersection(dataset_sf, buffer)
          coordinates <- data.frame(X = st_coordinates(to_erase)[,'X'],
                                    Y = st_coordinates(to_erase)[,'Y'])
          to_erase <- to_erase %>% as.data.frame() %>% dplyr::select(-geometry)
          
          # select all points except to_erase to train
          train <- anti_join(dataset, to_erase, by = colnames(dataset))
          test <- to_erase
          
          # set unique seed for each fold
          set.seed(i)
          
          # train the model
          model <- ranger(logTotal ~ . -Total -X -Y -Ano_CS,
                          data = train,
                          num.tree = ntree,
                          mtry = mtry,
                          min.node.size = nodesize,
                          importance='none')
          
          # generate training predictions 
          preds <- predict(model, data=train)
          train$pred <- preds[[1]]
          
          # calculate metrics 
          mae_train <- MAE(preds[[1]], train$logTotal)
          mse_train <- mean((preds[[1]] - train$logTotal)^2)
          rmse_train <- RMSE(pred = preds[[1]], obs = train$logTotal)
          Rsq_train <- 1 - (sum((train$logTotal - preds[[1]])^2) / sum((train$logTotal - mean(train$logTotal))^2))
          
          # generate testing predictions
          preds <- predict(model, data=test)
          test$pred <- preds[[1]]
          
          # calculate testing metrics (no Rsq)
          mae_test <- MAE(preds[[1]], test$logTotal)
          mse_test <- mean((preds[[1]] - test$logTotal)^2)
          rmse_test <- RMSE(pred = preds[[1]], obs = test$logTotal)
          Rsq_test <- 1 - (sum((test$logTotal - preds[[1]])^2) / sum((test$logTotal - mean(test$logTotal))^2))

          # store all these values in the metrics dataframe
          metrics <- c(i, mae_train, mse_train, rmse_train, Rsq_train,
                       mae_test, mse_test, rmse_test, Rsq_test)
          filepath <- paste0('results/SLOO/metrics_', i, '.csv')
          write.csv2(metrics, filepath, row.names=F)
        }
toc() # 1226.16 sec

# join all metrics rows into one dataframe
column_names <- c('obs','MAE_train','MSE_train','RMSE_train', 'Rsq_train',
                  'MAE_test','MSE_test','RMSE_test','Rsq_test')
metrics <- data.frame(matrix(ncol = length(column_names), nrow = 0))
file_list <- list.files("results/SLOO", pattern = 'metrics_', full.names = T)
for (file in file_list) {
  metrics_row <- t(read.csv2(file))
  metrics <- rbind(metrics, metrics_row)
}
colnames(metrics) <- column_names
write.csv2(metrics, "results/SLOO/metrics.csv", row.names=F)

summary(metrics$Rsq_train)
summary(metrics$Rsq_test[is.finite(metrics$Rsq_test)])
summary(metrics$RMSE_train)
summary(metrics$RMSE_test)

# very bad results with SLOO

# we continue with random CV