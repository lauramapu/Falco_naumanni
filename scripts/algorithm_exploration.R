
rm(list=ls())
source('scripts/utils.R')

# Total, default hyperpar, non-spatial

dataset <- read.csv("data/modeling_data_v3.csv")

response_var <- "Total" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry = NULL,
                    nodesize = NULL,
                    ntree=1000)

print_results(results)

# calculate spatial predictors with spatialRF and MEM method (moran mapping)

# coordinates of the cases
xy <- dataset[, c("X", "Y")]
# distance matrix
distance.matrix <- as.matrix(dist(xy))
# distance thresholds (same units as distance_matrix)
# same units as projection!!!!! (degrees)
distance.thresholds <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)

mem <- mem_multithreshold(
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10)
dataset <- cbind(dataset, mem)

response_var <- "Total" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry = NULL,
                    nodesize = NULL, 
                    ntree = 1000)

print_results(results)

# spatial predictors but with PCA method

pca <- pca_multithreshold(
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  max.spatial.predictors = 10)
dataset <- cbind(dataset, pca)

write.csv2(dataset, 'data/modeling_data_v2.csv', row.names=F)
# dataset <- read.csv2('data/modeling_data_v2.csv')

response_var <- "Total" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry = NULL,
                    nodesize = NULL,
                    ntree = 1000)

print_results(results)

# total, best hyperparameter combination and no spatial predictors

dataset <- read.csv2('data/modeling_data_v3.csv')

response_var <- "Total" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry=6,
                    nodesize=60,
                    ntree = 1000)

print_results(results)

# logtotal, default hyperpar, no spatial predictors

dataset <- read.csv2('data/modeling_data_v3.csv')

response_var <- "logTotal" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry=NULL,
                    nodesize=NULL,
                    ntree = 1000)

print_results(results)

# logtotal, best hyperpar, no spatial predictors

dataset <- read.csv2('data/modeling_data.csv')

response_var <- "logTotal" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y, -linear_infrastr)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry = 6,
                    nodesize = 40,
                    ntree = 9000)

print_results(results)

# logTotal, best hyperpar, pca spatial predictors:

dataset <- read.csv2('data/modeling_data_v2.csv')

response_var <- "logTotal" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry = 6,
                    nodesize = 40,
                    ntree = 9000)

print_results(results)

# Total, best hyperpar, pca spatial predictors:

dataset <- read.csv2('data/modeling_data_v2.csv')

response_var <- "Total" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry = 6,
                    nodesize = 60,
                    ntree = 1000)

print_results(results)

# logTotal, best hyperpar FOR SPATIAL (PCA), pca spatial predictors:

dataset <- read.csv2('data/modeling_data_v2.csv')

response_var <- "logTotal" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry = 5,
                    nodesize = 50,
                    ntree = 1000)

print_results(results)
# best result

##################################################################################
########## SUPER FAST DEFAULT XGBOOST
#####################################################################
###############################################
################################


rm(list=ls())
source('scripts/utils.R')
library(xgboost)
library(caret)

train_xgb <- function(data, 
                      response_var, 
                      predictors, 
                      n_folds = 10, 
                      nrounds = 10000,  # XGBoost hyperparameters
                      eta = 0.01,      # Learning rate
                      max_depth = 6,  # Tree depth
                      subsample = 0.8,
                      colsample_bytree = 0.8) {   
  
  # Parallelize
  cl <- makeCluster(detectCores() - 1) 
  registerDoParallel(cl)
  
  # Construct formula
  formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
  
  # Folds loop
  results <- foreach(i = 1:n_folds,
                     .packages = c("xgboost", "dplyr", "caret"),
                     .combine = 'c', 
                     .multicombine = TRUE, 
                     .inorder = TRUE,
                     .verbose = TRUE) %dopar% {
                       
                       # Set unique seed for each fold
                       set.seed(i)
                       
                       # Use caret to randomly split train/test data (80-20)
                       index <- createDataPartition(data[[response_var]], p = 0.8, list = FALSE)
                       train <- data[index, ]
                       test <- data[-index, ]
                       
                       # Prepare data for XGBoost
                       dtrain <- xgb.DMatrix(data = as.matrix(train[predictors]), label = train[[response_var]])
                       dtest <- xgb.DMatrix(data = as.matrix(test[predictors]), label = test[[response_var]])
                       
                       # Train the model
                       model <- xgboost(data = dtrain, 
                                        nrounds = nrounds, 
                                        eta = eta, 
                                        max_depth = max_depth,
                                        subsample = subsample,
                                        colsample_bytree = colsample_bytree,
                                        objective = "reg:squarederror",
                                        verbose = 0)  # XGBoost doesn't need `importance='none'`
                       
                       # Generate training predictions
                       preds_train <- predict(model, newdata = dtrain)
                       train$pred <- preds_train
                       
                       # Calculate training metrics
                       mae_train <- MAE(preds_train, train[[response_var]])
                       mse_train <- mean((preds_train - train[[response_var]])^2)
                       rmse_train <- RMSE(pred = preds_train, obs = train[[response_var]])
                       Rsq_train <- 1 - (sum((train[[response_var]] - preds_train)^2) / 
                                           sum((train[[response_var]] - mean(train[[response_var]]))^2))
                       
                       # Generate testing predictions
                       preds_test <- predict(model, newdata = dtest)
                       test$pred <- preds_test
                       
                       # Calculate testing metrics
                       mae_test <- MAE(preds_test, test[[response_var]])
                       mse_test <- mean((preds_test - test[[response_var]])^2)
                       rmse_test <- RMSE(pred = preds_test, obs = test[[response_var]])
                       Rsq_test <- 1 - (sum((test[[response_var]] - preds_test)^2) / 
                                          sum((test[[response_var]] - mean(test[[response_var]]))^2))
                       
                       # Return results as a list
                       list(
                         metrics = c(fold = i, mae_train, mse_train, rmse_train, Rsq_train,
                                     mae_test, mse_test, rmse_test, Rsq_test),
                         train = train[[response_var]], preds_train = preds_train,
                         test = test[[response_var]], preds_test = preds_test
                       )
                     }
  
  stopCluster(cl)
  return(results)
}

# function to print results (mean metrics and preds~obs plots)
print_results <- function(results) {
  
  # extract metrics in sequence (5 elements per fold: metrics, train, preds_train, test, preds_test)
  metrics_list <- results[seq(1, length(results), by = 5)]
  
  # to df
  metrics_df <- do.call(rbind, metrics_list)
  
  # assign column names
  columnnames <- c("fold", "mae_train", "mse_train", "rmse_train", "Rsq_train", 
                   "mae_test", "mse_test", "rmse_test", "Rsq_test")
  colnames(metrics_df) <- columnnames
  
  # calculate average
  metrics_means <- colMeans(metrics_df)
  
  # print
  print(metrics_means[-1])
  
  # extract observations and predicted per fold for both train and test
  train_list <- results[seq(2, length(results), by = 5)]  
  preds_train_list <- results[seq(3, length(results), by = 5)] 
  
  test_list <- results[seq(4, length(results), by = 5)]  
  preds_test_list <- results[seq(5, length(results), by = 5)]  
  
  # plots list (x = preds, y = obs)
  plot_list_train <- lapply(1:10, function(i) {
    ggplot(data.frame(real = train_list[[i]], pred = preds_train_list[[i]]), aes(x = pred, y = real)) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste("Fold", i, "- Train"), x = "Predicted", y = "Real") +
      theme_minimal()
  })
  
  plot_list_test <- lapply(1:10, function(i) {
    ggplot(data.frame(real = test_list[[i]], pred = preds_test_list[[i]]), aes(x = pred, y = real)) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
      labs(title = paste("Fold", i, "- Test"), x = "Predicted", y = "Real") +
      theme_minimal()
  })
  
  # combine
  grid.arrange(grobs = c(plot_list_train, plot_list_test), ncol = 5)
}

# Total, default hyperpar, non-spatial

dataset <- read.csv2("data/modeling_data.csv")

response_var <- "Total" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_xgb(data = dataset, 
                     response_var = response_var, 
                     predictors = predictors)

print_results(results)

# much worse than RF

response_var <- "logTotal" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_xgb(data = dataset, 
                     response_var = response_var, 
                     predictors = predictors)

print_results(results)

# still worse performance in default than RF

dataset <- read.csv2("data/modeling_data_v2.csv")

response_var <- "Total" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_xgb(data = dataset, 
                     response_var = response_var, 
                     predictors = predictors)

print_results(results)

# as bad as non-spatial

response_var <- "logTotal" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_xgb(data = dataset, 
                     response_var = response_var, 
                     predictors = predictors)

print_results(results)

# better than non-spatial but still worse than default log+spatial RF