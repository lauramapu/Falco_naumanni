# hyperparameter tuning

rm(list=ls())
source('scripts/utils.R')

# hyperparameter tuning on parallel

# we need to tune mtry, nodesize and ntree
# we need to cross validate to ensure all data is being used in test
# we're using random cross-validation in 5 folds

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)

dataset <- read.csv2("data/modeling_data.csv")
dataset[ , sapply(dataset, is.numeric)] <- lapply(dataset[ , sapply(dataset, is.numeric)], as.numeric)

# set grid to search hyperparameters
# define ranges
grid <- expand.grid(mtry = 2:ncol(dplyr::select(dataset, -Total, -Ano_CS, -logTotal, -X, -Y)), # variables in each tree
                    nodesize = c(1,2,3,4,5,6,7,8,9, seq(10, nrow(dataset), by = 10))) # minimum size of terminal nodes
write.csv2(grid, "data/grid_hyperpar.csv", row.names=F)
# grid <- read.csv2("data/grid_hyperpar.csv")

# run models iterating through each hyperpar combination
tic("mtry and nodesize tuning")
metrics <- foreach(i = 1:nrow(grid),
                   .packages=c("ranger", "dplyr", "caret"),
                   .combine='rbind') %:% # iterations through hyperpar combinations
  
           foreach(j = 1:10,
           .combine='rbind') %dopar% { # iterations through random cv folds (5)
            
            # set unique seed for each fold
            set.seed(j**2)
            
            # use caret to randomly split train/test data (80-20)
            index <- createDataPartition(dataset$Total, p = 0.8, list = FALSE)
            train <- dataset[index, ]
            test <- dataset[-index, ]  
            
            # train the model
            model <- ranger(logTotal ~ . -Total -X -Y -Ano_CS,
                            data = train,
                            mtry = grid[i,1],
                            num.tree = 1000,
                            min.node.size = grid[i,2],
                            importance='none')
            
            # generate training predictions 
            preds <- predict(model, data=train)
            train$pred <- preds[[1]]
            
            # calculate metrics (caret)
            mae_train <- MAE(preds[[1]], train$logTotal)
            mse_train <- mean((preds[[1]] - train$logTotal)^2)
            rmse_train <- RMSE(pred = preds[[1]], obs = train$logTotal)
            Rsq_train <- 1 - (sum((train$logTotal - preds[[1]])^2) / sum((train$logTotal - mean(train$logTotal))^2))
            
            # generate testing predictions 
            preds <- predict(model, data=test)
            test$pred <- preds[[1]]
            
            # calculate metrics (caret)
            mae_test <- MAE(preds[[1]], test$logTotal)
            mse_test <- mean((preds[[1]] - test$logTotal)^2)
            rmse_test <- RMSE(pred = preds[[1]], obs = test$logTotal)
            Rsq_test <- 1 - (sum((test$logTotal - preds[[1]])^2) / sum((test$logTotal - mean(test$logTotal))^2))
            
            # store all the relevant values in the metrics dataframe
            metrics <- c(i, j, grid[i,1], grid[i,2], mae_train, mse_train, rmse_train, Rsq_train,
                         mae_test, mse_test, rmse_test, Rsq_test)
            return(metrics)
          }
toc() # 12077.67 sec

# assign col names
metrics <- as.data.frame(metrics)
column_names <- c('combination', 'fold', 'mtry', 'nodesize', 'MAE_train', 'MSE_train', 'RMSE_train', 'Rsquared_train',
                  'MAE_test', 'MSE_test', 'RMSE_test', 'Rsquared_test')
colnames(metrics) <- column_names
write.csv2(metrics, "results/hyperparameters.csv", row.names = F)
# metrics <- read.csv2("results/hyperparameters.csv")

metrics_mean <- metrics %>%
  group_by(combination) %>%
  summarise(across(-fold, mean))
write.csv2(metrics_mean, "results/hyperparameters_mean.csv", row.names = F)
# metrics_mean <- read.csv2("results/hyperparameters_mean.csv")

# best combination should be the one that gets lower metrics sum (error measures)
# and lower difference between training and testing metrics
# because if test < train  = overfitting
# if test > train = underfitting
# equal metrics mean generalization, which is our goal

metrics_mean$sum <- rowSums(metrics_mean[,c('Rsquared_train', 'Rsquared_test')])
metrics_mean[which.max(metrics_mean$sum),]

metrics_mean$Rsquared_diff <- abs(metrics_mean$Rsquared_test-metrics_mean$Rsquared_train)
metrics_mean[which.min(metrics_mean$Rsquared_diff),]

plot(metrics_mean$Rsquared_diff, metrics_mean$sum)

# scatterplot of both derived metrics
p <- ggplot(metrics_mean, aes(x = sum, y = Rsquared_diff)) +
  geom_point(aes(text = paste("Sum:", sum, "<br>Rsquared_diff:", Rsquared_diff, "<br>mtry:", mtry, "<br>nodesize:", nodesize))) +
  labs(title = "Scatterplot Interactivo", x = "Sum", y = "Rsquared_diff")

# interactive plot
p_interactive <- ggplotly(p, tooltip = "text")
p_interactive

plot(metrics_mean$mtry, metrics_mean$Rsquared_test) # better the higher
plot(metrics_mean$nodesize, metrics_mean$Rsquared_test) # lower

plot(metrics_mean$Rsquared_train, metrics_mean$Rsquared_test) # better the higher

# since results are overall very good we're just choosing the comb with highest Rsq in test
# generalization in this comb is very good (loss of 0.02)
metrics_mean[which.max(metrics_mean$Rsquared_test),]
# mtry = 4; nodesize = 20
# still Rsq test is quite low (~0.1)

# now we need to check 1000 trees are actually enough

ntrees <- seq(1000, 10000, by = 500)

# run models iterating through each hyperpar combination
tic("ntrees tuning")
metrics <- foreach(i = 1:length(ntrees),
                   .packages=c("ranger", "dplyr", "caret"),
                   .combine='rbind') %:% # iterations through trees number
  
  foreach(j = 1:10,
          .combine='rbind') %dopar% { # iterations through random cv folds (10)
            
            # set unique seed for each fold
            set.seed(j**2)
            
            # use caret to randomly split train/test data (80-20)
            index <- createDataPartition(dataset$logTotal, p = 0.8, list = FALSE)
            train <- dataset[index, ]
            test <- dataset[-index, ]  
            
            # train the model
            model <- ranger(logTotal ~ . -X -Y -Ano_CS -Total,
                            data = train,
                            num.tree = ntrees[[i]],
                            mtry = metrics_mean$mtry[which.max(metrics_mean$Rsquared_test)],
                            min.node.size = metrics_mean$nodesize[which.max(metrics_mean$Rsquared_test)],
                            importance='none')
            
            # generate training predictions 
            preds <- predict(model, data=train)
            train$pred <- preds[[1]]
            
            # calculate metrics (caret)
            mae_train <- MAE(preds[[1]], train$logTotal)
            mse_train <- mean((preds[[1]] - train$logTotal)^2)
            rmse_train <- RMSE(pred = preds[[1]], obs = train$logTotal)
            Rsq_train <- 1 - (sum((train$logTotal - preds[[1]])^2) / sum((train$logTotal - mean(train$logTotal))^2))
            
            # generate testing predictions a
            preds <- predict(model, data=test)
            test$pred <- preds[[1]]
            
            # calculate metrics (caret)
            mae_test <- MAE(preds[[1]], test$logTotal)
            mse_test <- mean((preds[[1]] - test$logTotal)^2)
            rmse_test <- RMSE(pred = preds[[1]], obs = test$logTotal)
            Rsq_test <- 1 - (sum((test$logTotal - preds[[1]])^2) / sum((test$logTotal - mean(test$logTotal))^2))
            
            # store all the relevant values in the metrics dataframe
            metrics <- c(ntrees[i], j, mae_train, mse_train, rmse_train, Rsq_train,
                         mae_test, mse_test, rmse_test, Rsq_test)
            return(metrics)
          }
toc() # 341.59 sec

# assign col names
metrics <- as.data.frame(metrics)
column_names <- c('ntrees', 'fold', 'MAE_train', 'MSE_train', 'RMSE_train', 'Rsquared_train',
                  'MAE_test', 'MSE_test', 'RMSE_test', 'Rsquared_test')
colnames(metrics) <- column_names
write.csv2(metrics, "results/hyperparameters_ntrees.csv", row.names = F)

metrics_mean <- metrics %>%
  group_by(ntrees) %>%
  summarise(across(-fold, mean))
write.csv2(metrics_mean, "results/hyperparameters_ntrees_mean.csv", row.names = F)

plot(metrics_mean$ntrees, metrics_mean$Rsquared_test)
metrics_mean[which.max(metrics_mean$Rsquared_test),]

# highest Rsq in test at 1000 trees

registerDoSEQ()
