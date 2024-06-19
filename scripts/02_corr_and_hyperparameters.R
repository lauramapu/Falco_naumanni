# hyperparameter tuning

library(ranger)
library(caret)
library(dplyr)
library(tictoc)
library(doParallel)

# hyperparameter tuning on parallel

# we need to tune mtry, nodesize and ntree
# we need to cross validate to ensure all data is being used in test
# we're using random cross-validation in 5 folds

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

x <- foreach(
  i = 1:10, 
  .combine = 'c'
) %dopar% {
  sqrt(i)
}
x # check parallel execution is on

dataset <- read.csv2("data/modeling_data.csv")

# set grid to search hyperparameters
# define ranges
set.seed(123)
grid <- expand.grid(mtry = 2:ncol(dataset[,5:44])/2, # variables in each tree
                    nodesize = seq(10, nrow(dataset[dataset$Total>0,]), by = 10)) # minimum size of terminal nodes
write.csv2(grid, "data/grid_hyperpar.csv", row.names=F)
# grid <- read.csv2("data/grid_hyperpar.csv")

# run models iterating through each hyperpar combination
tic("random hyperparameters")
metrics <- foreach(i = 1:nrow(grid),
                   .packages=c("ranger", "sf", "dplyr", "caret"),
                   .combine='rbind') %:% # iterations through hyperpar combinations
  
  foreach(j = 1:5,
          .combine='rbind',
          verbose=T) %dopar% { # iterations through random cv folds (5)
            
            # set unique seed for each fold
            set.seed(j**2)
            
            # use caret to randomly split train/test data (80-20)
            index <- createDataPartition(dataset$Total, p = 0.8, list = FALSE)
            train <- dataset[index, ]
            test <- dataset[-index, ]  
            
            # calculate case weights in training
            # class frecs
            freq_zero <- sum(train$Total == 0)
            freq_nonzero <- sum(train$Total != 0)
            # calculate if frec zero > non_zero
            if (freq_zero > freq_nonzero) {
              # inversaly proportional weights to the frecs
              weight_zero <- 1 / freq_zero
              weight_nonzero <- 1 / freq_nonzero
              
              # Asignar los pesos a cada observación
              train <- train %>%
                mutate(case_weight = ifelse(response == 0, weight_zero, weight_nonzero))
              
              # normalize weights to sum up 1
              train$case_weight <- train$case_weight / sum(train$case_weight)
            } else {
              # assign same weight to all obs
              train <- train %>%
                mutate(case_weight = 1 / n())
            }
            
            # train the model
            model <- ranger(Total ~ . -X -Y -geometry -Ano_CS -case_weight,
                            data = train,
                            case.weights = train$case_weight,
                            mtry = grid[i,1],
                            num.tree = 1000,
                            min.node.size = grid[i,2],
                            importance='none')
            
            # generate training predictions (only select probs of presence)
            preds <- predict(model, data=train)
            train$pred <- preds[[1]]
            
            # calculate metrics (caret)
            mae_train <- MAE(preds, train$Total)
            mse_train <- mean((preds - train$Total)^2)
            rmse_train <- RMSE(preds, train$Total)
            Rsq_train <- 1 - (sum((train$Total - preds)^2) / sum((train$Total - mean(train$Total))^2))
            
            # generate testing predictions and binaries from maxSSS
            preds <- predict(model, data=test)
            test$pred <- preds[[1]]
            
            # calculate metrics (caret)
            mae_test <- MAE(preds, test$Total)
            mse_test <- mean((preds - test$Total)^2)
            rmse_test <- RMSE(preds, test$Total)
            Rsq_test <- 1 - (sum((test$Total - preds)^2) / sum((test$Total - mean(test$Total))^2))
            
            # store all the relevant values in the metrics dataframe
            # column_names <- c('mtry', 'ntree', 'nodesize','thr_value',
            # 'Sens_train', 'F1_train','B.Accuracy_train','TSS_train','AUC_train',
            # 'Sens_test', 'F1_test', 'B.Accuracy_test','TSS_test','AUC_test')
            
            metrics <- c(i, j, grid[i,1], grid[i,2], mae_train, mse_train, rmse_train, Rsq_train,
                         mae_test, mse_test, rmse_test, Rsq_test)
            return(metrics)
          }
toc() 

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

# best combination should be the one that gets higher sensitivity
# and lower difference between training and testing sensitivity
# because if test sens < train sens = overfitting
# if test sens > train sens = underfitting
# equal metrics mean generalization, which is our goal

metrics_mean$sum <- rowSums(metrics_mean[,6:11])
metrics_mean[which.min(metrics_mean$sum),]

metrics_mean$sens_diff <- abs(metrics_mean$Sens_test-metrics_mean$Sens_train)
metrics_mean[which.min(metrics_mean$sens_diff),]

# we choose nodesize=90 and mtry=6, since it's the combination with lowest sens diff (0)
# and still high metrics sum (7.62, max=7.96)
# this combination forms the most complex trees in terms of depth and still maximizes generalization
# better combinations obtain the same metrics from 1000 to 10000 trees, so we choose 1000
