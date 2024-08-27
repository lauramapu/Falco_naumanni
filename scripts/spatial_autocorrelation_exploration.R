# we're trying some methods to correct spatial autocorrelation
# we're be using blockCV to conduct some of the methods

library(raster)
library(sf)
library(MazamaSpatialUtils)
library(blockCV)
library(mapview)
library(spdep)
library(doParallel)
library(tictoc)

# load dataset
dataset <- read.csv2("data/modeling_data.csv")
# convert to sf
dataset_sf <- st_as_sf(dataset, coords=c('X','Y'), crs=4326)
dataset_sf <- cbind(dataset_sf, st_coordinates(dataset_sf))

# spanish provinces borders
provinces <- mapSpain::esp_get_prov()
# excluding Canarias and Baleares
provinces <- provinces[!provinces$iso2.prov.name.es %in% c("Las Palmas", "Santa Cruz de Tenerife", "Baleares"), ]
# dissolve
provinces$campo <- 1
mask_spain <- dissolve(provinces, field='campo')

# calculate moran autocorrelation index between observations

# select variable of interest
variable <- dataset_sf[dataset_sf$Total>0,]$Total
# generate matrix of spatial weights based on nearest neighbors
coords <- st_coordinates(dataset_sf[dataset_sf$Total>0,])
neighbors <- knearneigh(coords, k = 4) # k es el número de vecinos más cercanos
neighbors_list <- knn2nb(neighbors)
weights <- nb2listw(neighbors_list, style = "W")
# moran test
moran_result <- moran.test(variable, weights)
print(moran_result)
# although moran's index is low (0.16) results show super high spatial autocorrelation effect (p-value < 2.2e-16)

# since the study species is colonial and thus the spatial location of colonies is key in the distribution
# we're not gonna correct the effect of spatial autocorrelation, we're just adding it to our model as a predictor

# we're doing it with all points

# select variable of interest
variable <- dataset_sf$Total
# generate matrix of spatial weights based on nearest neighbors
coords <- st_coordinates(dataset_sf)
neighbors <- knearneigh(coords, k = 4) # k es el número de vecinos más cercanos
neighbors_list <- knn2nb(neighbors)
weights <- nb2listw(neighbors_list, style = "W")
# moran test
moran_result <- moran.test(variable, weights)
print(moran_result)

# calculate eigenvectors
weights_matrix <- listw2mat(weights)
# diagonal matrix in degress
D <- diag(rowSums(weights_matrix))
# Laplacian Normalized Matrix 
L <- D - weights_matrix
# calculate Laplacian Normalized Matriz (simetric)
D_inv_sqrt <- diag(1 / sqrt(diag(D)))
L_normalized <- D_inv_sqrt %*% L %*% D_inv_sqrt
# calculate eigenvalues and eigenvectors
eigen_result <- eigen(L_normalized)
# extract eigenvectors
eigenvectors <- eigen_result$vectors
# extract eigenvalues
eigenvalues <- eigen_result$values

# erase the imaginary part because it's insignificant and due to computation error with large float numbers
eigenvectors_real <- Re(eigenvectors)
eigenvalues_real <- Re(eigenvalues)

# now we need to choose a reasonable number of eigenvectors basing on eigenvalues
plot(eigenvalues_real, type="b", main="Codo de Eigenvalores", xlab="Índice del Eigenvector", ylab="Eigenvalor")
# it seems there's no reasonable number but we'll try to calculate with a 90% threshold

# total variance
total_variance <- sum(eigenvalues_real)
# accumulated explained variance
acc_variance <- Re(cumsum(eigenvalues_real) / total_variance)
# choose a threshold (percent of absorbed variance)
threshold <- 0.90
# find minimum number of eigenvectors that explain our threshold of explained variance
n <- which(acc_variance >= threshold)[1]

# get first n eigenvectors
eigenvectors_signif <- eigenvectors[, 1:n]

# k-means clustering
# need to transpose because kmeans groups by obs (rows) and eigenvectors are columns
kmeans_result <- kmeans(t(Re(eigenvectors_signif)), centers = 10)

# show clusters and frecs
table(kmeans_result$cluster)
# 3221 vectors are grouped into '4', but some must be super rare since they're classified all alone

View(kmeans_result$centers)
# we're adding centers as variables

kmeans_centers <- as.data.frame(scale(t(kmeans_result$centers)))
colnames(kmeans_centers) <- paste('eigen', colnames(kmeans_centers), sep = '')

dataset_spatial <- cbind(dataset, kmeans_centers)
write.csv2(dataset_spatial, "data/modeling_data_spatial.csv", row.names=F)

# TEMPORAL BLOCK CROSS VALIDATION
# we split data per years, so that we have 3 blocks to cross validate (folds)
# also year of census is also related to the spatial location of points because in 2017 and 2018
# only a few locations were visited

rm(list=ls())

# activate parallel process
# n.cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK")) # slurm cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

dataset <- read.csv2("data/modeling_data_spatial.csv")
dataset <- dataset %>% dplyr::select(-geometry)

tic("temporal cross validation")
metrics <- foreach(i = 1:length(unique(dataset$Ano_CS)),
           .packages=c("randomForest", "sf", "dplyr", "caret"),
           .combine='rbind') %dopar% { # iterations through folds
    
            # set unique seed for each fold
            set.seed(i+20)
            
            # set year of test
            year <- unique(dataset$Ano_CS)[i]
            
            # use caret to randomly split train/test data (80-20)
            test <- subset(dataset, Ano_CS == year)
            train <- subset(dataset, Ano_CS != year)
            
            # extracting sample size
            prNum <- nrow(train) - as.numeric(table(train$Total)['0']) # number of presences
            
            # train the model
            model <- randomForest(Total ~ . -X -Y -Ano_CS,
                                  data = train,
                                  sampsize = prNum,
                                  mtry = as.integer(ncol(train)*0.1),
                                  ntree = 1000,
                                  nodesize = as.integer(nrow(train)*0.1))
            
            # generate training predictions (only select probs of presence)
            preds <- predict(model)
            train$pred <- preds
            
            # calculate metrics (caret)
            mae_train <- MAE(preds, train$Total)
            mse_train <- mean((preds - train$Total)^2)
            rmse_train <- RMSE(preds, train$Total)
            
            # generate testing predictions and binaries from maxSSS
            preds <- predict(model, newdata = test)
            test$pred <- preds
            
            # calculate metrics (caret)
            mae_test <- MAE(preds, test$Total)
            mse_test <- mean((preds - test$Total)^2)
            rmse_test <- RMSE(preds, test$Total)
            
            # store all the relevant values in the metrics dataframe
            
            metrics <- c(year, mae_train, mse_train, rmse_train, mae_test, mse_test, rmse_test)
            return(metrics)
            }
toc() 

# assign col names
metrics <- as.data.frame(metrics)
column_names <- c('year', 'MAE_train', 'MSE_train', 'RMSE_train',
                  'MAE_test', 'MSE_test', 'RMSE_test')
colnames(metrics) <- column_names
write.csv2(metrics, "results/hyperparameters.csv", row.names = F)
# metrics <- read.csv2("results/hyperparameters.csv")
