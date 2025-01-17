library(readxl)
library(dplyr)
library(mapview)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggspatial) 
library(ggrepel)
library(sf)
library(units)
library(mapSpain)
library(cowplot)   # combine ggplot2 objects
library(gridExtra)
library(raster)
library(MazamaSpatialUtils)
library(biscale)
library(grDevices) # to create custom color palette
library(biscale)
library(terra)
library(purrr)
library(caret)
library(pdp)
library(dplyr)
library(reshape2) # to melt data for graphics
library(tictoc)
library(pROC) # roc curve, calculate thresholds
library(tidyr)
library(modEvA) # boyce index
library(car) # vif
library(spatialRF) # spatial autocorrelation tests
library(corrplot)
library(ranger)
library(plotly)
library(doParallel)
library(landscapemetrics)

# function to list and cite all the above packages
versionandcitation <- function() {
  # list
  packages <- c("dplyr", "mapview", "ggplot2", "RColorBrewer", "viridis",
                "ggspatial", "ggrepel", "sf", "units", "mapSpain", 
                "cowplot", "gridExtra", "raster", "MazamaSpatialUtils", 
                "biscale", "grDevices", "terra", "purrr", "randomForest", 
                "caret", "pdp", "reshape2", "tictoc", "pROC", "tidyr")
  
  # format the reference
  format_citation <- function(citation_info) {
    title <- citation_info$title
    authors <- sapply(citation_info$author, function(a) paste(a$given, a$family))
    year <- citation_info$year
    note <- citation_info$note
    url <- citation_info$url
    
    citation_text <- paste0(title, ". ", paste(authors, collapse = ", "), ". ", year, ". ", note, ". ", url)
    return(citation_text)
  }
  
  # versions and references
  package_info <- lapply(packages, function(pkg) {
    version <- as.character(packageVersion(pkg))
    citation_info <- citation(pkg)
    citation_text <- format_citation(citation_info)
    list(Package = pkg, Version = version, Citation = citation_text)
  })
  
  # to df
  package_info_df <- do.call(rbind, lapply(package_info, as.data.frame))
  rownames(package_info_df) <- NULL
  
  # print
  print(package_info_df)
  
  # save to txt
  write.table(package_info_df, file = "R-packages.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# function to load spanish borders
getpeninsularspain <- function() {
  provinces <- mapSpain::esp_get_prov()[!mapSpain::esp_get_prov()$iso2.prov.name.es %in%
                                          c("Las Palmas", "Santa Cruz de Tenerife", "Baleares", "Melilla", "Ceuta"), ]
  provinces$campo <- 1
  spain <- provinces %>%
    dissolve(field='campo') %>%
    st_transform(crs=4326)
}

# crop+mask function for all variables
mask_raster <- function(raster_layer) {
  maskk <- st_transform(spain, st_crs(raster_layer))
  cropped <- raster::crop(raster_layer, extent(maskk))
  masked <- raster::mask(cropped, maskk)
  return(masked)
}

# function to construct bivariated maps
bimap <- function(biclass_object, title, xlab, ylab, mypal) {
  # the plot itself
  p <- ggplot() +
    geom_raster(data = biclass_object, aes(x = x, y = y, fill = bi_class)) +
    bi_scale_fill(pal = mypal, dim=3, na.value="transparent") +
    labs(x = "Longitude", y = "Latitude", fill = "") +
    ggtitle(title) +
    theme_bw() + 
    theme( # legend modifications
      legend.position = 'none',  
      legend.background = element_rect(fill = "transparent"),  # transparent background
      legend.key = element_rect(fill = "transparent", color = NA)  # transparent background for legend keys
    ) +
    geom_sf(data = spain, fill = "transparent", color = "black", inherit.aes=F)
  p
  
  # bivariate legend 
  legend <- bi_legend(pal = mypal,
                      dim = 3,
                      xlab = xlab,
                      ylab = ylab,
                      size = 10,
                      arrows = T,
                      pad_width = 1.5) # distance between colors
  legend <- legend +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      axis.title.x = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 10, color = "black", angle = 90), # vertical text
      axis.text = element_blank(),  # hide exes values
      axis.ticks = element_blank()  # hide exes values
    )
  legend
  
  # combine plot and legend in same plot
  combined_plot <- ggdraw() +
    draw_plot(p) +
    draw_plot(legend, x = 0.68, y = 0.11, width = 0.25, height = 0.25)
  return(combined_plot)
}

# function to find most common value in row with sapply, else calculate mean integer
mayor <- function(x) {
  x <- na.omit(x)
  
  # if all are NA, then NA
  if (length(x) == 0) {
    return(NA)
  }
  # calculate frecs
  freq_table <- table(x)
  # extract max value
  max_value <- which.max(freq_table)
  # check if max value is unique
  if (length(which(freq_table == max(freq_table))) == 1) {
    value <- names(freq_table)[max_value]
    return(as.numeric(value))
  } else {
    # if it's not unique, calculate mean of all them and return nearest integer
    numeric_values <- as.numeric(as.character(x))
    mean_value <- mean(numeric_values)
    rounded_value <- round(mean_value)
    return(rounded_value)
  }
}

solar_overlap <- function(sp_distribution) {
  nsolar <- table(topcells$solar)[[2]]
  superposition_area <- sum(topcells$solar == 1 & topcells[,sp_distribution] == 10, na.rm = TRUE)
  superposition_percent <- sum(topcells$solar == 1 & topcells[,sp_distribution] == 10, na.rm = TRUE) / nsolar *100
  vector_final <- c(superposition_percent, superposition_area)
}

# custom color palette for biscale
custom_bipal <- function(min.xy, max.y, max.x, max.xy, dim) {
  
  # Convierte los colores a valores RGB
  tl <- col2rgb(min.xy)
  tr <- col2rgb(max.y)
  bl <- col2rgb(max.x)
  br <- col2rgb(max.xy)
  
  # Función para interpolar colores
  interpolate_color <- function(tl, tr, bl, br, x, y) {
    top <- (1 - x) * tl + x * tr
    bottom <- (1 - x) * bl + x * br
    color <- (1 - y) * top + y * bottom
    return(rgb(color[1], color[2], color[3], maxColorValue = 255))
  }
  
  # Inicializa la matriz para almacenar los valores de color
  color_matrix <- matrix(NA, nrow = dim, ncol = dim)
  
  # Genera la cuadrícula de colores
  for (i in seq_len(dim)) {
    for (j in seq_len(dim)) {
      x <- (i - 1) / (dim - 1)
      y <- (j - 1) / (dim - 1)
      color_matrix[i, j] <- interpolate_color(tl, tr, bl, br, x, y)
    }
  }
  
  # Genera el vector de nombres basado en las dimensiones
  custom_pal_names <- as.vector(outer(1:dim, 1:dim, function(x, y) paste(x, y, sep = "-")))
  
  # Convierte la matriz de colores en un vector nombrado
  color_vector <- as.vector(t(color_matrix))
  names(color_vector) <- custom_pal_names
  
  return(color_vector)
}

# function to train ranger model with 10 folds and random CV
train_rf <- function(data, 
                     response_var, 
                     predictors, 
                     n_folds = 10, 
                     ntree,
                     mtry,            
                     nodesize) {   
  
  # paralellize
  cl <- makeCluster(detectCores() - 1) 
  registerDoParallel(cl)
  
  # construct formula
  formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
  
  # folds loop
  results <- foreach(i = 1:n_folds,
                     .packages = c("ranger", "dplyr", "caret"),
                     .combine = 'c', 
                     .multicombine = TRUE, 
                     .inorder = TRUE,
                     .verbose = TRUE) %dopar% {
                       
                       # set unique seed for each fold
                       set.seed(i)
                       
                       # use caret to randomly split train/test data (80-20)
                       index <- createDataPartition(data[[response_var]], p = 0.8, list = FALSE)
                       train <- data[index, ]
                       test <- data[-index, ]
                       
                       # train the model
                       model <- ranger(formula = formula,
                                       data = train,
                                       mtry = mtry,
                                       min.node.size = nodesize,
                                       ntree = ntree, 
                                       importance='none')
                       
                       # generate training predictions
                       preds_train <- predict(model, data = train)
                       train$pred <- preds_train$predictions
                       
                       # calculate training metrics
                       mae_train <- MAE(preds_train$predictions, train[[response_var]])
                       mse_train <- mean((preds_train$predictions - train[[response_var]])^2)
                       rmse_train <- RMSE(pred = preds_train$predictions, obs = train[[response_var]])
                       Rsq_train <- 1 - (sum((train[[response_var]] - preds_train$predictions)^2) / 
                                           sum((train[[response_var]] - mean(train[[response_var]]))^2))
                       
                       # generate testing predictions
                       preds_test <- predict(model, data = test)
                       test$pred <- preds_test$predictions
                       
                       # calculate testing metrics
                       mae_test <- MAE(preds_test$predictions, test[[response_var]])
                       mse_test <- mean((preds_test$predictions - test[[response_var]])^2)
                       rmse_test <- RMSE(pred = preds_test$predictions, obs = test[[response_var]])
                       Rsq_test <- 1 - (sum((test[[response_var]] - preds_test$predictions)^2) / 
                                          sum((test[[response_var]] - mean(test[[response_var]]))^2))
                       
                       # Return results as a list
                       list(
                         metrics = c(fold = i, mae_train, mse_train, rmse_train, Rsq_train,
                                     mae_test, mse_test, rmse_test, Rsq_test),
                         train = train[[response_var]], preds_train = preds_train$predictions,
                         test = test[[response_var]], preds_test = preds_test$predictions
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