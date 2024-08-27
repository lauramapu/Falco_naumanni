# fitting process with spatialrf (Blas Benito)

rm(list=ls())

library(spatialRF)
library(dplyr)

# load data
dataset <- read.csv2("data/modeling_data.csv")
dataset <- dataset %>%
  rename(x=X, y=Y)
dataset[ , sapply(dataset, is.numeric)] <- lapply(dataset[ , sapply(dataset, is.numeric)], as.numeric)

# # vamos a hacer un apaño cutre para no tener coordenadas duplicadas hasta que sepa que hacer con los duplicados
# # sumamos 0.001 a la x y la y de las coords duplicadas
# dataset <- dataset %>%
#   group_by(x, y) %>%
#   mutate(
#     row_id = row_number(),  # Crear un identificador de fila para cada grupo de duplicados
#     x = if_else(row_id == 1, x, x + 0.1 * (row_id - 1)),  # Sumar 0.001 por cada duplicado adicional
#     y = if_else(row_id == 1, y, y + 0.1 * (row_id - 1))
#   ) %>%
#   ungroup() %>%
#   dplyr::select(-row_id) # Eliminar la columna de identificador de fila
  
any(duplicated(dataset[,c('x','y')])) # check

# names of the response variable and the predictors
dependent.variable.name <- "logTotal"
predictor.variable.names <- colnames(dataset)[6:49]

# coordinates of the cases
xy <- dataset[, c("x", "y")]

# distance matrix
distance.matrix <- as.matrix(dist(xy))

# distance thresholds (same units as distance_matrix)
# same units as projection!!!!! (degrees)
distance.thresholds <- c(0.05, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4, 0.5, 1)

# random seed for reproducibility
random.seed <- 21

# plot scatterplots of predictors
spatialRF::plot_training_df(
  data = dataset,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  ncol = 10,
  point.color = viridis::viridis(100, option = "F"),
  line.color = "gray30"
)

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

# interactions <- spatialRF::the_feature_engineer(
#   data = dataset,
#   dependent.variable.name = dependent.variable.name,
#   predictor.variable.names = predictor.variable.names,
#   xy = xy,
#   importance.threshold = 0.50, # uses 50% best predictors
#   cor.threshold = 0.75, # max corr between interactions and predictors
#   seed = random.seed,
#   repetitions = 100,
#   verbose = TRUE
# )
# # esto tarda un huevo, 190 posibles interacciones

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
# still correlation at low distances but it says there's no corr between residuals

spatialRF::plot_importance(
  model.spatial, 
  verbose = FALSE) + 
  ggplot2::ggtitle("Non-spatial model") 

kableExtra::kbl(
  head(model.spatial$importance$per.variable, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

# since there's no spatial correlation we continue with no spatial predictors
