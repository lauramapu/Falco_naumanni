
# calculate interactions

dataset <- read.csv2("data/modeling_data.csv")

# Seleccionar las columnas de interés
variables <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) # Cambia por las variables que quieras usar

# Crear las interacciones usando model.matrix()
interaction_terms <- model.matrix(~ .^2, data = dataset[variables])

# Eliminar la columna del intercepto (si es necesario)
interaction_terms <- interaction_terms[, -1]

# Agregar las interacciones al dataframe original
dataset <- cbind(dataset[,1:5], interaction_terms)

response_var <- "logTotal" 
predictors <- colnames(dplyr::select(dataset, -Total, -logTotal, -Ano_CS, -X, -Y)) 

results <- train_rf(data = dataset, 
                    response_var = response_var, 
                    predictors = predictors,
                    mtry = 5,
                    nodesize = 50,
                    ntree = 10000)

print_results(results)
