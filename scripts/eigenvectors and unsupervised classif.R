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

# Obtener los primeros n eigenvectores
eigenvectores_significativos <- eigenvectors[, 1:n]

# Imprimir el número de eigenvectores seleccionados y los eigenvalores correspondientes
print(n)
print(eigenvalues[1:n])

# now we now how many vectors we have to use to absorb 90% of the variance
# next step is simplifying these vectors with a PCA

# Realizar el análisis de componentes principales
pca <- prcomp(eigenvectors_real, scale. = TRUE)

# Verificar la varianza explicada por cada componente
varianza_explicada <- pca$sdev^2
cum_var <- cumsum(varianza_explicada)
cum_var_pct <- cum_var / sum(varianza_explicada)

# Encontrar el número de componentes que explican al menos el 90% de la varianza
n_components <- which(cum_var_pct >= 0.90)[1]

# number of needed components to explain 90% of total variance is 494 so this won't help us
# we're conducting a kmeans clustering to extract only 10 vectors from all the eigenvectors

# Realizar el clustering K-Means
kmeans_result <- kmeans(eigenvectors_real, centers = 10)

# Mostrar los clusters
table(kmeans_result$cluster)
# 3221 vectors are grouped into '6', but some must be super rare since they're classified all alone

# Agrupar los eigenvectores según los clusters
eigenvectors_clusters <- tapply(as.data.frame(pca_selected), kmeans_result$cluster, function(x) colSums(x))
# Mostrar los eigenvectores agrupados
eigenvectors_clusters

# now we need to incorporate these data as variables in our model
# we're using the distances approach, in which we calculate the distances of each observation to its 
# corresponding cluster center

# Función para calcular la distancia euclidiana
calc_distance <- function(row, center) {
  sqrt(sum((row - center)^2))
}

# Asegúrate de que los nombres de las columnas coincidan entre eigenvectors_real y kmeans_result$centers
colnames(kmeans_result$centers) <- colnames(eigenvectors_real)

# Añadir las distancias al centro asignado como nueva columna en dataset_sf
dataset_sf$DistToAssignedCluster <- NA

for (i in 1:nrow(eigenvectors_real)) {
  cluster <- kmeans_result$cluster[i]
  center <- kmeans_result$centers[cluster, ]
  # Calcular la distancia al centro asignado
  distance <- calc_distance(eigenvectors_real[i, ], center)
  # Añadir la distancia al data frame original dataset_sf
  dataset_sf$DistToAssignedCluster[i] <- distance
}

# Verificar si hay NA en la nueva columna
sum(is.na(dataset_sf$DistToAssignedCluster))

# Mostrar el resultado final
print(head(dataset_sf))

dataset_sf$AssignedCluster <- kmeans_result$cluster

########################
# MORAN EIGENVECTORS SPATIAL FILTERING (MESF)

# create list of neighbors using all coords
coords <- as.matrix(dataset[, c("X", "Y")])
nb <- dnearneigh(coords, d1 = 0, d2 = max(dist(coords)))

# Crear una lista de pesos espaciales
listw <- nb2listw(nb, style="W")

# Definir la fórmula del modelo
formula <- as.formula(paste("Total ~", paste(names(dataset)[5:44], collapse=" + ")))

# Aplicar la función ME
me_model <- ME(formula, data = dataset, family = gaussian, listw = listw, alpha = 0.05, nsim = 99, verbose = TRUE)

# Inspeccionar el modelo ajustado
summary(me_model)

# Ver los autovectores seleccionados
me_model$fit
