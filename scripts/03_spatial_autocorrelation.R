# we're trying some methods to correct spatial autocorrelation
# we're be using blockCV to conduct some of the methods

library(raster)
library(sf)
library(MazamaSpatialUtils)
library(blockCV)
library(mapview)
library(spdep)

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

# TEMPORAL BLOCK CROSS VALIDATION

