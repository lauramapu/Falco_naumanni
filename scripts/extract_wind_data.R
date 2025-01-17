
rm(list=ls())
source('scripts/utils.R')

spain <- getpeninsularspain()

gpkg_files <- list.files("../../SDMs-Conflicts/spatial_data/siose_ar", pattern = "\\.gpkg$",
                         full.names = TRUE, recursive = TRUE)
# with recursive=T it searches in all subfolders

# extract the USOS layer
extract_wind <- function(file_path) {
  # read all layers
  layers <- st_layers(file_path)
  # find the one that matches _USOS in its name
  layer_name <- layers$name[grep("_T_USOS$", layers$name)]
  if (length(layer_name) == 0) {
    # if not found return a warning
    stop("'_T_USOS' layer not found")
  }
  # read layer
  layer <- st_read(file_path, layer = layer_name[1])
  # select only solar and wind polygons
  layer <- layer[layer$ID_USO_MAX %in% c(2441), ]
  return(layer)
}

n.cores <- detectCores() - 1  # Puedes usar todos los núcleos menos uno para no saturar tu máquina
cl <- makeCluster(n.cores)
clusterExport(cl, varlist = c("extract_wind", "gpkg_files"))
clusterEvalQ(cl, library(sf))

# erase baleares [6], las palmas [30], tenerife [32], ceuta [45],
# melilla [46] and posesiones [53]
capas <- parLapply(cl, gpkg_files[-c(6, 30, 32, 45, 46, 53)], extract_wind)

# project all to LAEA
capas2 <- parLapply(cl, capas, function(x) {st_transform(x, crs = 3035)})

stopCluster(cl)

# join layers

# initialize object to store layers
sf_merged <- capas2[[1]]
# iterate through each layer
for (i in 2:length(capas2)) {
  if (ncol(capas2[[i]]) > 17) {
    capas2[[i]] <- capas2[[i]][,-1] # erase the id column
    sf_merged <- rbind(sf_merged, capas2[[i]]) # add layer
  } else {  # ncol = 17 ok
    sf_merged <- rbind(sf_merged, capas2[[i]])
  }
}
mapview(sf_merged)

write_sf(sf_merged, "spatial_data/energy_infrastr/wind_turbines/siose_wind.shp")


