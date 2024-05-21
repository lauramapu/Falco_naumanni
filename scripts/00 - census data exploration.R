# presence points and environmental variables

library(readxl)
library(sf)
library(mapview)
library(mapSpain)
library(MazamaSpatialUtils)
library(raster)
library(dplyr)

# read census data
census <- read_xls("data/primilla_Nacional_2016-2018_V5.xls", sheet = 'CENSO NACIONAL')
View(census)

# total individuals is called Total_CS
presence <- subset(census, census$Total_CS > 0)
absence <- subset(census, census$Total_CS == 0)

# reproyect to 4326 from corresponding huso (29, 30 and 31)
to_point <- function(x) {
  points29 <- st_transform(st_as_sf(subset(x, x$Huso == 29), coords = c('Coord.X', 'Coord.Y'), crs = 25829), crs = 4326)
  points30 <- st_transform(st_as_sf(subset(x, x$Huso == 30), coords = c('Coord.X', 'Coord.Y'), crs = 25830), crs = 4326)
  points31 <- st_transform(st_as_sf(subset(x, x$Huso == 31), coords = c('Coord.X', 'Coord.Y'), crs = 25831), crs = 4326)
  puntos <- rbind(points29, points30, points31)
  return(puntos)
}

presence_points <- to_point(presence)
mapview(presence_points)

absence_points <- to_point(absence)
mapview(absence_points)
# same areas overall, so maybe we'll need to use pseudoabsences

unique(presence$Ano_CS) # years 2016, 2017 and 2018

allpoints <- rbind(presence_points, absence_points)

# random pseudoabsence generation

# spanish provinces borders
provinces <- mapSpain::esp_get_prov()
# excluding Canarias and Baleares
provinces <- provinces[!provinces$iso2.prov.name.es %in% c("Las Palmas", "Santa Cruz de Tenerife", "Baleares"), ]
# dissolve
provinces$campo <- 1
mask_spain <- dissolve(provinces, field='campo')

# environmental variables

# CLIMATE

# processing of present tas
processing <- function(mi_pattern, mi_calc){
  listadefiles <- list.files("spatial_data/climate/present",
                             pattern = mi_pattern,
                             recursive = TRUE,
                             full.names = TRUE)
  print(listadefiles)
  stackdefiles <- stack(listadefiles)
  maskk <- st_transform(mask_spain, crs=crs(stackdefiles))
  cortado <- raster::crop(stackdefiles, extent(maskk))
  masqueado <- raster::mask(cortado, maskk)
  raster_final <- calc(masqueado, mi_calc)
  return(raster_final)
}

mean_temp_2016 <- processing('2016', mi_calc = mean)
mean_temp_2017 <- processing('2017', mi_calc = mean)
mean_temp_2018 <- processing('2018', mi_calc = mean)

# processing of future tas
processing <- function(mi_pattern, mi_calc){
  listadefiles <- list.files("spatial_data/climate/future",
                             pattern = mi_pattern,
                             recursive = TRUE,
                             full.names = TRUE)
  print(listadefiles)
  stackdefiles <- stack(listadefiles)
  maskk <- st_transform(mask_spain, crs=crs(stackdefiles))
  cortado <- raster::crop(stackdefiles, extent(maskk))
  masqueado <- raster::mask(cortado, maskk)
  raster_final <- calc(masqueado, mi_calc)
  return(raster_final)
}

tas_gfdl_1 <- processing("^CHELSA_gfdl.*ssp126_tas_.*\\.tif$", mi_calc = mean) 
tas_ipsl_1 <- processing("^CHELSA_ipsl.*ssp126_tas_.*\\.tif$", mi_calc = mean) 
tas_mpi_1 <- processing("^CHELSA_mpi.*ssp126_tas_.*\\.tif$", mi_calc = mean) 
tas_mri_1 <- processing("^CHELSA_mri.*ssp126_tas_.*\\.tif$", mi_calc = mean) 
tas_uke_1 <- processing("^CHELSA_uke.*ssp126_tas_.*\\.tif$", mi_calc = mean) 

tas_1 <- calc(stack(tas_gfdl_1,tas_ipsl_1,tas_mpi_1,tas_mri_1,tas_uke_1), mean)

tas_gfdl_5 <- processing("^CHELSA_gfdl.*ssp585_tas_.*\\.tif$", mi_calc = mean) 
tas_ipsl_5 <- processing("^CHELSA_ipsl.*ssp585_tas_.*\\.tif$", mi_calc = mean) 
tas_mpi_5 <- processing("^CHELSA_mpi.*ssp585_tas_.*\\.tif$", mi_calc = mean) 
tas_mri_5 <- processing("^CHELSA_mri.*ssp585_tas_.*\\.tif$", mi_calc = mean) 
tas_uke_5 <- processing("^CHELSA_uke.*ssp585_tas_.*\\.tif$", mi_calc = mean) 

tas_5 <- calc(stack(tas_gfdl_5,tas_ipsl_5,tas_mpi_5,tas_mri_5,tas_uke_5), mean)

# z-score formula: (x-mu)/sigma
# we're doing it for 5 rasters: present (2016, 2017, 2018), scenario 126 and scenario 585
manual_zscore <- function(raster1, raster2, raster3, raster4, raster5) {
  raster1v <- na.omit(as.vector(raster1))
  raster2v <- na.omit(as.vector(raster2))
  raster3v <- na.omit(as.vector(raster3))
  raster4v <- na.omit(as.vector(raster4))
  raster5v <- na.omit(as.vector(raster5))
  raster6 <- c(raster1v, raster2v, raster3v, raster4v, raster5v)
  assign(deparse(substitute(raster1)), (raster1-mean(raster6))/sd(raster6), envir = .GlobalEnv)
  assign(deparse(substitute(raster2)), (raster2-mean(raster6))/sd(raster6), envir = .GlobalEnv)
  assign(deparse(substitute(raster3)), (raster3-mean(raster6))/sd(raster6), envir = .GlobalEnv)
  assign(deparse(substitute(raster4)), (raster4-mean(raster6))/sd(raster6), envir = .GlobalEnv)
  assign(deparse(substitute(raster5)), (raster5-mean(raster6))/sd(raster6), envir = .GlobalEnv)
}

manual_zscore(mean_temp_2016, mean_temp_2017, mean_temp_2018, tas_1, tas_5)
plot(tas_5) # done

# resample to 10x10km

mean_temp_2016 <- aggregate(mean_temp_2016, fact = 10, fun='mean')
mean_temp_2017 <- aggregate(mean_temp_2017, fact = 10, fun='mean')
mean_temp_2018 <- aggregate(mean_temp_2018, fact = 10, fun='mean')
tas_1 <- aggregate(tas_1, fact = 10, fun='mean')
tas_5 <- aggregate(tas_5, fact = 10, fun='mean')

# aggregate observations per cells

base <- mean_temp_2016
base[] <- 1:ncell(base)
allpoints$cell_number <- raster::extract(base, allpoints)

allpoints <- allpoints %>% 
  group_by(cell_number) %>%
  summarise(geometry = geometry, Total_CS = sum(Total_CS))

# extract env values to points

env_df <- data.frame(st_coordinates(allpoints), Total_CS = allpoints$Total_CS, mean_temp = NA)

for (i in 1:nrow(allpoints)) {
  if (allpoints[i,'Ano_CS']==2016) {
    env_df[i, 'mean_temp'] <- raster::extract(mean_temp_2016, allpoints[i,])
  }
  else if (allpoints[i,'Ano_CS']==2017) {
    env_df[i, 'mean_temp'] <- raster::extract(mean_temp_2017, allpoints[i,])
  }
  else if (allpoints[i,'Ano_CS']==2018) {
    env_df[i, 'mean_temp'] <- raster::extract(mean_temp_2018, allpoints[i,])
  }
}

