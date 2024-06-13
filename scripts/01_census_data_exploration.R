# presence points and environmental variables

library(readxl)
library(sf)
library(mapview)
library(mapSpain)
library(MazamaSpatialUtils)
library(raster)
library(dplyr)
library(randomForest)
library(units)
library(terra)
library(corrplot)
library(car) # vif

# read census data
census <- read_xls("data/primilla_Nacional_2016-2018_V5.xls", sheet = 'CENSO NACIONAL')
View(census)

# filter presences
presence <- census %>%
  subset(census$Total_CS > 0 | census$Parejas_CE > 0) %>% # presence either in simple or exhaustive census
  filter(Provincia != 'Melilla')
  
# absence <- subset(census, census$Total_CS == 0) # true absences?

# reproyect to 4326 from corresponding huso (29, 30 and 31)
to_point <- function(x) {
  points29 <- st_transform(st_as_sf(subset(x, x$Huso == 29), coords = c('Coord.X', 'Coord.Y'), crs = 25829), crs = 4326)
  points30 <- st_transform(st_as_sf(subset(x, x$Huso == 30), coords = c('Coord.X', 'Coord.Y'), crs = 25830), crs = 4326)
  points31 <- st_transform(st_as_sf(subset(x, x$Huso == 31), coords = c('Coord.X', 'Coord.Y'), crs = 25831), crs = 4326)
  puntos <- rbind(points29, points30, points31)
  return(puntos)
}

raw_presence_points <- to_point(presence)
mapview(raw_presence_points, zcol='Total_CS')

# random pseudoabsence generation
# spanish provinces borders
provinces <- mapSpain::esp_get_prov()
# excluding Canarias and Baleares
provinces <- provinces[!provinces$iso2.prov.name.es %in% c("Las Palmas", "Santa Cruz de Tenerife", "Baleares"), ]
# dissolve
provinces$campo <- 1
mask_spain <- dissolve(provinces, field='campo')

set.seed(1)
pseudoabsences <- as.data.frame(st_sample(mask_spain, size = 1000, type = "random"))

# convert to single spatial objects to later extract raster values
PA_to_points <- function(original_points) {
  coordinates_obj <- as.data.frame(st_coordinates(original_points$geometry)) 
  points_object <- st_as_sf(coordinates_obj, coords = c("X", "Y"), crs = 4326)
  return(points_object)
}
pseudoabsences_sf <- PA_to_points(pseudoabsences)

unique(presence$Ano_CS) # years 2016, 2017 and 2018

# allpoints <- rbind(presence_points, absence_points) # true absences

# handle observations

presence_points <- raw_presence_points %>%
  mutate(Ano_CS = ifelse(is.na(Ano_CS), 2016, Ano_CS)) %>% # all year NA belong to 2016 (see census documentation)
  mutate(Parejas_CE = ifelse(is.na(Parejas_CE), 0, Parejas_CE)) %>% # all NA in Parejas is 0
  mutate(Total_CS = ifelse(is.na(Total_CS), 0, Total_CS)) %>% # all NA in Total is 0
  mutate(Total = ifelse(2 * Parejas_CE >= Total_CS, 2 * Parejas_CE, Total_CS)) %>% # if Parejas*2 >= Total then parejas, else Total
  select(Total, Ano_CS, geometry) %>% 
  cbind(st_coordinates(raw_presence_points))

# same format as presence
set.seed(21)
pseudoabsences_sf <- pseudoabsences_sf %>%
  mutate(Total = 0, Ano_CS = sample(presence_points$Ano_CS, n(), replace = TRUE)) %>%
  cbind(st_coordinates(pseudoabsences_sf))

allpoints <- rbind(presence_points, pseudoabsences_sf) 

# absence_points <- to_point(absence)
# mapview(absence_points)
# same areas overall

# distances between colonies (presence points)
# we need a function to first calculate nearest point and then extract distance

(nearest = st_nearest_feature(presence_points))
(dist = st_distance(presence_points, presence_points[nearest,], by_element=TRUE))
presence_points$min_dist <- drop_units(dist) # meters

zero_distance <- presence_points[presence_points$min_dist<1,]
mapview(zero_distance)

# environmental variables

##### CLIMATE #####

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

tas_2016 <- processing(".*tas_.*2016.*", mi_calc = mean)
tas_2017 <- processing(".*tas_.*2017.*", mi_calc = mean)
tas_2018 <- processing(".*tas_.*2018.*", mi_calc = mean)

tasmax_2016 <- processing(".*tasmax_.*2016.*", mi_calc = mean)
tasmax_2017 <- processing(".*tasmax_.*2017.*", mi_calc = mean)
tasmax_2018 <- processing(".*tasmax_.*2018.*", mi_calc = mean)

tasmin_2016 <- processing(".*tasmin_.*2016.*", mi_calc = mean)
tasmin_2017 <- processing(".*tasmin_.*2017.*", mi_calc = mean)
tasmin_2018 <- processing(".*tasmin_.*2018.*", mi_calc = mean)

# temp seasonality = stddev * 100
tasseas_2016 <- processing(".*tas_.*2016.*", mi_calc = sd) * 100
tasseas_2017 <- processing(".*tas_.*2017.*", mi_calc = sd) * 100
tasseas_2018 <- processing(".*tas_.*2018.*", mi_calc = sd) * 100

# present pr
pr_2016 <- processing(".*pr.*2016.*", mi_calc = sum)
pr_2017 <- processing(".*pr.*2017.*", mi_calc = sum)
pr_2018 <- processing(".*pr.*2018.*", mi_calc = sum)

# prec seasonality (cv = (stdv / (1 + mean)) * 100)
prseas_2016 <- processing(".*pr.*2016.*", mi_calc = sd) / (1 + processing(".*pr.*2016.*", mi_calc = mean)) * 100
prseas_2017 <- processing(".*pr.*2017.*", mi_calc = sd) / (1 + processing(".*pr.*2017.*", mi_calc = mean)) * 100
prseas_2018 <- processing(".*pr.*2018.*", mi_calc = sd) / (1 + processing(".*pr.*2018.*", mi_calc = mean)) * 100

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

tassd_gfdl_1 <- processing("^CHELSA_gfdl.*ssp126_tas_.*\\.tif$", mi_calc = sd) 
tassd_ipsl_1 <- processing("^CHELSA_ipsl.*ssp126_tas_.*\\.tif$", mi_calc = sd) 
tassd_mpi_1 <- processing("^CHELSA_mpi.*ssp126_tas_.*\\.tif$", mi_calc = sd) 
tassd_mri_1 <- processing("^CHELSA_mri.*ssp126_tas_.*\\.tif$", mi_calc = sd) 
tassd_uke_1 <- processing("^CHELSA_uke.*ssp126_tas_.*\\.tif$", mi_calc = sd) 

tasseas_1 <- calc(stack(tassd_gfdl_1,tassd_ipsl_1,tassd_mpi_1,tassd_mri_1,tassd_uke_1), mean) *100

tassd_gfdl_5 <- processing("^CHELSA_gfdl.*ssp585_tas_.*\\.tif$", mi_calc = sd) 
tassd_ipsl_5 <- processing("^CHELSA_ipsl.*ssp585_tas_.*\\.tif$", mi_calc = sd) 
tassd_mpi_5 <- processing("^CHELSA_mpi.*ssp585_tas_.*\\.tif$", mi_calc = sd) 
tassd_mri_5 <- processing("^CHELSA_mri.*ssp585_tas_.*\\.tif$", mi_calc = sd) 
tassd_uke_5 <- processing("^CHELSA_uke.*ssp585_tas_.*\\.tif$", mi_calc = sd) 

tasseas_5 <- calc(stack(tassd_gfdl_5,tassd_ipsl_5,tassd_mpi_5,tassd_mri_5,tassd_uke_5), mean) *100

tasmax_gfdl_1 <- processing("^CHELSA_gfdl.*ssp126_tasmax_.*\\.tif$", mi_calc = mean) 
tasmax_ipsl_1 <- processing("^CHELSA_ipsl.*ssp126_tasmax_.*\\.tif$", mi_calc = mean) 
tasmax_mpi_1 <- processing("^CHELSA_mpi.*ssp126_tasmax_.*\\.tif$", mi_calc = mean) 
tasmax_mri_1 <- processing("^CHELSA_mri.*ssp126_tasmax_.*\\.tif$", mi_calc = mean) 
tasmax_uke_1 <- processing("^CHELSA_uke.*ssp126_tasmax_.*\\.tif$", mi_calc = mean) 

tasmax_1 <- calc(stack(tasmax_gfdl_1,tasmax_ipsl_1,tasmax_mpi_1,tasmax_mri_1,tasmax_uke_1), mean)

tasmax_gfdl_5 <- processing("^CHELSA_gfdl.*ssp585_tasmax_.*\\.tif$", mi_calc = mean) 
tasmax_ipsl_5 <- processing("^CHELSA_ipsl.*ssp585_tasmax_.*\\.tif$", mi_calc = mean) 
tasmax_mpi_5 <- processing("^CHELSA_mpi.*ssp585_tasmax_.*\\.tif$", mi_calc = mean) 
tasmax_mri_5 <- processing("^CHELSA_mri.*ssp585_tasmax_.*\\.tif$", mi_calc = mean) 
tasmax_uke_5 <- processing("^CHELSA_uke.*ssp585_tasmax_.*\\.tif$", mi_calc = mean) 

tasmax_5 <- calc(stack(tasmax_gfdl_5,tasmax_ipsl_5,tasmax_mpi_5,tasmax_mri_5,tasmax_uke_5), mean)

tasmin_gfdl_1 <- processing("^CHELSA_gfdl.*ssp126_tasmin_.*\\.tif$", mi_calc = mean) 
tasmin_ipsl_1 <- processing("^CHELSA_ipsl.*ssp126_tasmin_.*\\.tif$", mi_calc = mean) 
tasmin_mpi_1 <- processing("^CHELSA_mpi.*ssp126_tasmin_.*\\.tif$", mi_calc = mean) 
tasmin_mri_1 <- processing("^CHELSA_mri.*ssp126_tasmin_.*\\.tif$", mi_calc = mean) 
tasmin_uke_1 <- processing("^CHELSA_uke.*ssp126_tasmin_.*\\.tif$", mi_calc = mean) 

tasmin_1 <- calc(stack(tasmin_gfdl_1,tasmin_ipsl_1,tasmin_mpi_1,tasmin_mri_1,tasmin_uke_1), mean)

tasmin_gfdl_5 <- processing("^CHELSA_gfdl.*ssp585_tasmin_.*\\.tif$", mi_calc = mean) 
tasmin_ipsl_5 <- processing("^CHELSA_ipsl.*ssp585_tasmin_.*\\.tif$", mi_calc = mean) 
tasmin_mpi_5 <- processing("^CHELSA_mpi.*ssp585_tasmin_.*\\.tif$", mi_calc = mean) 
tasmin_mri_5 <- processing("^CHELSA_mri.*ssp585_tasmin_.*\\.tif$", mi_calc = mean) 
tasmin_uke_5 <- processing("^CHELSA_uke.*ssp585_tasmin_.*\\.tif$", mi_calc = mean) 

tasmin_5 <- calc(stack(tasmin_gfdl_5,tasmin_ipsl_5,tasmin_mpi_5,tasmin_mri_5,tasmin_uke_5), mean)

# future pr

pr_gfdl_1 <- processing("^CHELSA_gfdl.*ssp126_pr_.*\\.tif$", mi_calc = sum) 
pr_ipsl_1 <- processing("^CHELSA_ipsl.*ssp126_pr_.*\\.tif$", mi_calc = sum) 
pr_mpi_1 <- processing("^CHELSA_mpi.*ssp126_pr_.*\\.tif$", mi_calc = sum) 
pr_mri_1 <- processing("^CHELSA_mri.*ssp126_pr_.*\\.tif$", mi_calc = sum) 
pr_uke_1 <- processing("^CHELSA_uke.*ssp126_pr_.*\\.tif$", mi_calc = sum) 

pr_1 <- calc(stack(pr_gfdl_1,pr_ipsl_1,pr_mpi_1,pr_mri_1,pr_uke_1), mean)

pr_gfdl_5 <- processing("^CHELSA_gfdl.*ssp585_pr_.*\\.tif$", mi_calc = sum) 
pr_ipsl_5 <- processing("^CHELSA_ipsl.*ssp585_pr_.*\\.tif$", mi_calc = sum) 
pr_mpi_5 <- processing("^CHELSA_mpi.*ssp585_pr_.*\\.tif$", mi_calc = sum) 
pr_mri_5 <- processing("^CHELSA_mri.*ssp585_pr_.*\\.tif$", mi_calc = sum) 
pr_uke_5 <- processing("^CHELSA_uke.*ssp585_pr_.*\\.tif$", mi_calc = sum) 

pr_5 <- calc(stack(pr_gfdl_5,pr_ipsl_5,pr_mpi_5,pr_mri_5,pr_uke_5), mean)

prsd_gfdl_1 <- processing("^CHELSA_gfdl.*ssp126_pr_.*\\.tif$", mi_calc = sd) 
prsd_ipsl_1 <- processing("^CHELSA_ipsl.*ssp126_pr_.*\\.tif$", mi_calc = sd) 
prsd_mpi_1 <- processing("^CHELSA_mpi.*ssp126_pr_.*\\.tif$", mi_calc = sd) 
prsd_mri_1 <- processing("^CHELSA_mri.*ssp126_pr_.*\\.tif$", mi_calc = sd) 
prsd_uke_1 <- processing("^CHELSA_uke.*ssp126_pr_.*\\.tif$", mi_calc = sd)

prsd_1 <- calc(stack(prsd_gfdl_1,prsd_ipsl_1,prsd_mpi_1,prsd_mri_1,prsd_uke_1), mean)

prmean_gfdl_1 <- processing("^CHELSA_gfdl.*ssp126_pr_.*\\.tif$", mi_calc = mean) 
prmean_ipsl_1 <- processing("^CHELSA_ipsl.*ssp126_pr_.*\\.tif$", mi_calc = mean) 
prmean_mpi_1 <- processing("^CHELSA_mpi.*ssp126_pr_.*\\.tif$", mi_calc = mean) 
prmean_mri_1 <- processing("^CHELSA_mri.*ssp126_pr_.*\\.tif$", mi_calc = mean) 
prmean_uke_1 <- processing("^CHELSA_uke.*ssp126_pr_.*\\.tif$", mi_calc = mean)

prmean_1 <- calc(stack(prmean_gfdl_1,prmean_ipsl_1,prmean_mpi_1,prmean_mri_1,prmean_uke_1), mean)

prseas_1 <- (prsd_1/(1-prmean_1))*100

prsd_gfdl_5 <- processing("^CHELSA_gfdl.*ssp585_pr_.*\\.tif$", mi_calc = sd) 
prsd_ipsl_5 <- processing("^CHELSA_ipsl.*ssp585_pr_.*\\.tif$", mi_calc = sd) 
prsd_mpi_5 <- processing("^CHELSA_mpi.*ssp585_pr_.*\\.tif$", mi_calc = sd) 
prsd_mri_5 <- processing("^CHELSA_mri.*ssp585_pr_.*\\.tif$", mi_calc = sd) 
prsd_uke_5 <- processing("^CHELSA_uke.*ssp585_pr_.*\\.tif$", mi_calc = sd)

prsd_5 <- (prsd_gfdl_5+prsd_ipsl_5+prsd_mpi_5+prsd_mri_5+prsd_uke_5)/5

prmean_gfdl_5 <- processing("^CHELSA_gfdl.*ssp585_pr_.*\\.tif$", mi_calc = mean) 
prmean_ipsl_5 <- processing("^CHELSA_ipsl.*ssp585_pr_.*\\.tif$", mi_calc = mean) 
prmean_mpi_5 <- processing("^CHELSA_mpi.*ssp585_pr_.*\\.tif$", mi_calc = mean) 
prmean_mri_5 <- processing("^CHELSA_mri.*ssp585_pr_.*\\.tif$", mi_calc = mean) 
prmean_uke_5 <- processing("^CHELSA_uke.*ssp585_pr_.*\\.tif$", mi_calc = mean)

prmean_5 <- (prmean_gfdl_5+prmean_ipsl_5+prmean_mpi_5+prmean_mri_5+prmean_uke_5)/5

prseas_5 <- (prsd_5/(1-prmean_5))*100

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

manual_zscore(tas_2016, tas_2017, tas_2018, tas_1, tas_5)
manual_zscore(tasmax_2016, tasmax_2017, tasmax_2018, tasmax_1, tasmax_5)
manual_zscore(tasmin_2016, tasmin_2017, tasmin_2018, tasmin_1, tasmin_5)
manual_zscore(tasseas_2016, tasseas_2017, tasseas_2018, tasseas_1, tasseas_5)
manual_zscore(pr_2016, pr_2017, pr_2018, pr_1, pr_5)
manual_zscore(prseas_2016, prseas_2017, prseas_2018, prseas_1, prseas_5)

##### ELEVATION AND SLOPE ##### 

# crop+mask function for all variables
mask_raster <- function(raster_layer) {
  maskk <- st_transform(mask_spain, st_crs(raster_layer))
  cropped <- raster::crop(raster_layer, extent(maskk))
  masked <- raster::mask(cropped, maskk)
  return(masked)
}

# load rasters
dem_list <- list.files("spatial_data/dem", pattern = "", full.names = TRUE)
# rasters have different extent (dont know why) so we have to merge instead of stack
dem_list <- lapply(dem_list, raster) # read as list of rasters (list.files only reads filepaths)
dem <- do.call(merge, dem_list) # do call to merge with dem_list

# calculate slope
slope <- raster::terrain(dem, opt = "slope", unit = "degrees", neighbors = 8)

# list topography variables and crop/mask by reprojected mask
topo_variables <- list(dem,slope)
topo_variables <- lapply(topo_variables, mask_raster)
topo_variables <- lapply(topo_variables, scale)

### LAND USES

raw_landuses <- raster("spatial_data/land_uses/EU_landSystem.tif")
# extract unique values from corine
unique_values <- unique(na.omit(values(raw_landuses)))

# loop through each value and extract binary raster
binary_rasters <- list()
for (value in unique_values) {
  binary_raster <- raw_landuses == value
  binary_rasters[[as.character(value)]] <- binary_raster
}

binary_rasters$"43" 

# now we sum up the following rasters:
# settlements: 21, 22, 23
# forest: 41, 42, 43
# water: 11, 12, 13
# basing on the known ecological preferences of steppe birds and for simplification purposes

elements_to_erase <- c('0','11', '12', '13', '41', '42', '43')
binary_rasters_filtered <- binary_rasters[setdiff(names(binary_rasters), elements_to_erase)]
filtered_rasters <- stack(binary_rasters_filtered)

forest <- calc(stack(binary_rasters$'41', binary_rasters$'42', binary_rasters$'43'), sum)
water <- calc(stack(binary_rasters$'11', binary_rasters$'12', binary_rasters$'13'), sum)

landuses <- c(forest, water, binary_rasters_filtered)

# crop/mask landuse stack by reprojected mask

land_variables <- lapply(landuses, mask_raster)
# if the scipen error occurs do options(scipen = 99)

# ANTHROPIC PERTURBATION
# road and train length per pixel

# load layers
linear_infrastr <- st_transform(st_join(st_read("spatial_data/roads/RT_VIARIA_CARRETERA/rt_tramo_vial.shp"),
                                        st_read("spatial_data/roads/RT_FFCC/rt_tramofc_linea.shp")), crs=4326)

# transform to terra format
v <- linear_infrastr %>% vect() %>% as.lines()
r <- rast(tas_2016)
 
# calculate length
x <- rasterizeGeom(v, r, "length")    
plot(x)

# transform to raster package format
x <- raster(x)
# log transform
linear_infrastr <- calc(x, fun = function(x) log10(x + 1)) 
linear_infrastr <- mask_raster(linear_infrastr)

# resample to same extent and res than climate variables
# bilinear for continuous and ngb for categorical

resampling <- function(x, method) {
  resampled <- projectRaster(x, tas_1, method = method)
  return(resampled)
}

topo_variables <- lapply(topo_variables, function(x) resampling(x, method = 'bilinear'))
land_variables <- lapply(land_variables, function(x) resampling(x, method = 'ngb'))
linear_infrastr <- resampling(linear_infrastr, method='bilinear')

# now we need to calculate environmental values in buffers of two different radius for each point and assign those values to each point
# we will generate two different dataset and the bind: env values in a 15km buffer, and values in a 60km buffer
# for all continuous values we will calculate the mean value
# for land uses we will calculate 

env2016 <- c(tas_2016, tasmin_2016, tasmax_2016, tasseas_2016, pr_2016, prseas_2016,
           topo_variables, linear_infrastr, land_variables)
env2017 <- c(tas_2017, tasmin_2017, tasmax_2017, tasseas_2017, pr_2017, prseas_2017,
           topo_variables, linear_infrastr, land_variables)
env2018 <- c(tas_2018, tasmin_2018, tasmax_2018, tasseas_2018, pr_2018, prseas_2018,
           topo_variables, linear_infrastr, land_variables)
envf1 <- c(tas_1, tasmin_1, tasmax_1, tasseas_1, pr_1, prseas_1,
           topo_variables, liear_infrastr, land_variables)
envf5 <- c(tas_5, tasmin_5, tasmax_5, tasseas_5, pr_5, prseas_5,
           topo_variables, linear_infrastr, land_variables)

# save and load these objects as rds
saveRDS(env2016, "objects/env2016.rds")
saveRDS(env2017, "objects/env2017.rds")
saveRDS(env2018, "objects/env2018.rds")
saveRDS(envf1, "objects/envf1.rds")
saveRDS(envf5, "objects/envf5.rds")
# env2016 <- readRDS("objects/env2016.rds")
# env2017 <- readRDS("objects/env2017.rds")
# env2018 <- readRDS("objects/env2018.rds")
# envf1 <- readRDS("objects/envf1.rds")
# envf5 <- readRDS("objects/envf5.rds")

# also save a set of rasters to later explore spatial autocorrelation
columnnames <- c('tmean', 'tmax', 'tmin', 'tseas', 'pr', 'prseas', 'elev', 'slope', 'linear_infrastr',
                 'forest', 'water', 'shrub', 'bare', 'for_shr_bare', 'settle_med', 'settle_low', 'for_shr_grass',
                 'for_shr_crop', 'for_shr_agric', 'mosaic_low', 'settle_high', 'crop_med', 'grass_low', 
                 'crop_high', 'grass_med', 'grass_high', 'mosaic_med', 'ext_perm_crop', 'crop_low', 'int_perm_crop', 
                 'mosaic_high')
writeraster <- function(x, name) {
  filename <- paste0('spatial_data/env_asciis/', name, '.asc')
  writeRaster(x, filename, format='ascii', overwrite=TRUE)
}
names(env2016) <- columnnames
mapply(writeraster, env2016, names(env2016))

# convert to stack
stack2016 <- stack(env2016)
stack2017 <- stack(env2017)
stack2018 <- stack(env2018)

env_pres <- data.frame(matrix(ncol=length(env2016), nrow=0))

for (i in 1:nrow(allpoints)) {
  
  # generate buffer
  buffer <- st_buffer(allpoints[i,], dist = 15000)
  
  # extract values from each stack per year
  if (allpoints$Ano_CS[i]==2016) {
    values <- na.omit(as.data.frame(extract(stack2016, buffer)))
  }
  else if (allpoints$Ano_CS[i]==2017) {
    values <- na.omit(as.data.frame(extract(stack2017, buffer)))
  }
  else if (allpoints$Ano_CS[i]==2018) {
    values <- na.omit(as.data.frame(extract(stack2018, buffer)))
  }
  
  # calculate mean for continuous and proportion for categorical
  mean_cont <- colMeans(values[,1:9]) # continuous variables
  mean_cat <- colSums(values[,10:31]) / nrow(values) # land uses (proportion of ones vs total)
  mean_values <- c(mean_cont, mean_cat)
  env_pres <- rbind(env_pres, mean_values)
}

env_df_15 <- cbind(allpoints, env_pres)
colnames(env_df_15)[5:35] <- columnnames
write.csv2(env_df_15, "data/fnaumanni_15km.csv", row.names=F)

# repeat the latter process with a 60km buffer

env_pres <- data.frame(matrix(ncol=length(env2016), nrow=0))

for (i in 1:nrow(allpoints)) {
  buffer <- st_buffer(allpoints[i,], dist = 60000)
  if (allpoints$Ano_CS[i]==2016) {
    values <- na.omit(as.data.frame(extract(stack2016, buffer)))
  }
  else if (allpoints$Ano_CS[i]==2017) {
    values <- na.omit(as.data.frame(extract(stack2017, buffer)))
  }
  else if (allpoints$Ano_CS[i]==2018) {
    values <- na.omit(as.data.frame(extract(stack2018, buffer)))
  }
  mean_cont <- colMeans(values[,1:9]) 
  mean_cat <- colSums(values[,10:31]) / nrow(values) 
  mean_values <- c(mean_cont, mean_cat)
  env_pres <- rbind(env_pres, mean_values)
}

env_df_60 <- cbind(allpoints, env_pres)
colnames(env_df_60)[5:35] <- columnnames
write.csv2(env_df_60, "data/fnaumanni_60km.csv", row.names=F)

# add suffix to col names in each dataset and join
colnames(env_df_15)[5:35] <- paste(colnames(env_df_15)[5:35], "15", sep = "_")
colnames(env_df_60)[5:35] <- paste(colnames(env_df_60)[5:35], "60", sep = "_")
env_df <- cbind(env_df_15, env_df_60[,5:35])
env_df <- env_df %>% select(-geometry.1) # duplicated geometry

# save datasets
write.csv2(env_df, "data/fnaumanni_envdata.csv", row.names=F)

# CORRELATION TEST

rm(list = ls())

# load env df
env_df <- read.csv2("data/fnaumanni_envdata.csv")

# select only env cols
corrdata <- env_df[,5:66] # same variables but with different buffer size
st_geometry(corrdata) <- NULL

# Calculate the correlation matrix
cor_matrix <- cor(corrdata, use = "complete.obs")
# save following plot
png("results/correlation_plot.png", width = 3000, height = 2500)
# corrplot
corrplot(cor_matrix, type = "full", method = "number")
# finish burn
dev.off()
# Convertir la matriz de correlación en un data frame
cor_df <- as.data.frame(as.table(cor_matrix))
# filter pairs with cor higher than 0.80
high_corr <- subset(cor_df, abs(Freq) > 0.75 & abs(Freq) < 1)
# order per the magnitude of the correlation
high_corr <- high_corr[order(-abs(high_corr$Freq)), ]
# print
print(high_corr)

# high correlations in:
vifdata <- corrdata[, -which(names(corrdata) %in% c('tmin_15', 'tmin_60', 'tmax_15', 'tmax_60', 'tseas_60', 'pr_60', 'tmean_60',
                                                    'settle_high_60', 'slope_60', 'crop_med_60', 'int_perm_crop_60', 'settle_high_15',
                                                    'elev_15', 'for_shr_grass_60', 'grass_low_60', 'mosaic_high_15', 'forest_60', 
                                                    'linear_infrastr_60', 'bare_60', 'for_shr_crop_60'))]
# aqui nos sale mucha correlacion entre settle_med y settle_high entonces igual estaria bien juntarlas para no perderlas

### VIF ANALYSIS
model <- lm(tmean_15 ~ ., data = vifdata)
vifmodel <- as.data.frame(car::vif(model))
vifmodel # crop_med_15 = 348
vifmodel

vifdata <- vifdata[, -which(names(vifdata) %in% c('crop_med_15'))]
model <- lm(tmean_15 ~ ., data = vifdata)
vifmodel <- as.data.frame(car::vif(model))
vifmodel # prseas_60 = 10.4

vifdata <- vifdata[, -which(names(vifdata) %in% c('prseas_60'))]
model <- lm(tmean_15 ~ ., data = vifdata)
vifmodel <- as.data.frame(car::vif(model))
vifmodel # all under 10 (acceptable; Ringim et al., 2012)

write.csv2(cor_matrix, "results/cor_matrix.csv", row.names=T)
write.csv2(vifmodel, "results/vif.csv", row.names=T)

# so the data we're using for modeling is the following
final_envdf <- env_df[, -which(names(env_df) %in% c('tmin_15', 'tmin_60', 'tmax_15', 'tmax_60', 'tseas_60', 'pr_60', 'tmean_60',
                                                    'settle_high_60', 'slope_60', 'crop_med_60', 'int_perm_crop_60', 'settle_high_15',
                                                    'elev_15', 'for_shr_grass_60', 'grass_low_60', 'mosaic_high_15', 'forest_60', 
                                                    'linear_infrastr_60', 'bare_60', 'for_shr_crop_60', 'crop_med_15', 'prseas_60'))]
write.csv2(final_envdf, "data/modeling_data.csv", row.names=F)
