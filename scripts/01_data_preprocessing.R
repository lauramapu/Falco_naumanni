# presence points and environmental variables

rm(list=ls())
source('scripts/utils.R')

# read census data
census <- read_xls("data/primilla_Nacional_2016-2018_V5.xls", sheet = 'CENSO NACIONAL')
# View(census)

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

unique(presence$Ano_CS) # years 2016, 2017 and 2018

# spanish provinces borders
spain <- getpeninsularspain()

# handle observations

presence_points <- raw_presence_points %>%
  mutate(Ano_CS = ifelse(is.na(Ano_CS), 2016, Ano_CS)) %>% # all year NA belong to 2016 (see census documentation)
  mutate(Parejas_CE = ifelse(is.na(Parejas_CE), 0, Parejas_CE)) %>% # all NA in Parejas is 0
  mutate(Total_CS = ifelse(is.na(Total_CS), 0, Total_CS)) %>% # all NA in Total is 0
  mutate(Total = ifelse(2 * Parejas_CE >= Total_CS, 2 * Parejas_CE, Total_CS)) %>% # if Parejas*2 >= Total then parejas*2, else Total
  dplyr::select(Total, Ano_CS, geometry) %>% 
  cbind(st_coordinates(raw_presence_points))

# absence_points <- to_point(absence)
# mapview(absence_points)
# same areas overall so we're not including them

# join obs per raster cells (sum up values per 1x1km cells)
# load base raster
base <- raster('spatial_data/climate/present/CHELSA_tas_01_2016_V.2.1.tif') %>%
  crop(extent(spain)) %>%
  mask(spain)
# assign unique value to each pixel
base[] <- 1:ncell(base)
# raster::extract centroid coordinates
basedf <- base %>%
  as.data.frame(xy=T)
# rename to later join
basedf <- basedf %>%
  rename(ncell = colnames(basedf)[3])
# raster::extract values to points
presence_points$ncell <- raster::extract(base, presence_points)
# group by, summarize, and join with centroids coordinates
modeling_points <- presence_points %>%
  st_drop_geometry() %>%
  group_by(ncell) %>%
  summarize(Total = sum(Total), Ano_CS = modal(Ano_CS)) %>%
  left_join(basedf, by = 'ncell') %>%
  st_as_sf(coords=c('x','y'), crs=4326)
# add coords
modeling_points <- cbind(modeling_points, st_coordinates(modeling_points))
modeling_points <- modeling_points[,-1]

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
  maskk <- st_transform(spain, crs=crs(stackdefiles))
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
  maskk <- st_transform(spain, crs=crs(stackdefiles))
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
  maskk <- st_transform(spain, st_crs(raster_layer))
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
# raster::extract unique values from corine
unique_values <- unique(na.omit(values(raw_landuses)))

# loop through each value and raster::extract binary raster
binary_rasters <- list()
for (value in unique_values) {
  binary_raster <- raw_landuses == value
  binary_rasters[[as.character(value)]] <- binary_raster
}

binary_rasters$"43" 

# now we sum up the following rasters:
# forest: 41, 42, 43
# water: 11, 12, 13
# basing on the known ecological preferences of steppe birds and for simplification purposes

# in further analysis we can see high correlation between all settlements
# so we join them all together

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
linear_infrastr <- raster(x) %>%
  calc(fun = function(x) log10(x + 1)) %>%
  mask_raster()

# WATER STREAMS
# stream length per pixel

# load layers
streams <- st_transform(rbind(st_read("spatial_data/streams/A_RiosCompletosv2.shp"),
                                        st_read("spatial_data/streams/M_RiosCompletosv2.shp")), crs=4326)

# transform to terra format
v <- streams %>% vect() %>% as.lines()
r <- rast(tas_2016)

# calculate length
x <- rasterizeGeom(v, r, "length")    
plot(x)

# log transform
streams <- raster(x) %>%
  calc(fun = function(x) log10(x + 1)) %>%
  mask_raster()

# resource avalability: preys

base <- tas_2016
base[] <- 1:ncell(base) # unique values

files <- list.files('spatial_data/preys', recursive=T, full.names = T, pattern='.shp')

preys <- base
preys[] <- 0

for (file in files) {
  shp <- st_read(file) %>% st_set_crs(25830) %>% st_transform(crs(base))
  values <- na.omit(unlist(raster::extract(base, shp))) # extract pixels in each cell
  binary <- calc(base, function(x) { ifelse(x %in% values, 1, 0) }) # generate binary raster
  preys <- preys + binary # accumulate binaries in the raster
}

preys <- mask_raster(scale(preys))

# resample to same extent and res than climate variables
# bilinear for continuous and ngb for categorical

resampling <- function(x, method) {
  resampled <- projectRaster(x, tas_2016, method = method)
  return(resampled)
}

topo_variables <- lapply(topo_variables, function(x) resampling(x, method = 'bilinear'))
land_variables <- lapply(land_variables, function(x) resampling(x, method = 'ngb'))
linear_infrastr <- resampling(linear_infrastr, method='bilinear')
streams <- resampling(streams, method='bilinear')
preys <- resampling(preys, method='bilinear') # not categorical but integers

# now we need to calculate environmental values in buffers of two different radius for each point and assign those values to each point
# we will generate two different dataset and the bind: env values in a 15km buffer, and values in a 60km buffer
# for all continuous values we will calculate the mean value
# for land uses we will calculate 

env2016 <- c(tas_2016, tasmin_2016, tasmax_2016, tasseas_2016, pr_2016, prseas_2016,
           topo_variables, linear_infrastr, streams, preys, land_variables)
env2017 <- c(tas_2017, tasmin_2017, tasmax_2017, tasseas_2017, pr_2017, prseas_2017,
           topo_variables, linear_infrastr, streams, preys, land_variables)
env2018 <- c(tas_2018, tasmin_2018, tasmax_2018, tasseas_2018, pr_2018, prseas_2018,
           topo_variables, linear_infrastr, streams, preys, land_variables)
envf1 <- c(tas_1, tasmin_1, tasmax_1, tasseas_1, pr_1, prseas_1,
           topo_variables, linear_infrastr, streams, preys, land_variables)
envf5 <- c(tas_5, tasmin_5, tasmax_5, tasseas_5, pr_5, prseas_5,
           topo_variables, linear_infrastr, streams, preys, land_variables)

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
columnnames <- c('tmean', 'tmax', 'tmin', 'tseas', 'pr', 'prseas', 'elev', 'slope', 'linear_infrastr', 'streams', 'preys',
                 'forest', 'water', 'shrub', 'bare', 'for_shr_bare', 'settle_med', 'settle_low', 'for_shr_grass',
                 'for_shr_crop', 'for_shr_agric', 'mosaic_low', 'settle_high', 'crop_med', 'grass_low', 
                 'crop_high', 'grass_med', 'grass_high', 'mosaic_med', 'ext_perm_crop', 'crop_low', 'int_perm_crop', 
                 'mosaic_high')

writeraster <- function(x, name, folder) {
  filename <- paste0('spatial_data/',folder,'/', name, '.asc')
  writeRaster(x, filename, format='ascii', overwrite=TRUE)
}

dir.create('spatial_data/env2016')
dir.create('spatial_data/env2017')
dir.create('spatial_data/env2018')
dir.create('spatial_data/envf1')
dir.create('spatial_data/envf5')

names(env2016) <- columnnames
names(env2017) <- columnnames
names(env2018) <- columnnames
names(envf1) <- columnnames
names(envf5) <- columnnames

mapply(writeraster, env2016, names(env2016), 'env2016')
mapply(writeraster, env2017, names(env2017), 'env2017')
mapply(writeraster, env2018, names(env2018), 'env2018')
mapply(writeraster, envf1, names(envf1), 'envf1')
mapply(writeraster, envf5, names(envf5), 'envf5')

# convert to stack
stack2016 <- stack(env2016)
stack2017 <- stack(env2017)
stack2018 <- stack(env2018)

# no buffer (cell extraction)
# empty df
env_pres <- data.frame(matrix(ncol=length(env2016), nrow=0))

for (i in 1:nrow(modeling_points)) {
  
  # extract values from each stack per year
  if (modeling_points$Ano_CS[i]==2016) {
    values <- as.data.frame(raster::extract(stack2016, modeling_points[i,]))
  }
  else if (modeling_points$Ano_CS[i]==2017) {
    values <- as.data.frame(raster::extract(stack2017, modeling_points[i,]))
  }
  else if (modeling_points$Ano_CS[i]==2018) {
    values <- as.data.frame(raster::extract(stack2018, modeling_points[i,]))
  }
  
  # bind
  env_pres <- rbind(env_pres, values)
}

env_df_0 <- cbind(modeling_points, env_pres)
colnames(env_df_0)[5:37] <- columnnames

# find nearest cell to 1898 (tarifa) and 1884 (conil de la frontera), both 2016
env2016_df <- na.omit(as.data.frame(stack2016, xy=T))
# extract coords
coord_df_0 <- env_df_0[1884, c("X", "Y")]
# calculate abs diffs between x and y
env2016_df$diff_x <- abs(env2016_df$x - coord_df_0$X)
env2016_df$diff_y <- abs(env2016_df$y - coord_df_0$Y)
# sum diffs
env2016_df$total_diff <- env2016_df$diff_x + env2016_df$diff_y
# find the lowest
min_diff_index <- which.min(env2016_df$total_diff)
# raster::extract row
closest_row <- env2016_df[min_diff_index, ]
# print
print(closest_row)
# assign
env_df_0[1884, 5:37] <- closest_row[,3:35]

# extract coords
coord_df_0 <- env_df_0[1898, c("X", "Y")]
# calculate abs diffs between x and y
env2016_df$diff_x <- abs(env2016_df$x - coord_df_0$X)
env2016_df$diff_y <- abs(env2016_df$y - coord_df_0$Y)
# sum diffs
env2016_df$total_diff <- env2016_df$diff_x + env2016_df$diff_y
# find the lowest
min_diff_index <- which.min(env2016_df$total_diff)
# raster::extract row
closest_row <- env2016_df[min_diff_index, ]
# print
print(closest_row)
# assign
env_df_0[1898, 5:37] <- closest_row[,3:35]

# save
write.csv2(st_drop_geometry(env_df_0), "data/fnaumanni_0km.csv", row.names=F)

# 2km buffer (~12km2)
# empty df
env_pres <- data.frame(matrix(ncol=length(env2016), nrow=0))
# here we also calculate shannon diversity index
for (i in 1:nrow(modeling_points)) {
  
  # generate buffer
  buffer <- st_buffer(modeling_points[i,], dist = 2000)
  
  # extract values from each stack per year
  if (modeling_points$Ano_CS[i]==2016) {
    values <- na.omit(as.data.frame(raster::extract(stack2016, buffer)))
  }
  else if (modeling_points$Ano_CS[i]==2017) {
    values <- na.omit(as.data.frame(raster::extract(stack2017, buffer)))
  }
  else if (modeling_points$Ano_CS[i]==2018) {
    values <- na.omit(as.data.frame(raster::extract(stack2018, buffer)))
  }
  
  # calculate mean for continuous and proportion for categorical
  mean_cont <- colMeans(values[,1:11]) # continuous variables
  landuses <- colSums(values[,12:33]) # presence of each land sue
  total <- sum(landuses)
  mean_cat <- landuses / total # proportion of ones vs total
  # calculate shannon diversity index (H') for land uses
  # max H' is log(22) = 3.091042
  shannon <- -sum(mean_cat[mean_cat > 0] * log(mean_cat[mean_cat > 0]))
  # bind
  mean_values <- c(mean_cont, mean_cat, shannon)
  env_pres <- rbind(env_pres, mean_values)
}

env_df_15 <- cbind(modeling_points, env_pres)
colnames(env_df_15)[5:37] <- columnnames
colnames(env_df_15)[38] <- 'shannon_index'
# # scale landuses (many zeros)
# env_df_15[,16:37] <- scale(st_drop_geometry(env_df_15)[,16:37]) # zero variance variables get NA
# env_df_15[is.na(env_df_15)] <- 0
# save
write.csv2(st_drop_geometry(env_df_15), "data/fnaumanni_15km.csv", row.names=F)

# 4.5km buffer (~64km2)
# empty df
env_pres <- data.frame(matrix(ncol=length(env2016), nrow=0))

for (i in 1:nrow(modeling_points)) {
  buffer <- st_buffer(modeling_points[i,], dist = 4500)
  if (modeling_points$Ano_CS[i]==2016) {
    values <- na.omit(as.data.frame(raster::extract(stack2016, buffer)))
  }
  else if (modeling_points$Ano_CS[i]==2017) {
    values <- na.omit(as.data.frame(raster::extract(stack2017, buffer)))
  }
  else if (modeling_points$Ano_CS[i]==2018) {
    values <- na.omit(as.data.frame(raster::extract(stack2018, buffer)))
  }
  mean_cont <- colMeans(values[,1:11]) 
  mean_cat <- colSums(values[,12:33]) / nrow(values) 
  mean_values <- c(mean_cont, mean_cat)
  env_pres <- rbind(env_pres, mean_values)
}

env_df_60 <- cbind(modeling_points, env_pres)
colnames(env_df_60)[5:37] <- columnnames
# # scale buffer results
# env_df_60[,16:37] <- scale(st_drop_geometry(env_df_60)[,16:37])
# env_df_60[is.na(env_df_60)] <- 0
# save
write.csv2(st_drop_geometry(env_df_60), "data/fnaumanni_60km.csv", row.names=F)

# add suffix to col names in each dataset and join
colnames(env_df_0)[5:37] <- paste(colnames(env_df_0)[5:37], "0", sep = "_")
colnames(env_df_15)[5:37] <- paste(colnames(env_df_15)[5:37], "15", sep = "_")
colnames(env_df_60)[5:37] <- paste(colnames(env_df_60)[5:37], "60", sep = "_")
env_df <- cbind(env_df_0, env_df_15[,5:37], env_df_60[,5:37], env_df_15[,'shannon_index'])
env_df <- env_df %>% dplyr::select(-geometry.1, -geometry.2, -geometry.3) %>% st_drop_geometry() # erase geometry

# add log column of response because highly imbalanced
env_df$logTotal <- log10(env_df$Total)

# save complete database
write.csv2(env_df, "data/fnaumanni_envdata.csv", row.names=F)

# plot histograms of all variables
all_cols <- colnames(env_df[5:105])
plot_list <- list()

for (col in all_cols) {
  p <- ggplot(env_df, aes_string(x = col)) +
    geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = col, x = col, y = "Frequency")
  
  plot_list[[col]] <- p
}

all_hists <- grid.arrange(grobs = plot_list, ncol = 11) 
ggsave('data/all_variables_histograms.jpg', all_hists, width=30, height=15)

# modelling observations versus shannon index
# just to preeliminarly check if there's any relation
plot(env_df$logTotal, env_df$shannon_index)
model <- lm(logTotal ~ shannon_index, data = env_df)
summary(model)
# shannon index has no significant relation with abundance (p-value>0.1)
# but still there's some kind of positive relation (more abundance on more heterogeneity)

# we're calculating correlation between zero distance and both buffers to check if there are significant corrs
# first with the smallest buffer
cor_matrix <- cor(env_df[,5:70], use = "complete.obs")
# to df
cor_df <- as.data.frame(as.table(cor_matrix))
# filter pairs with cor higher than 0.80
high_corr <- subset(cor_df, abs(Freq) > 0.75 & abs(Freq) < 1)
# order per the magnitude of the correlation
high_corr <- high_corr[order(-abs(high_corr$Freq)), ]
# print
print(high_corr)
# there's high corr between all climatic variables, elevation, preys, and some landuses, including:
# crop_high, crop_med, crop_low, grass_med (extense landuses in general)

# now we're calculating corrs between zero distance and the largest buffer
cor_matrix <- cor(env_df[,c(5:37, 71:103)], use = "complete.obs")
cor_df <- as.data.frame(as.table(cor_matrix))
high_corr <- subset(cor_df, abs(Freq) > 0.75 & abs(Freq) < 1)
high_corr <- high_corr[order(-abs(high_corr$Freq)), ]
print(high_corr)
# same happens, high corr between all climatic variables, elevation and preys, but no landuses!

# basing on this, we're choosing only zero distance values for climatic variables and elevation
# for crop_high, crop_med, crop_low and grass_med we're using zero distance and largest buffer 
# for the rest of variables, we're using zero distance and selecting buffer size with the following criteria:

# we're choosing 15km2 buffer when H' > mean(range(H')) and 60km2 when lower
# (lower habitat area when diversity is high and viceversa)
medium_shannon <- mean(range(env_df$shannon_index))
mean_shannon <- mean(env_df$shannon_index)
median_shannon <- median(env_df$shannon_index)
# we're choosing median to cut obs in half
# mean and median are very similar

# add these variables to the final dataset
shannon_env <- matrix(ncol=length(env2016), nrow=0)
  
# select buffer size per observation by shannon index
for (i in 1:nrow(env_df)) {
  if (env_df[i, 'shannon_index'] >= median_shannon) {
    shannon_env <- rbind(shannon_env, as.numeric(st_drop_geometry(env_df_15[i,5:37])))
  }
  else if (env_df[i, 'shannon_index'] < median_shannon) {
    shannon_env <- rbind(shannon_env, as.numeric(st_drop_geometry(env_df_60[i,5:37])))
  }
}

#shannon_df <- as.data.frame(shannon_env)
colnames(shannon_env) <- columnnames

env_final <- cbind(env_df[, c('Total', 'logTotal', 'Ano_CS', 'X', 'Y', 'shannon_index')],
                          shannon_env)

write.csv2(env_final, 'data/fnaumanni_byshannon.csv', row.names=F)

# CORRELATION TEST

rm(list = ls())

# load env df
env_df <- read.csv2("data/fnaumanni_byshannon.csv")
# all cols as numeric
env_df[ , sapply(env_df, is.numeric)] <- lapply(env_df[ , sapply(env_df, is.numeric)], as.numeric)

# select only env cols
corrdata <- env_df[,5:39] # same variables but with different buffer size
# # some landuse_0 variables get 0 in every point so can't be computed and thus we erase
# zero_cols <- sapply(corrdata, function(col) all(col == 0))
# corrdata <- corrdata[, !zero_cols]

# Calculate the correlation matrix
cor_matrix <- cor(corrdata, use = "complete.obs")
# save following plot
png("results/correlation_plot.png", width = 3000, height = 2500)
# corrplot
corrplot(cor_matrix, type = "full", method = "number")
# finish burn
dev.off()
# to df
cor_df <- as.data.frame(as.table(cor_matrix))
# filter pairs with cor higher than 0.80
high_corr <- subset(cor_df, abs(Freq) > 0.75 & abs(Freq) < 1)
# order per the magnitude of the correlation
high_corr <- high_corr[order(-abs(high_corr$Freq)), ]
# print
print(high_corr)
write.csv2(cor_matrix, "results/cor_matrix.csv", row.names=T)

# we're erasing: tmin, tmax and elev
# and also mosaic_high and grass_high because zero std deviation

env_df <- env_df %>% select(-tmin, -tmax, -elev, -grass_high, -mosaic_high)
write.csv2(env_df, "data/modeling_data.csv", row.names=F)
