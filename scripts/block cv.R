# spatial block cross validation

# create grid and folds

# spanish provinces borders
provinces <- mapSpain::esp_get_prov()
# excluding Canarias and Baleares
provinces <- provinces[!provinces$iso2.prov.name.es %in% c("Las Palmas", "Santa Cruz de Tenerife",
                                                           "Baleares", "Melilla", "Ceuta"), ]
# dissolve
provinces$campo <- 1
mask_spain <- dissolve(provinces, field='campo')

# 50km grid
grid_50 <- st_make_grid(
  mask_spain,
  offset = st_bbox(mask_spain)[c("xmin", "ymin")],
  cellsize = c(0.4166666, 0.4166666), # 50km in degrees
  crs = 4326, 
  what = "polygons",
  square = TRUE,
)
plot(grid_50)
mask_spain <- st_transform(mask_spain, crs=st_crs(grid_50))
plot(mask_spain, add=T)

# intersect to erase out of border cells
grid_mask <- st_intersection(grid_50, mask_spain)
plot(grid_mask)

# grid to sf
grid_df <- st_sf(id = seq_along(grid_mask), geometry = grid_mask)
# add column to group into 10 spatial blocks
grid_df$fold <- cut(seq(nrow(grid_df)), breaks = 5, labels = FALSE)
# visualize
mapview(grid_df, zcol = 'fold')

# erase polygons 4 and 157 (two super little rocks in the sea)
grid_df <- grid_df %>% filter(!(id %in% c(4, 157)))
plot(grid_df)

# now we need to spatially join this grid with the points to assign cell ids 
intersect <- st_join(dataset_sf, grid_df)
blocks <- sort(unique(intersect$id)) # a total of 319 blocks
table(intersect$fold) # number of points between blocks differs much sometimes

# Filtrar los puntos con NA en fold
na_points <- intersect %>% filter(is.na(fold))

# Filtrar los puntos con valores no NA en fold
non_na_points <- intersect %>% filter(!is.na(fold))

# Encontrar el índice del punto más cercano con valor no NA
nearest_idx <- st_nearest_feature(na_points, non_na_points)

# Extraer los valores de fold correspondientes a estos índices
nearest_folds <- non_na_points$fold[nearest_idx]

# Asignar los valores de fold a los puntos NA
na_points$fold <- nearest_folds

# Combinar los puntos con NA corregidos con los puntos originales sin NA
intersect <- bind_rows(non_na_points, na_points)

# Verificar que ya no hay NAs en fold
sum(is.na(intersect$fold))