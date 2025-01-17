
# density kernels of variables variability in the modeling points vs all the territory
# to check if there's actually enough variability in our points to model anything

rm(list=ls())
source('scripts/utils.R')

dataset <- read.csv2('data/modeling_data.csv')
envpres <- read.csv2('data/env2018_buffered.csv')

variables <- colnames(dataset)[6:39]
#variables <- setdiff(variables, "preys")

plot_list <- lapply(variables, function(variable) {
  ggplot() +
    geom_density(data = dataset, aes_string(x = variable), color = "blue", fill = "blue", alpha = 0.3) +  
    geom_density(data = envpres, aes_string(x = variable), color = "red", fill = "red", alpha = 0.3) +  
    labs(x = variable, y = "Density", title = variable) +
    theme_minimal()
})

variability <- grid.arrange(grobs = plot_list, ncol = 7)
ggsave('data/variables_variability.jpg', variability, width=25, height=12)




