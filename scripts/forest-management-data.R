# data
# Art 12 & 17 spp - distribution models current & potential (spatial)
# Select forest prefered or suitable (from EEA habitat preferences for each spp)
# Spp to EEA threats 
# Threats to forest management intensity (see bioclima work)
# Forest management data (spatial)

# gist:
# For each species get the relative value of a pixel's forest for their habitat area 
# 1 = all undisturbed forest
# 0 = no forest
# forest area penalized by threats associated with a given spp and a given management type
# sum across all spp and normalize to 0-1


# Load required libraries
library(terra)        # For spatial data manipulation
library(raster)       # For handling raster data
library(tidyverse)    # For data manipulation and visualization
library(fasterize)    # For fast rasterization of spatial data
library(exactextractr) # For fast exact extraction of raster values using polygons
library(sf)           # For handling simple features (spatial vector data)
rm(list = ls())       # Clear all objects from the workspace

## Map species to each layer

# Load a species distribution model (SDM) as a template raster
spp_template <- rast("data/CurrentSDMs/Amphibians/EnsembleNormThreshold__Alytes_cisternasii.tif")

# Load European Union (EU) forest raster data
EU_forests <- rast("data/Nature_map_EU/Nature_EU_all1.tif")

# Reproject the species template to match the CRS (coordinate reference system) of the EU forests raster
spp_template <- project(spp_template, crs(EU_forests))

# Get the number of cells in the species template raster
ncells <- ncell(spp_template)

# Create a sequence of cell IDs from 1 to the number of cells and assign it to the raster values
values <- 1:ncells
values(spp_template) <- values
names(spp_template) <- "ID"

# Convert the raster to polygons without dissolving the boundaries between them
polygons <- as.polygons(spp_template, dissolve = FALSE)

# Convert the polygons to a simple feature (sf) object
spp_template_shp <- st_as_sf(polygons) 

# Check if a classified forest raster already exists; if not, create it
if (!file.exists("data/Nature_map_EU/Nature_EU_all_rc.tif")) {
  EU_forests <- rast("data/Nature_map_EU/Nature_EU_all1.tif")
  
  # Create a classification matrix
  m <- c(10, 12, 1,
         19, 21, 2,
         30, 33, 3)
  rclmat <- matrix(m, ncol = 3, byrow = TRUE)
  
  # Classify the forest raster based on the matrix
  EU_forests <- classify(EU_forests, rclmat, others = NA)
  
  # Save the classified raster
  writeRaster(EU_forests, "data/Nature_map_EU/Nature_EU_all_rc.tif")
  
  # Further classify the raster into natural, multi-use, and production forests
  forest_natural <- classify(EU_forests, 
                             rcl = matrix(c(0.5, 1.5, 1,
                                            1.5, 5, NA), 
                                          ncol = 3, byrow = TRUE), 
                             others = NA)
  forest_multi <- classify(EU_forests, 
                           rcl = matrix(c(1.5, 2.5, 1,
                                          2.5, 5, NA,
                                          0, 1.5, NA), 
                                        ncol = 3, byrow = TRUE), 
                           others = NA)
  forest_prod <- classify(EU_forests, 
                          rcl = matrix(c(2.5, 3.5, 1,
                                         0, 2.5, NA), 
                                       ncol = 3, byrow = TRUE), 
                          others = NA)
  
  # Save the different types of forest rasters
  writeRaster(forest_natural, "data/Nature_map_EU/Nature_EU_forest_natural.tif")
  writeRaster(forest_multi, "data/Nature_map_EU/Nature_EU_forest_multi.tif")
  writeRaster(forest_prod, "data/Nature_map_EU/Nature_EU_forest_prod.tif")
}

# Check if a 10km forest raster exists; if not, create it
if (!file.exists("data/Nature_map_EU/10km_forest.tif")) {
  forest_multi <- rast("data/Nature_map_EU/Nature_EU_forest_multi.tif")
  forest_prod <- rast("data/Nature_map_EU/Nature_EU_forest_prod.tif")
  forest_natural <- rast("data/Nature_map_EU/Nature_EU_forest_natural.tif")
  
  pixels <- EU_forests
  values(pixels) <- 1
  plot(pixels)
  
  # Extract the sum of binary raster values within each polygon
  forest_10km <- spp_template_shp |>
    mutate(forest_natural = exact_extract(forest_natural, spp_template_shp, fun = "sum"),
           forest_multi = exact_extract(forest_multi, spp_template_shp, fun = "sum"),
           forest_prod = exact_extract(forest_prod, spp_template_shp, fun = "sum"),
           pixels = exact_extract(pixels, spp_template_shp, fun = "sum"))
  
  # Save the forest data at 10km resolution as a shapefile
  st_write(forest_10km, "data/Nature_map_EU/10km_forest.shp")
  
  # Transform the data into a longer format and rasterize it
  forest_10km <- forest_10km |>
    pivot_longer(-c("ID", "geometry"))
  
  forest_10km_rast <- fasterize(forest_10km, raster(spp_template), field = "value", fun = "first", by = "name")
  writeRaster(forest_10km_rast, "data/Nature_map_EU/10km_forest.tif")
}

# Read in the 10km forest raster data
forest_10km_rast <- rast("data/Nature_map_EU/10km_forest.tif")
names(forest_10km_rast) <- c("natural", "multi", "prod", "pixels")

## Grab species data

# List all species distribution model (SDM) files
filelist_temp <- list.files("data/CurrentSDMs//", pattern = "*.tif$", recursive = TRUE)

# Function to split the string and keep the part after "__"
split_and_keep_after <- function(x) {
  split_string <- strsplit(x, split = "__")[[1]]
  return(substr(split_string[2], 1, nchar(split_string[2]) - 4))
}

# Apply the function to each element in the list to extract species names
names_sdms <- lapply(filelist_temp, split_and_keep_after)
names_sdms <- unlist(names_sdms)

### MAES (Mapping and Assessment of Ecosystems and their Services)

# Load habitat preferences data and filter for "Preferred" and "Suitable" habitats
MAES <- readRDS("data/EEA_habitatpreferences.rds") |>
  filter(typeasso %in% c("Preferred", "Suitable")) |>
  mutate(maes_label = recode(maes_label,
                             `Sparserly vegetated land` = "SparseVeg",
                             `Heathland and shrub` = "HeathlandShrub",
                             `Rivers and lakes` = "RiversLakes",
                             `Woodland and forest` = "WoodlandForest",
                             `Marine inlets and transitional waters` = "MarineTransitional"))

# Fix the naming of species to match the format above
MAES$speciesname <- gsub(" ", "_", MAES$speciesname, fixed = TRUE)
MAES_forest <- MAES |>
  filter(maes_label == "WoodlandForest")

# Join species names with forest habitat preferences and remove unnecessary columns
names_forests <- as.data.frame(names_sdms) |>
  rename(speciesname = names_sdms) |>
  left_join(MAES_forest) |>
  drop_na(maes_label) |> 
  dplyr::select(-c(season, typeasso)) |>
  unique()

# Filter the list of species SDM files to only include those matching the forest habitat preferences
filelist_temp_forest <- 
  as.data.frame(filelist_temp) |>
  mutate(speciesname = names_sdms) |>
  left_join(names_forests) |>
  drop_na()

# Load the SDM rasters, project them to match the template, and set non-habitat cells to NA
spp <- rast(paste0("data/CurrentSDMs/", filelist_temp_forest$filelist_temp))
spp <- project(spp, crs(spp_template))
spp[spp != 1] <- NA
names(spp) <- filelist_temp_forest$speciesname

# Combine forest and species data into a single raster stack
stack <- c(forest_10km_rast, spp)

# Organize planning units (PU) by combining raster data with spatial coordinates
pu_features_data <- 
  bind_cols(as_tibble(raster::as.data.frame(stack, xy = TRUE))) 

# Reshape data to longer format, keeping only necessary columns
pu_features_data <- pu_features_data |>
  pivot_longer(-c(x, y, natural, multi, prod, pixels))

##################### Threats by Land Production Intensity #################
## Species threats data

# Load species ID data and clean the species names
species_id <- readRDS("data/MinSpeciesToCover.rds") |>
  dplyr::select(speciesname, taxon_id)
species_id$speciesname <- gsub(" ", "_", species_id$speciesname, fixed = TRUE)

# Load species threat data from two different directives and combine them
threats_17 <- read_csv("data/HabitatsDirectiveSpeciesThreats.csv") |>
  dplyr::select(taxon_id, pressure_code)
threats_12 <- read_csv("data/BirdsDirectiveSpeciesThreats.csv") |>
  dplyr::select(taxon_id, pressure_code)
spp_threats <- rbind(threats_17, threats_12)

# Load land cover threats data and reshape it to wide format
lc_threats <- read_csv("data/pressures_threats.csv") |>
  rename(
    WoodlandForest_primary = G4Mconservation,
    WoodlandForest_multi = G4Mmultipurpose,
    WoodlandForest_prod = G4Mproduction) |>
  pivot_longer(-Code) |>
  drop_na(value)

lc_threats_wide <- lc_threats |>
  pivot_wider(values_from = value, 
              names_from = name) |>
  rename(pressure_code = Code) 

# Join species threats data with land cover threats, summarizing by species
spp_threats <- spp_threats |> 
  left_join(lc_threats_wide) |>
  dplyr::select(taxon_id, pressure_code, 
                WoodlandForest_primary,
                WoodlandForest_multi,
                WoodlandForest_prod) |>
  group_by(taxon_id) |>
  summarise(WoodlandForest_multi = sum(WoodlandForest_multi, na.rm = TRUE),
            WoodlandForest_prod = sum(WoodlandForest_prod, na.rm = TRUE))

# Calculate threat weightings based on quantiles
threat_weightings <- spp_threats |> 
  pivot_longer(-taxon_id) |>
  filter(value > 0) |>
  mutate(forest = "forest") |>
  group_by(forest) |>
  summarise(low_threat = quantile(value, 0.10),
            mid_threat = quantile(value, 0.50),
            high_threat = quantile(value, 0.90)) 

low_threat <- threat_weightings$low_threat[[1]]
mid_threat <- threat_weightings$mid_threat[[1]]
high_threat <- threat_weightings$high_threat[[1]]

# Adjust species threats based on the calculated threat weightings
spp_threats_weighted <- spp_threats |>
  pivot_longer(-taxon_id) |>
  mutate(
    value = case_when(
      value < low_threat ~ 1,
      value >= low_threat & value < mid_threat ~ 0.75,
      value >= mid_threat & value < high_threat ~ 0.25,
      TRUE ~ 0 # Default case, optional
    )
  ) |>
  pivot_wider(names_from = name, values_from = value)

# Calculate forest area adjustments based on the weighted threats
pu_forest_spp_data <- pu_features_data |>
  drop_na(value) |>
  filter(value >0) |>
  #mutate(value = replace_na(value,0)) |>
  filter(pixels > 19000) |>
  mutate(area_forest = (natural + multi + prod)) |>
  filter(area_forest > 0) |>
  rename(speciesname = name) |>
  left_join(species_id) |> 
  left_join(spp_threats_weighted) |>
  drop_na(taxon_id) |>
  dplyr::select(-c(speciesname)) |>
  mutate(natural_spp = natural*value,
         multi_spp = multi *value,
         prod_spp = prod*value) |>
  mutate(multi_spp_adj = multi * WoodlandForest_multi,
         prod_spp_adj = prod * WoodlandForest_prod) |> 
  group_by(taxon_id) |>
  # Group by species and calculate area of habitat and area of forest adjusted by threats
  mutate(aoh_spp = sum(natural_spp + multi_spp_adj + prod_spp_adj)) |>
  ungroup() |>
  mutate(aoh_spp = aoh_spp / pixels) 

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

## richness 
richness_natural <- pu_forest_spp_data |>
  group_by(x, y) |>
  summarise(richness_natural = sum(value, na.rm = T))|>
  ungroup() 

richness_multi <- pu_forest_spp_data |>
  mutate(value = value*WoodlandForest_multi) |>
  group_by(x, y) |>
  summarise(richness_multi = sum(value, na.rm = T)) |>
  ungroup() 

richness_prod <- pu_forest_spp_data |>
  mutate(value = value*WoodlandForest_prod) |>
  group_by(x, y) |>
  summarise(richness_prod = sum(value, na.rm = T)) |>
  ungroup() 


richness_management <- richness_natural |> 
  left_join(richness_multi) |> left_join(richness_prod) |>
  pivot_longer(-c(x,y)) |>
  drop_na(value) |>
  rename(richness = value) |>
  mutate(richness_importance = range01(richness))

richness_importance <- richness_management |>
  dplyr::select(-richness) |>
  pivot_wider(names_from = name, 
              values_from = richness_importance) 

plot(rast(richness_importance))

## rwr 
rwr_natural <- pu_forest_spp_data |>
  mutate(value = value/aoh_spp) |>
  group_by(x, y) |>
  summarise(rwr_natural = sum(value, na.rm = T))|>
  ungroup() 

rwr_multi <- pu_forest_spp_data |>
  mutate(value = value*WoodlandForest_multi/aoh_spp) |>
  group_by(x, y) |>
  summarise(rwr_multi = sum(value, na.rm = T)) |>
  ungroup() 

rwr_prod <- pu_forest_spp_data |>
  mutate(value = value*WoodlandForest_prod/aoh_spp) |>
  group_by(x, y) |>
  summarise(rwr_prod = sum(value, na.rm = T)) |>
  ungroup() 


rwr_management <- rwr_natural |> 
  left_join(rwr_multi) |> left_join(rwr_prod) |>
  pivot_longer(-c(x,y)) |>
  drop_na(value) |>
  rename(rwr = value) |>
  mutate(rwr_importance = range01(rwr))

rwr_importance <- rwr_management |>
  dplyr::select(-rwr) |>
  pivot_wider(names_from = name, 
              values_from = rwr_importance) 


forest_biodiversity <- rwr_importance |> left_join(richness_importance)


# Calculate forest species richness and importance based on adjusted forest area
forest_spp_richness <- pu_forest_spp_data |>
  drop_na() |>
  #filter(aoh_spp > 50) |>
  mutate(rw = (natural_spp + multi_spp_adj + prod_spp_adj) / aoh_spp,
         richness = (natural_spp + multi_spp_adj + prod_spp_adj)) |>
  group_by(x, y) |>
  summarise(rwr_importance = sum(rw, na.rm = TRUE),
            richness_importance = sum(richness, na.rm = TRUE),
            area_forest = mean(area_forest)/pixels,
            area_natural = mean(natural) / pixels,
            area_multi = mean(multi) / pixels,
            area_prod = mean(prod) / pixels) |> 
  ungroup() |> drop_na() |>
  mutate(rwr_importance_areaadj = range01(rwr_importance),
         richness_importance_areaaj = range01(richness_importance)) |>
  dplyr::select(-c(rwr_importance,richness_importance)) |>
  unique()

forest_biodiversity <- forest_spp_richness  |> left_join(forest_biodiversity)

plot(rast(forest_biodiversity))

write_csv(forest_biodiversity, "data/forest_spp_importance.csv")

writeRaster(rast(forest_biodiversity), "data/forest_spp_importance.tif")
