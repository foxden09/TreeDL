# Libraries
library(terra) 
library(ggplot2)
library(sf)
library(lidR)
library(rgl)
library(raster)

## PREPROCESSING ---------------------------------------------------------------

# Read and filter LiDAR data
las_on <- readLAS("LiDAR/PSU_072020_580800e_0148300n_leafon.las", filter = "-keep_first -drop_z_below 0")
las_off <- readLAS("LiDAR/PSU_032020_580800e_0148300n_leafoff.las", filter = "-keep_first -drop_z_below 0")

# Set coordinate system of leaf-off to match leaf-on 
crs_leaf_on <- st_crs(las_on)
st_crs(las_off) <- crs_leaf_on

# Crop to make processing more manageable
xmin <- 580800
xmax <- 581000
ymin <- 147400
ymax <- 147600
cropped_las_on <- clip_rectangle(las_on, xmin, ymin, xmax, ymax)
cropped_las_off <- clip_rectangle(las_off, xmin, ymin, xmax, ymax)

## HEIGHT NORMALIZATION AND FILTERING ---------------------------------------------------------------

# Normalize heights + denoise
#las_on_denoised <- filter_poi(cropped_las_on, Classification != 1)
# codes contained: 1 (unclassified) 2 (ground) 5 (high veg) 7 
nlas_on <- normalize_height(las_on_denoised, knnidw())
nlas_off <- normalize_height(las_off_denoised, knnidw())

# Plot height distribution
hist(nlas_off@data$Z, breaks = 100, main = "Distribution of Heights", xlab = "Height (meters)")

# Keep only high vegetation
las_on_5 <- filter_poi(nlas_on, Classification == 5)
las_off_5 <- filter_poi(nlas_off, Classification == 5)

## CANOPY MODELING ---------------------------------------------------------------

# Create CHM (canopy height models)
chm_on <- rasterize_canopy(nlas_on, res = 0.5, algorithm = pitfree(subcircle = 0.15), pkg = "terra")
chm_off <- rasterize_canopy(nlas_off, res = 0.5, algorithm = pitfree(subcircle = 0.15), pkg = "terra")
w <- matrix(1, 3, 3)
smoothed <- terra::focal(chm_on, w, fun = mean, na.rm = TRUE)

# Locate trees 
f <- function(x) {x * 0.1 + 5}
f2 <- function(x) {
  y <- 2.6 * (-(exp(-0.08*(x-2)) - 1)) + 3
  y[x < 10] <- 2
  y[x >= 2 & x < 10] <- 3
  y[x > 20] <- 5
  return(y)
}
ttops_on <- locate_trees(nlas_on, lmf(5, hmin = 10))
ttops_off_constant <- locate_trees(nlas_off, lmf(5, hmin = 10))
ttops_off_f <- locate_trees(nlas_off, lmf(f, hmin = 10))
ttops_off_f2 <- locate_trees(nlas_off, lmf(f2, hmin = 10))

# Plot tree tops
plot(chm_on, col = height.colors(50))
plot(chm_off)
plot(sf::st_geometry(ttops_on), add = TRUE, pch = 4, col = "red", lwd =1.52)
plot(sf::st_geometry(ttops_off), add = TRUE, pch = 4, col = "white", lwd =1.52)
plot(sf::st_geometry(ttops_off_f), add = TRUE, pch = 4, col = "red", lwd =1.52)
plot(sf::st_geometry(ttops_off_f2), add = TRUE, pch = 4, col = "white", lwd =1.52)
plot(sf::st_geometry(ttops_off_constant), add = TRUE, pch = 4, col = "white", lwd =1.52)
  
## SEGMENTATION ---------------------------------------------------------------

# Segment
ttops_off <- ttops_off_f
algo_dalponte <- dalponte2016(smoothed, ttops_off, th_tree = 10, th_seed = 0.65, th_cr = 0.2, max_cr = 20)
algo_silva <- silva2016(chm_on, ttops_off)
algo_li <- li2012(hmin = 10)
trees_on <- segment_trees(nlas_on, algo_silva) 
trees_off <- segment_trees(nlas_off, silva2016(chm_off, ttops_off)) 
#trees_on@data$treeID <- trees_on@data$tree_segmentation
#head(trees_on@data$treeID, 1000)

x <- plot(trees_on, bg = "white", size = 4, color = "treeID", height = TRUE, treetops = TRUE)
add_treetops3d(x, ttops_off)

## SAVE PROCESSED DATA ---------------------------------------------------------------

# Convert to ArcGIS-compatible form
trees_on_sf <- st_as_sf(segment_trees(nlas_on, dalponte2016(chm_on, ttops_off)))
st_write(trees_on_sf, "trees_on_sf_dalponte.shp", layer_options = "SHPT=POINTZ")


