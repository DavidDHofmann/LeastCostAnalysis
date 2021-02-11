################################################################################
#### Example for Connectivity Analysis from Hofmann et al. 2021
################################################################################
# Description: This script exemplifies the modelling approach described Hofmann
# et al. 2021. The intention of this script is to help interested researchers to
# reproduce a similar analysis for their own data. Note that we assume that the
# integrated step selection analysis (iSSF) has already been completed. For
# further detail on iSSF analysis, we highly encourage you to read Fieberg et
# al. (2020), the code examples given in Muff et al. 2020, and the "amt" package
# vignette.

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)       # To handle spatial data
library(rgeos)        # To manipulate spatial data
library(gdistance)    # To conduct least-cost analysis
library(NLMR)         # To simulate spatial layers
library(viridis)      # For nice colors in visualizations
library(rasterVis)    # For plotting
library(arrangements) # To identify all pairwise combinations
library(dismo)        # To compute voronoi polygons

# Specify a seed for reproducability
set.seed(1289)

################################################################################
#### Some Useful Functions
################################################################################
# Function to make a color transparent
colTrans <- function(color, percent = 50){
  rgb <- col2rgb(color)
  col_new <- rgb(
      rgb[1]
    , rgb[2]
    , rgb[3]
    , max = 255
    , alpha = (100 - percent) * 255 / 100
  )
  return(col_new)
}

# Function to normalize a raster to a range between 0 and 1
normalizeMap <- function(x){
  (x - minValue(x))/(maxValue(x) - minValue(x))
}

# Function to truncate permeability surface and remove lower and upper
# percentiles (similar to Squires et al. 2010)
truncRaster <- function(x, percentile = 0.01){
  upper <- quantile(values(x), 1 - percentile, na.rm = TRUE)
  lower <- quantile(values(x), percentile, na.rm = TRUE)
  values(x)[values(x) > upper] <- upper
  values(x)[values(x) < lower] <- lower
  return(x)
}

# Function to create centroids that lie within a polygon
gCentroidWithin <- function(pol){
  pol$.tmpID <- 1:length(pol)
  initialCents <- gCentroid(pol, byid = T)
  centsDF <- SpatialPointsDataFrame(initialCents, pol@data)
  centsDF$isCentroid <- TRUE
  centsInOwnPoly <- sapply(1:length(pol), function(x){
    gIntersects(pol[x,], centsDF[x, ])
  })
  if (all(centsInOwnPoly) == TRUE){
        return(centsDF)
  } else {
    newPoints <- SpatialPointsDataFrame(
        gPointOnSurface(
            pol[!centsInOwnPoly, ]
          , byid = T
        )
      , pol@data[!centsInOwnPoly,]
    )
    newPoints$isCentroid <- FALSE
    centsDF <- rbind(centsDF[centsInOwnPoly, ], newPoints)
    centsDF <- centsDF[order(centsDF$.tmpID), ]
    centsDF@data <- centsDF@data[, - which(names(centsDF@data) == ".tmpID")]
    return(centsDF)
  }
}

################################################################################
#### Simulation of Covariates and Prediction of Permeability Surface
################################################################################
# For simplicity, we will only simulate three covariates: Shrub/Grassland cover
# (continuous), Human Influence (continuous), and Water (binary).
n <- 300
shrubs <- nlm_gaussianfield(ncol = n, nrow = n)
humans <- nlm_gaussianfield(ncol = n, nrow = n)
water <- nlm_gaussianfield(ncol = n, nrow = n, autocorr_range = 50)
water <- water < 0.45

# Make nice layernames
names(shrubs) <- "Shrubs"
names(humans) <- "Humans"
names(water)  <- "Water"

# Plot the covariates
levelplot(stack(shrubs, humans, water))

# Now we use the coefficients from the iSSF model to predict a permeability
# surface. Let's assume that individuals avoid water and humans but prefer
# shrubs/grassland.
perm <- exp(-0.6 * water - 0.4 * humans + 0.4 * shrubs)

# Truncate predicted values (if desired)
perm <- truncRaster(perm, percentile = 0.01)

# Plot the resulting surface map
plot(perm, col = viridis(50), main = "Permeability Surface")

# We can now turn the permeability surface into a graph using the gdistance
# package. Here, we'll take into account each cell's 8 neighbours. However, you
# could also go for 4 or 16 neighbors only (see gdistance vignette)
trans <- transition(perm, transitionFunction = mean, directions = 8)

# Apply geographic correction
trans <- geoCorrection(trans, type = "c", multpl = FALSE)

################################################################################
#### Selection of Source Points - Approach 1: Along Map Border
################################################################################
# Retrieve polygon of map extent
ext1 <- as(extent(perm), "SpatialLines")

# Distribute desired number of source points along border
points1 <- spsample(ext1, n = 20, type = "regular")

# Visualize the points
plot(perm, col = viridis(50))
plot(points1, add = T, col = "red", pch = 20, cex = 2)

################################################################################
#### Selection of Source Points - Approach 2: Inside Protected Areas
################################################################################
# Alternatively, we could distribute source points within protected areas. Let's
# therefore simulate some protected areas. For this, well create some random
# shapes along the map borders.
ext1 <- as(extent(perm), "SpatialPolygons")
ext2 <- gBuffer(ext1, width = -50)
ext <-  gDifference(ext1, ext2)
rand <- spsample(ext, n = 150, type = "random")
rand <- voronoi(rand, ext = extent(ext))
prot <- rand[sample(x = 1:nrow(rand), size = 10), ]
prot$ID <- LETTERS[1:nrow(prot)]

# Visualize the protected areas
plot(perm, col = viridis(50))
plot(prot, add = T, border = "red")

# Overlay the study area with a regular grid of source points, but keep only
# those that lie within protected areas.
points2 <- spsample(ext1, type = "regular", n = 200)
keep <- gContains(prot, points2, byid = TRUE) %>% rowSums(.)
points2 <- points2[keep > 0, ]

# Visualize
plot(perm, col = viridis(50))
plot(prot, add = T, border = "red")
plot(points1, add = T, col = "red", pch = 20, cex = 2)

# Find protected areas which do not yet contain any source points yet
missing <- is.na(over(prot, points2))

# Create centroid for those protected areas
centroids <- gCentroidWithin(prot[missing, ])
points2 <- rbind(SpatialPoints(points2), SpatialPoints(centroids))

# Visualize them on top of the permeability surface
plot(perm, col = viridis(50))
plot(prot, add = T, border = "red")
plot(points1, add = T, col = "red", pch = 20, cex = 2)

################################################################################
#### Calculate Least-Cost Paths
################################################################################
# Select the points you want to use
# points <- points1
points <- points2

# Identify all possible combinations of start and end points
combis <- as.data.frame(combinations(length(points), 2))
names(combis) <- c("Origin", "Destin")

# How many combinations are there?
nrow(combis)

# Find the shortest paths for all possible connections for the first set of
# source points.
paths <- lapply(1:length(unique(combis$Origin)), function(x){
  shortestPath(trans
    , origin  = points[x, ]
    , goal    = points[combis$Destin[combis$Origin == x], ]
    , output  = "SpatialLines"
  )
}) %>% do.call(rbind, .)

# Visualize paths
plot(perm, col = viridis(50))
plot(paths, add = T, col = "white")
plot(prot, add = T, border = "red")
plot(points, add = T, col = "red", pch = 20, cex = 2)

################################################################################
#### Calculate Least-Cost Corridors
################################################################################
# To calculate least-cost corridors, we first need to compute cumulative cost
# maps for each of the source points
cost <- lapply(1:length(points), function(x){
  accCost(trans, points[x, ])
}) %>% stack()

# Can visualize the cumulative cost map for each source point
plot(cost[[1]], col = viridis(50))
plot(points[1, ], add = T, col = "red", pch = 20, cex = 5)

# To identify a corridor between two maps, we need to sum their cost-maps and
# identify the cheapest set of pixels. Let's write a function to do this.
getCorr <- function(x, y, percent = 0.01){
  summed <- x + y
  threshold <- minValue(summed) * (1 + percent)
  corr <- reclassify(summed, c(threshold, Inf, NA))
  corr <- normalizeMap(corr)
  corr <- reclassify(corr, c(NA, NA, 1), right = FALSE)
  return(corr)
}

# Identify a corridors
corrs <- lapply(1:nrow(combis), function(x){
  getCorr(cost[[combis$Origin[x]]], cost[[combis$Destin[x]]], percent = 0.02)
}) %>% stack()
corrs <- corrs * (-1) + 1
corrs <- calc(corrs, sum)

# Visualize corridors
plot(corrs, col = magma(50))
plot(points, add = T, col = "red", pch = 20, cex = 2)
plot(paths, add = T, col = colTrans("white", percent = 80))

################################################################################
#### Visualize everything
################################################################################
# Visualize LCPs and LCCs
plot(corrs, col = magma(50))
plot(points, add = T, col = "red", pch = 20, cex = 2)
plot(prot, add = T, border = colTrans("red", percent = 20))
plot(paths, add = T, col = colTrans("white", percent = 80))
