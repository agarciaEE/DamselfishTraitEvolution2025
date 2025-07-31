### distance matrices
# Load necessary libraries
library(geosphere)
library(oce)
library(raster)
library(gdistance)
library(BioGeoBEARS)
library(cladoRcpp)
library(tidyverse)

setwd("~/Unil/Research/5_Damselfish_evo/")

# Define the centroids of each geographic area
centroids <- data.frame(
  Longitude = c(70, 120, -150, -90, -20),  # center longitudes
  Latitude = c(-20, -5, 10, 0, 0),  # center latitudes
  row.names = c("Indian Ocean", "Indo-Australian Archipelago", "Central Pacific", "East Pacific", "Atlantic")
)

# Print the table
print(centroids)
# plot
maps::map("world")
points(centroids[,1:2], pch = 19, cex = 1.2)

#################
# Function to go from current lonlat coords to the corresponding lonlat coords at a specific Age
#################
# given current lon lat input origin
xyA <- as.numeric(centroids[1,])
# given current lon lat input end
xyB <- as.numeric(centroids[3,])
# given target Age to locate the equivalent coordinates
target_age <- 7
# given the file path to coral paleo data
paleodata_file<- "./BioGeo/Theo_paleo/Matrix_steps.rda"
# given the dir path to distance matrices files between Ages
path_ij <- "./BioGeo/Theo_paleo/Distance_matrices_external"
# given the prefix of distance matrices files between Ages
prefix_ij <- "Distance_"
# given the dir path to distance matrices files within Ages
path_ii <- "./BioGeo/Theo_paleo/Distance_matrices_internal"
# given the prefix of distance matrices files within Ages
prefix_ii <- "Distance_"
# select type of distance
dist_type <- c("geodist", "Eculidean")
# select type of output
output <- c("all", "coords", "dists", "lastcoords", "lastdists")

get_paleodist_fromA2B <- function(xyA, xyB, target_age,
                                  paleodata_file = "./BioGeo/Theo_paleo/Matrix_steps.rda",
                                  path_ij = "./BioGeo/Theo_paleo/Distance_matrices_external",
                                  prefix_ij = "Distance_",
                                  path_ii = "./BioGeo/Theo_paleo/Distance_matrices_internal",
                                  prefix_ii = "Distance_",
                                  dist_type = c("geodesic", "Euclidean"),
                                  output = c("all", "coords", "dists", "lastcoords", "lastdist")){
  
  output <- output[1]
  
  # load coral paleo data
  load(paleodata_file)
  
  # subset coral reef coordinates in current Age
  # value of 3 represents temperate coral reefs
  # value of 30 represents tropical coral reefs
  coords_0 <-  Matrix_steps[Matrix_steps[,"Coralb_0"] %in% c(3, 30), c("x", "y")]
  
  # first find the closest coral coordinates of the data set to locations A and B
  idxA_i <- which.min(sqrt((xyA[1] - coords_0[,1])^2 + (xyA[2] - coords_0[,2])^2))
  xyA <- as.numeric(coords_0[idxA_i,])
  
  idxB_i <- which.min(sqrt((xyB[1] - coords_0[,1])^2 + (xyB[2] - coords_0[,2])^2))
  xyB <- as.numeric(coords_0[idxB_i,])
  
  if (output %in% c("all", "dists")){
    i = 0 # Current Age
    dm_file_ii <- file.path(path_ii, paste0(prefix_ii, i, ".rda"))
    # load distance matrix of Age i
    message("Loading distance matrix of Age ", i, "...")
    if (!file.exists(dm_file_ii)){
      stop("File '", dm_file_ii, "' not found.\nDistance matrix cannot be loaded.")
    }
    load(dm_file_ii)
    
    # get distance from A to B in Age i
    if ("geodesic" %in% dist_type) { geod <- oce::geodDist(xyA[1], xyA[2], xyB[2], xyB[2]) } else { geod <- NULL }
    if ("Euclidean" %in% dist_type) { d <- dist1[idxA_i, idxB_i] } else { d <- NULL }
  }
  
  if (output %in% c("all", "coords")){
    # create matrix to trace evolution of coordinates back in time
    xyA_trace <- xyB_trace <- matrix(NA, nrow = target_age+1, ncol = 2, 
                                     dimnames = list(as.character(0:target_age),
                                                     c("x", "y")))
    # assign first entry to current Age 0
    xyA_trace[1,] <- xyA
    xyB_trace[1,] <- xyB
  }

  # loop til get to the target age
  for (i in 1:target_age){
  
    # get distance matrices file paths
    dm_file_ij <- file.path(path_ij, paste0(prefix_ij, i, "_", i-1, ".rda"))
    dm_file_ii <- file.path(path_ii, paste0(prefix_ii, i, ".rda"))
    
    # load distance matrix from Age i-1 to Age i
    message("Loading distance matrix from Age ", i-1, " to Age ", i, "...")
    if (!file.exists(dm_file_ij)){
      stop("File '", dm_file_ij, "' not found.\nDistance matrix cannot be loaded.")
    }
    load(dm_file_ij)
    
    # get row indexes for A and B in Age i
    idxA_i <- which.min(dist2[, idxA_i])
    idxB_i <- which.min(dist2[, idxB_i])
    
    # get coordinates in Age i
    coords_i <- Matrix_steps[Matrix_steps[, paste0("Coralb_", i)] %in% c(3, 30), c("x", "y")]
  
    if (output %in% c("all", "coords")){
      # get equivalent coordinates of A and B
      xyA_i <- as.numeric(coords_i[idxA_i,])
      xyB_i <- as.numeric(coords_i[idxB_i,])
      
      # add coordinates of A and B in Age i to the tracking matrix
      xyA_trace[i+1,] <- xyA_i
      xyB_trace[i+1,] <- xyB_i
    }
    else if (output == "lastcoords" & i == target_age) {
      # get equivalent coordinates of A and B
      xyA_i <- as.numeric(coords_i[idxA_i,])
      xyB_i <- as.numeric(coords_i[idxB_i,])
    }
    
    if (output %in% c("all", "dists", "lastdist")){
      # load distance matrix of Age i
      message("Loading distance matrix of Age ", i, "...")
      if (!file.exists(dm_file_ii)){
        stop("File '", dm_file_ii, "' not found.\nDistance matrix cannot be loaded.")
      }
      load(dm_file_ii)
      
      if (output == "lastdist" & i == target_age) {
        # get distance from A to B in Age i
        if ("geodesic" %in% dist_type) { geod <- oce::geodDist(xy0_i[1], xy0_i[2], xy1_i[1], xy1_i[2]) } else {geod <- NULL }
        if ("Euclidean" %in% dist_type) { d <- dist1[idx0_i, idx1_i] } else { d <- NULL }
      } 
      else if (output != "lastdist") {
        # get distance from A to B in Age i
        if ("geodesic" %in% dist_type) { geod <- c(geod, oce::geodDist(xy0_i[1], xy0_i[2], xy1_i[1], xy1_i[2])) }
        if ("Euclidean" %in% dist_type) { d <- c(d, dist1[idx0_i, idx1_i]) }
      }
    }
  }
  dists <- list(geodesic = geod, Euclidean = d)
  dists <- dists[sapply(dists, length) > 0]
  
  if (output == "all") {
    obj <- list(coordsA = xyA_trace, 
                coordsB = xyB_trace,
                dists = dists)
  } 
  else if (output == "coords") {
    obj <- list(coordsA = xyA_trace, 
                coordsB = xyB_trace)
  }
  else if (output %in% c("dists", "lastdist")) {
    obj <- dists
  }
  else if (output == "lastcoords") {
    obj <- list(coordA = xyA_i, 
                coordB = xyB_i)
  }

  return(obj)
}

##
## Function to get equivalent coordinates from species Age
#########################
get_paleocoords <- function(xy.coords, target_age, res = 5,
                            extrapolate = FALSE, res.extra = 5,
                            paleodata_file = "./BioGeo/Theo_paleo/Matrix_steps.rda",
                            path_ij = "./BioGeo/Theo_paleo/Distance_matrices_external",
                            prefix_ij = "Distance_"){
  
  if (is.numeric(xy.coords)) {
    xy.coords <- data.frame(t(as.data.frame(xy.coords)))
  }
  
  if (ncol(xy.coords) > 0){
    extra_cols <- xy.coords[,3:ncol(xy.coords)]
    extra_colnames <- colnames(xy.coords)[3:ncol(xy.coords)]
    extra.cols <- TRUE
  } else {
    xy.coords$extracol = 1
  }
  # load coral paleo data
  load(paleodata_file)
  # subset coral reef coordinates in current Age
  # value of 3 represents temperate coral reefs
  # value of 30 represents tropical coral reefs
  coords_0 <-  Matrix_steps[Matrix_steps[,"Coralb_0"] %in% c(3, 30), c("x", "y", "Coralb_0")]
  rownames(coords_0) <- 1:nrow(coords_0)
  # first find the closest coral coordinates of the data set to locations base on a maximum distance (res)
  coords_0 <- ecospat::ecospat.sample.envar(coords_0, 1:2, 1:3, xy.coords, 1:2, 3, res = res)
  #xy.coords <- xy.coords[!is.na(xy.coords[,3]),] # remove non-assigned locations (this can be fine-tune modifying res)
  ncoords <- nrow(coords_0[,])
  # get indexes
  idx_i <- as.integer(rownames(coords_0[,]))
  extra_cols <- coords_0[,4]

  # loop til get to the target age
  for (i in 1:target_age){
    
    # get distance matrices file paths
    dm_file_ij <- file.path(path_ij, paste0(prefix_ij, i, "_", i-1, ".rda"))

    # load distance matrix from Age i-1 to Age i
    message("Loading distance matrix from Age ", i-1, " to Age ", i, "...")
    if (!file.exists(dm_file_ij)){
      stop("File '", dm_file_ij, "' not found.\nDistance matrix cannot be loaded.")
    }
    load(dm_file_ij)
    
    coords_i <- Matrix_steps[Matrix_steps[, paste0("Coralb_", i)] %in% c(3, 30), c("x", "y")]

    idx_i <- sapply(1:ncoords, function(n) { 
      d <- dist2[,idx_i[n]]
      d[!is.finite(d)] = NA
      id <- suppressWarnings(which(d == min(d, na.rm = T)))
      if (length(id)>0){
        dc <- sqrt((coords_0[n,1] - coords_i[id,1])^2 + (coords_0[n,2] - coords_i[id,2])^2)
        id[which.min(dc)]
      } else { NA }
      })
  }
  # get equivalent coordinates of Age i-1 in Age i
  xy.coords <- coords_i[idx_i,]

  if (extra.cols){
    xy.coords <- cbind(xy.coords, extra_cols)
    colnames(xy.coords)[3:ncol(xy.coords)] <- extra_colnames
  }
  if (extrapolate){
    xy.coords <- xy.coords[rowSums(is.na(xy.coords[,1:2])) == 2,]
    xy.coords <- ecospat::ecospat.sample.envar(xy.coords, 1:2, 1:2, na.exclude(xy.coords), 1:2, 3, res = res.extra)
  }
  
  return(na.exclude(xy.coords))
}

#########################
# get paleo maps and ages
#########################
age_lim <- 60

## From Cao_etal_2017
#########################
path <- "./BioGeo/Cao_etal_2017_BG_Supplement/bg-2017-94-supplement/Reconstructed_Paleogeog_Matthews2016_Revised_402_2Ma_GeoTiffs"
paleo_map_files <- list.files(path, pattern = "\\.tiff$")
ages <- as.numeric(sapply(paleo_map_files, function(i) sub("Ma.tiff", "", strsplit(i, "_")[[1]][3])))
paleo_map_files <- paleo_map_files[order(ages)]
ages <- sort(ages)
# paleo_map_all <- stack(lapply(paleo_map_files, function(i) raster(paste0(path,i))))
# plot(paleo_map_all)
paleo_map_files <- paleo_map_files[ages < age_lim]
ages <- ages[ages < age_lim]

# get maps
paleo_map_ <- stack(lapply(paleo_map_files, function(i) raster(file.path(path,i))))
paleo_map[paleo_map > 0] = NA
paleo_map[paleo_map == 0] = 1
current <- raster("~/Unil/Research/1_NicheCompetition/snakemake_wd/data/environmental_data/bedtemp.tif")
current[!is.na(current)] = 1
current <- raster::resample(current, paleo_map, method = "ngb")
paleo_map <- stack(current, paleo_map)
ages <- c(0, ages)
names(paleo_map) <- paste0(ages, "Ma")

# write time periods
time_periods <- sapply(1:(length(ages)), function(i) mean(c(ages, last(ages) + last(ages) - ages[length(ages)-1])[i:(i+1)]))
write.table(time_periods, file = "./input/damsel_biogeo_timeperiods.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


paleo_map_orig <- paleo_map # save original
# resample to reduce resolution and speed up computation
# res.template <- raster(ncol = 180, nrow = 90, ext = extent(paleo_map))
# paleo_map_res <- raster::resample(paleo_map, res.template, method = "ngb")
# plot(paleo_map[[7]])

paleo_map <- aggregate(paleo_map, 20)
plot(paleo_map_agg20)

# crop pole areas
ext <- c(-180, 180, -70, 70)
paleo_map <- crop(paleo_map, ext)
plot(paleo_map)

########################

## From Theo (Descombes et al. 2018)
###################################

path <- "./BioGeo/Theo_paleo/Coral_rasters"
paleo_map_files <- list.files(path, pattern = "\\.asc$")
paleo_map_files <- grep("Coral", paleo_map_files, value = TRUE)
ages <- as.numeric(sapply(paleo_map_files, function(i) sub(".asc", "", strsplit(i, "_")[[1]][2])))
paleo_map_files <- paleo_map_files[order(ages)]
ages <- sort(ages)
# paleo_map_all <- stack(lapply(paleo_map_files, function(i) raster(paste0(path,i))))
# plot(paleo_map_all)
paleo_map_files <- paleo_map_files[ages < age_lim]
ages <- ages[ages < age_lim]

# get maps
paleo_map_all <- stack(lapply(paleo_map_files, function(i) raster(file.path(path,i))))
load(paleodata_file) # use matrix_steps to get current map
current <- NINA::raster_projection(Matrix_steps[, c("x", "y", "Coralb_0")])
paleo_map_all <- stack(current, paleo_map_all)
ages <- c(0, ages)
#plot(paleo_map)
names(paleo_map_all) <- paste0(ages, "Ma")

par(mfrow = c(1,2))
plot(paleo_map_Cao$X53Ma)
plot(paleo_map_all$X53Ma)

# subset paleo maps selecting the most important Ages based on epochs
epochs <- rev(setNames(c(66, 56, 33.9, 23.03, 5.33, 2.58, 0),
                    c("Paleocene", "Eocene", "Oligocene", "Miocene", 
                      "Pliocene", "Pleistocene", "Holocene")))
epochs <- epochs[epochs <= max(ages)]
paleo_map <- paleo_map_all[[paste0("X",round(epochs), "Ma")]]
ages <- round(epochs)

# write time periods
#time_periods <- sapply(1:(length(ages)), function(i) mean(c(ages, last(ages) + last(ages) - ages[length(ages)-1])[i:(i+1)]))
write.table(ages[-1], file = "./BioGeo/input/damsel_biogeo_timeperiods_Descombes.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# transform values in paleo raster stack
paleo_map_orig <- paleo_map # save original
plot(paleo_map)
paleo_map[paleo_map==1]<-NA # land values to NA
paleo_map[paleo_map==10]<-NA # land tropical values to NA
paleo_map[paleo_map==2]<-10 # temperate deep waters (<-1000) to 10
paleo_map[paleo_map==20]<-20 # tropical deep waters  (<-1000) to 20 
paleo_map[paleo_map==3]<-20 # temperate coral habitats (>-1000) to 20 
paleo_map[paleo_map==30]<-30 # tropical coral habitats (>-1000) to 30 
plot(paleo_map)

#################################
# get data points for each region
#################################
# get data.frame
paleo_data <- as.data.frame(paleo_map, xy = T)

# remove common land masses
paleo_data <- paleo_data[rowSums(!is.na(paleo_data[,-c(1:2)])) > 0,]

# get distances from each location to each region
dists <- sapply(1:nrow(centroids), function(i) 
                sapply(1:nrow(paleo_data), function(j) 
                       geodDist(centroids[i,1], centroids[i,2],
                                paleo_data[j,1], paleo_data[j,2])))

colnames(dists) <- rownames(centroids)

# get neighboring locations to centroids of each region
reg_locs <- apply(dists, 2, function(i) which(i < 5000))
names(reg_locs) <- rownames(centroids)

# scale distances to max dist of each region
dists <- apply(dists, 2, function(i) 1-(i/max(i, na.rm = T)))

# get approximate geograhical areas (great circle data points) of each time period
reg_age_coords <- list()
for (age in names(paleo_data)[-c(1:2)]){
  reg_age_coords[[age]] <- list()
  for (region in names(reg_locs)){
    reg_age_coords[[age]][[region]] <- na.exclude(cbind(paleo_data[reg_locs[[region]],c("x", "y", age)], 
                                                        dist = dists[reg_locs[[region]], region]))[,c(1:2, 4)]
  }
  reg_age_coords[[age]] <- plyr::ldply(reg_age_coords[[age]], .id = "region")
}

plot(paleo_map_res[[age]])     
points(reg_age_coords[[age]][,2:3] , pch = 19, cex = 1.2, 
       col = rep(viridis::viridis(5))[as.factor(reg_age_coords[[age]][,1])])

num_areas <- length(CB2013_regions)
n_points <- 50
area_code <- setNames(names(CB2013_regions),
                      LETTERS[1:num_areas])
transLayerList <- list()
dist_matrix <- list()
age <- names(paleo_map)[1]
# for each time period
for (age in names(paleo_map)){
  
  # Create an empty distance matrix
  dist_matrix[[age]] <- matrix(0, nrow = num_areas, ncol = num_areas, 
                               dimnames = list(rownames(centroids), rownames(centroids)))
  
  tr.layer = paleo_map[[age]]
  tr.layer[tr.layer == 0] = NA
  tr.layer[!is.na(tr.layer)] = 1
  tr.layer = gdistance::geoCorrection(gdistance::transition(tr.layer, function(x) mean(x), 8))
  transLayerList[[age]] <- tr.layer
  
  plot(paleo_map[[age]])
  for (i in 1:(num_areas - 1)) {
    area1 <- rownames(centroids)[i]
    xy1 <- as.matrix(centroids[i,])
    coordsA1 <- reg_age_coords[[age]][reg_age_coords[[age]]$region == area1,2:4]
    samp_coordsA1 <- coordsA1[sample(1:nrow(coordsA1), n_points, prob = coordsA1[,3]),1:2]
    
    for (j in (i + 1):num_areas) {
      area2 <- rownames(centroids)[j]
      xy2 <- as.matrix(centroids[j,])
      coordsA2 <- reg_age_coords[[age]][reg_age_coords[[age]]$region == area2,2:4]
      samp_coordsA2 <- coordsA2[sample(1:nrow(coordsA2), n_points, prob = coordsA2[,3]),1:2]
      
      #dist <- sapply(1:nrow(samp_coordsA1), function(o)
      #                    sapply(1:nrow(samp_coordsA2), function(d)
      #                       tryCatch(lengthLine(tryCatch(gdistance::shortestPath(tr.layer, 
      #                                                                            as.matrix(samp_coordsA1[o,]), 
      #                                                                            as.matrix(samp_coordsA1[d,]),
      #                                                                            output="SpatialLines"), 
      #                                                    error = function(e) NA)), 
      #                                error = function(e) NA) / 1000))
      
      #hist(dist)
      #dist <- mean(dist, na.rm = T)
      
      sp <- sapply(1:n_points, function(i) tryCatch(gdistance::shortestPath(tr.layer, 
                                                                            as.numeric(coordsA1[sample(1:nrow(coordsA1), 1, prob = coordsA1[,3]),1:2]), 
                                                                            as.numeric(coordsA2[sample(1:nrow(coordsA2), 1, prob = coordsA2[,3]),1:2]), 
                                                                            output="SpatialLines"),
                                                    error = function(e) NA))
      
      sp <- sp[!sapply(sp, is.na)]
      lapply(sp, lines)
      dist <- unlist(lapply(sp, lengthLine)) / 1000
      dist <- mean(dist, na.rm = TRUE)
      cat(area1, "-", area2, ":", dist, "\n")
      
      dist_matrix[[age]][i, j] <- dist
      dist_matrix[[age]][j, i] <- dist_matrix[[age]][i, j]  # Distance matrix is symmetric
      
    }
  }
  rownames(dist_matrix[[age]]) <- colnames(dist_matrix[[age]]) <- sapply(colnames(dist_matrix[[age]]), function(i) names(area_code[area_code == i]))
  
}
#################################

#######################################################################
# use pomcacentridae CB2013 classified distributions to get data points
#######################################################################
CB2013_regions <- setNames(c("IO", "IAA", "CPO", "EPO", "AO", "TS"), 
                           c("Indian", "Indo-Australian Archipelago", "Central Pacific", 
                             "East Pacific", "Atlantic", "Tethys"))

num_areas <- length(CB2013_regions)
area_code <- setNames(CB2013_regions,
                      LETTERS[1:num_areas])

load("./BioGeo/input/CB2013_regions_spdf.rda")
pomacentridae_occ_dis_CB2013_regs <- read.csv("./BioGeo/input/pomacentridae_dis_occ_CB2013_regs.csv")

CB2013reg_coords <- unique(pomacentridae_occ_dis_CB2013_regs[,c(1:2, 4)])
CB2013reg_coords$CB2013[CB2013reg_coords$CB2013 == "MS"] <- "TS"

# check occurrences on CB2013 regions map
plot(meow_regions_sp, col = rep(viridis::viridis(6))[factor(meow_regions_sp$CB2013,
                                                            levels = CB2013_regions)], border = NA)
maps::map("world", add = T, col ="grey90", fill = T)
points(CB2013reg_coords[,1:2], pch = 21, 
       bg = rep(viridis::viridis(num_areas))[factor(CB2013reg_coords$CB2013, levels = CB2013_regions)])

# get coordinates of equivalent CB2013 regions in the past
paleo_CB2013_reg_coords <- list()
for (age in ages) {
  if (age == 0){
    df <- CB2013reg_coords
    df$CB2013 <- factor(df$CB2013, levels = CB2013_regions)
    paleo_CB2013_reg_coords[[paste0("X", age, "Ma")]] <- df
  } else {
    df <- get_paleocoords(CB2013reg_coords, target_age = age, res = 5, extrapolate = F, res.extra = 10)
    df$CB2013 <- factor(df$CB2013, levels = CB2013_regions)
    plot(paleo_map[[paste0("X", age, "Ma")]])
    points(df[,1:2], pch = 21, 
           bg = rep(viridis::viridis(num_areas))[df$CB2013])
    paleo_CB2013_reg_coords[[paste0("X", age, "Ma")]] <- df
  }
}

n_points <- 100

shift_line_longitudes <- function(spline) {
  coords <- spline@lines[[1]]@Lines[[1]]@coords
  coords[coords[,1] > 180, 1] = coords[coords[,1] > 180,1] - 360
  spline@lines[[1]]@Lines[[1]]@coords <- coords
  spline
}

transLayerList <- list()
dist_matrix <- list()
# for each time period
for (age in names(paleo_map)){
  message("Computing distance between CB20213 regions ", sub("X", "", age), "...")
  reg_coords_age <- paleo_CB2013_reg_coords[[age]]
  
  # Create an empty distance matrix
  dist_matrix[[age]] <- matrix(0, nrow = num_areas, ncol = num_areas, 
                               dimnames = list(CB2013_regions, CB2013_regions))

  # make transition layer of the species age
  tr.layer = 1/paleo_map[[age]]
  ## adjust coordinates to focus the map on the indo-pacific territory
  x1 <- crop(tr.layer, extent(-180.5, 0, -90.5, 90.5))
  x2 <- crop(tr.layer, extent(0, 180.5, -90.5, 90.5))   
  extent(x1) <- c(180, 360.5, -90, 90)
  tr.layer_360 <- merge(x1, x2, tolerance = 0.2)
  names(tr.layer_360) = names(tr.layer)
  
  #tr.layer = paleo_map_Cao[[age]]
  tr.layer = gdistance::geoCorrection(gdistance::transition(tr.layer, function(x) mean(x), 8))
  tr.layer_360 = gdistance::geoCorrection(gdistance::transition(tr.layer_360, function(x) mean(x), 8))
  transLayerList[[age]] <- tr.layer
  
  for (i in 1:num_areas) {
    area1 <- CB2013_regions[i]
    coordsA1 <- na.exclude(reg_coords_age[reg_coords_age$CB2013 == area1, 1:2])

    # if using Cao 2007 data
    #coordsA1 <- na.exclude(pomacentridae_occ_dis_CB2013_regs[pomacentridae_occ_dis_CB2013_regs$CB2013 == area1,1:2])

    for (j in i:num_areas) {
      area2 <- CB2013_regions[j]
      coordsA2 <- na.exclude(reg_coords_age[reg_coords_age$CB2013 == area2, 1:2])
      # if using Cao 2007 data
      #coordsA2 <- na.exclude(pomacentridae_occ_dis_CB2013_regs[pomacentridae_occ_dis_CB2013_regs$CB2013 == area2,1:2])
      plot(raster(tr.layer), main = paste(area1, "-", area2))
      
      sp <- NULL
      # get shortest paths
      while(length(sp) < n_points) {
        sp <- c(sp, sapply(1:n_points, function(i) {
          xy1 <- as.numeric(coordsA1[sample(1:nrow(coordsA1), 1),1:2])
          xy2 <- as.numeric(coordsA2[sample(1:nrow(coordsA2), 1),1:2])
          if (xy1[1] > 100 & xy2[1] < (-50)) {
            xy2[1] <-  xy2[1] + 360
            tr360 = TRUE
          }
          else if (xy2[1] > 100 & xy1[1] < (-50)) {
            xy1[1] <-  xy1[1] + 360
            tr360 = TRUE
          } 
          else if (xy1[1] > 20 & xy2[1] < (-50)) {
            if (runif(1)>0.5){
              xy2[1] <-  xy2[1] + 360
              tr360 = TRUE
            }
          }
          else if (xy2[1] > 20 & xy1[1] < (-50)) {
            if (runif(1)>0.5){
              xy1[1] <-  xy1[1] + 360
              tr360 = TRUE
            }
          }
          else { tr360 = FALSE }
          
          if (tr360){
            tryCatch(gdistance::shortestPath(tr.layer_360, xy1, xy2, output="SpatialLines"),
                     error = function(e) NA)
          } else {
            tryCatch(gdistance::shortestPath(tr.layer, xy1, xy2, output="SpatialLines"),
                     error = function(e) NA)
          }
          })
          )

        sp <- sp[!sapply(sp, is.na)]
        if (length(sp) == 0){ break }
        sp <- sp[sapply(1:length(sp), function(v) nrow(coordinates(sp[[v]])[[1]][[1]]) > 1)]
        # shift lines that go over 180 degrees
        sp <- sapply(sp, function(s) if (extent(s)[2] > 180) { shift_line_longitudes(s) } else { s })
      }
      if (length(sp) > 0){ 
        sp <- sp[sample(1:length(sp), n_points)]
        # plot paths
        lapply(sp, lines)
        # get distances from each path
        dist <- unlist(lapply(sp, lengthLine)) / 1000
        # get the mean distance
        dist <- mean(dist, na.rm = TRUE)
      } else {
        dist <- Inf
      }
      cat(area1, "-", area2, ":", dist, "\n")
      
      dist_matrix[[age]][i, j] <- dist
      dist_matrix[[age]][j, i] <- dist_matrix[[age]][i, j]  # Distance matrix is symmetric
      
    }
  }
  # rename  dist matrix rownames and colnames to area codes for BioGeoBEARS 
  rownames(dist_matrix[[age]]) <- colnames(dist_matrix[[age]]) <- sapply(colnames(dist_matrix[[age]]), 
                                                                         function(i) names(area_code[area_code == i]))
}
#saveRDS(transLayerList, file = "./BioGeo/input/Cao2007_transitionLayers.rds")
#saveRDS(transLayerList, file = "./BioGeo/input/Descombes2018_transitionLayers.rds")
mindist <- min(unlist(dist_matrix)[unlist(dist_matrix) > 0], na.rm = T)

dist_matrix_rescaled <- lapply(dist_matrix, function(i) round(i/mindist, 2))

maxdist_rescaled <- max(unlist(dist_matrix_rescaled)[is.finite(unlist(dist_matrix_rescaled)) & unlist(dist_matrix_rescaled) > 0], na.rm = T)
for (i in 1:length(dist_matrix)) {
  dist_matrix_rescaled[[i]][!is.finite(dist_matrix_rescaled[[i]])] = 999
}
dist_matrix_rescaled <- dist_matrix_rescaled[names(paleo_map)][-1]

lapply(1:length(dist_matrix_rescaled), function(i) {
  write.table(dist_matrix_rescaled[[i]], file = "./BioGeo/input/damsel_biogeo_distmatrix_Descombes.txt", 
              quote = FALSE, row.names = FALSE, sep = "\t", append = ifelse(i == 1, FALSE, TRUE))
  cat("\n", file = "./BioGeo/input/damsel_biogeo_distmatrix_Descombes.txt", append = TRUE)})
cat("END\n", file = "./BioGeo/input/damsel_biogeo_distmatrix_Descombes.txt", append = TRUE)

######################
# get adjacency matrix
######################
adjacent_matrix <- list()
for (age in names(paleo_map)){
  
  adjacent_matrix[[age]] <- matrix(0, nrow = num_areas, ncol = num_areas, 
                            dimnames = list(LETTERS[1:num_areas], LETTERS[1:num_areas]))
  for (i in rownames(adjacent_matrix[[age]])){
    idx <- which(colnames(adjacent_matrix[[age]]) == i)
    adj_cells <- c(idx-1, idx, idx+1)
    adj_cells[adj_cells < 1] = num_areas
    adj_cells[adj_cells > num_areas] = 1
    
    dist_sorted <- sort(dist_matrix[[age]][i,])
    adj_cells <- adj_cells[dist_matrix[[age]][i, adj_cells] < 20000] # anything farther than 20,000 km is not adjacent
    adjacent_matrix[[age]][i,adj_cells] = 1
  }
}

# correct for the formation of the IAA area (not until 5-10 Mya)
for (i in names(adjacent_matrix)){
  if (as.numeric(sub("Ma", "", sub("X", "", i))) >= 10) {
    idxc <- which(adjacent_matrix[[i]][, names(area_code)[area_code == "IAA"]] == 1)

    adjacent_matrix[[i]][idxc, which(area_code == "IAA") + 1] = 1 # pass the adjacency area to neighbour
    adjacent_matrix[[i]][which(area_code == "IAA") + 1,idxc] = 1 # pass the adjacency area to neighbour
    
    adjacent_matrix[[i]][, which(area_code == "IAA") ] = 0 # remove IAA
    adjacent_matrix[[i]][which(area_code == "IAA") ,] = 0 # remove IAA
  } 
}
# correct for the closure and disappearance of the Tethys area (more than 23 Mya)
for (i in names(adjacent_matrix)){
  if (as.numeric(sub("Ma", "", sub("X", "", i))) <= 23) {
    idxc <- which(adjacent_matrix[[i]][, names(area_code)[area_code == "TS"]] == 1)
    idxr <- which(area_code == "TS") + 1
    if (idxr > num_areas) { idxr = 1 }
    
    adjacent_matrix[[i]][idxc,idxr] = 1 # pass the adjacency area to neighbour
    adjacent_matrix[[i]][idxr,idxc] = 1 # pass the adjacency area to neighbour
    
    adjacent_matrix[[i]][, which(area_code == "TS") ] = 0 # remove IAA
    adjacent_matrix[[i]][which(area_code == "TS") ,] = 0 # remove IAA
  } 
}
adjacent_matrix <- adjacent_matrix[names(paleo_map)][-1]

# save adjacent matrix of each time period
lapply(1:length(adjacent_matrix), function(i) {
  write.table(adjacent_matrix[[i]], file = "./BioGeo/input/damsel_biogeo_adjmatrix_Descombes.txt", 
              quote = FALSE, row.names = FALSE, sep = "\t", append = ifelse(i == 1, FALSE, TRUE))
  cat("\n", file = "./BioGeo/input/damsel_biogeo_adjmatrix_Descombes.txt", append = TRUE)})
cat("END\n", file = "./BioGeo/input/damsel_biogeo_adjmatrix_Descombes.txt", append = TRUE)

######################
# Allowed areas matrix
######################
allowed_areas <- adjacent_matrix
for (i in 1:length(allowed_areas)) {
  idx <- which(apply(allowed_areas[[i]], 1, sum) > 0)
  allowed_areas[[i]][idx,idx] = 1
}

# save allowed areas matrix of each time period
lapply(1:length(allowed_areas), function(i) {
  write.table(allowed_areas[[i]], file = "./BioGeo/input/damsel_biogeo_allwdmatrix_Descombes.txt", 
              quote = FALSE, row.names = FALSE, sep = "\t", append = ifelse(i == 1, FALSE, TRUE))
  cat("\n", file = "./BioGeo/input/damsel_biogeo_allwdmatrix_Descombes.txt", append = TRUE)})
cat("END\n", file = "./BioGeo/input/damsel_biogeo_allwdmatrix_Descombes.txt", append = TRUE)

######################
# remove some non-adjacent ranges
######################
geogfn <- "./BioGeo/input/CB2021_pomacentridae_geodata_input.txt"

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

class(tipranges@df)

# Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
max_range_size = max(apply(tipranges@df, 1, function(i) sum(as.numeric(i), na.rm = T)))
areas = getareas_from_tipranges_object(tipranges)
num_areas <- length(areas)
#areas = c("A", "B", "C", "D", "E", "F")

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# How many states/ranges, by default: 26
length(states_list_0based)

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based)) {    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) ) {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

# Look at the ranges list
ranges_list

# How many states/ranges, by default: 163
length(ranges_list)

nonadjacent <- list()
for (age in names(paleo_map)[-1]){
  
  nonadjacent[[age]] <- sapply(ranges_list, function(x) {
                        adj <- TRUE
                        if (nchar(x) > 1){
                          splt <- strsplit(x, "")[[1]]
                          if (first(splt) == colnames(adjacent_matrix[[age]])[1] &&
                              last(splt) == colnames(adjacent_matrix[[age]])[ncol(adjacent_matrix[[age]])]) {
                            if (splt[2] == colnames(adjacent_matrix[[age]])[2]) {
                              splt <- c(last(splt), splt[-length(splt)])
                            } else {
                              splt <- c(splt[-1], first(splt))
                            }
                          }
                          for (i in 1:(length(splt)-1)){
                            j = i+1
                            adj <- adjacent_matrix[[age]][splt[i],splt[j]] == 1
                            if (!adj){
                              break
                            }
                          }
                        }
                        !adj
                      })
}

ranges_list = sapply(1:length(nonadjacent), function(i) ranges_list[!nonadjacent[[i]]])
sapply(ranges_list, length)  # now 14 14 13 13 16 
unique_states <- unique(unlist(ranges_list))
# save the UNIQUE states list into the BioGeoBEARS_run_object
saveRDS(unique_states, "./BioGeo/input/CB2013_biogeo_unique_states_Descombes.rds")

new_states_list = sapply(1:length(nonadjacent), function(i) states_list_0based[!nonadjacent[[i]]])
sapply(new_states_list, length)  # now 14 14 13 13 16 
# save the NEW states list into the BioGeoBEARS_run_object
saveRDS(new_states_list, "./BioGeo/input/CB2013_biogeo_selected_states_Descombes.rds")


