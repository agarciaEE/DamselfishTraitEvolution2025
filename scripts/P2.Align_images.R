####################
##### COLOR ANALYSES USING Damselfishes webdata pictures
###
# make sure these dependencies are installed:
###
library(patternize)
library(factoextra)
library(rasterVis)
library(vegan)
library(reshape)
library(rasterVis)
library(raster)
library(devtools)
library(Morpho)
library(sp)
library(spatialEco)
###
# The viridis package provides colour blind friendly colour schemes
###
library(viridis)
library(colorspace)
## Additional and handy packages
library(plyr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

## comparative phylogenetic methods packages
library(ape)
library(phytools)
library(geiger)
library(phangorn)
library(OUwie)
library(nlme)
library(BAMMtools)
library(picante)
library(coda)
library(caper)
require(mvtnorm)
library(rfishbase)
library(tidyverse)
library(taxize)
library(fasterize)
#load ColorAR package
#devtools::install_github("agarciaEE/ColorAR", auth_token = "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")
#install.packages("~/Unil/Research/Rpackages/ColAR", type = "source", repos = NULL)
library(ColorAR)

## set workfolder
path = "~/Unil/Research/5_Damselfish_evo/"
setwd(path) ###########

# load custom functions
source("./Rscripts/custom_functions.R")

################################################################################
##### -------------------------------------------------------------------- #####
#####                       Load image and landmarks                       #####
##### -------------------------------------------------------------------- #####
################################################################################
### NOTE:
### Images have been previously processed on Adobe Photoshop to:
#     - Flip picture to be left body sided
#     - Color correction
#     - W/B Balance
#     - Set landmarks on processed pictures

## Load images filenames
img_filenames <- file.path("images/all_selected_images", list.files("images/all_selected_images/"))
img_names <- sub("\\.[a-z]*$", "", basename(img_filenames))
length(img_names)

## load landmarks filenames
lmk_filenames <- file.path("images/all_selected_landmarks", list.files("images/all_selected_landmarks/"))
lmk_names <- sub("\\.[a-z]*$", "", basename(lmk_filenames))
length(lmk_names)

# check match between landmarks and images
img_names[!img_names %in% lmk_names]
lmk_names[!lmk_names %in% img_names]

# subset images with landmarks
img_filenames <- img_filenames[img_names %in% lmk_names]
img_names <- img_names[img_names %in% lmk_names]

## Load images
imageList = makeList(img_filenames, 'image', prepath = NULL, extension = "")
names(imageList) <- img_names

## load landmarks
landmarkList = makeList(lmk_filenames, 'landmark', prepath = NULL, extension = "")
names(landmarkList) <- lmk_names

# order landmarkList to match imageList
landmarkList <- landmarkList[names(imageList)]

# check number of landmarks are the correct (46)
if (any(sapply(landmarkList, nrow) != 46)) {
  message("Following landmarked images do not have exactly 46 landmark points and will be removed:",
      paste("\n\t- ", 
            names(landmarkList)[sapply(landmarkList, nrow) != 46]))
  wrong_lndmks <- names(landmarkList)[sapply(landmarkList, nrow) != 46]
  landmarkList <- landmarkList[sapply(landmarkList, nrow) == 46]
  imageList <- imageList[names(landmarkList)]
  img_names <- names(imageList)
  lmk_names <- names(landmarkList)
}

# set internal landmarks to drop when building the mask
lndmks_2drop <- c(1:10, 45, 46)

# set store raster on memory to avoid issues
rasterOptions(todisk = F)

# get species names and genera
sp_names <- sub("-", " ", img_names)
sp_names <- sapply(sp_names, function(i) strsplit(i, "_")[[1]][1])
table(sp_names)

sp_genera <- sapply(sp_names, function(i) strsplit(i, " ")[[1]][1])
table(sp_genera)

################################################################################
##### -------------------------------------------------------------------- #####
#####      Transform images using landmark based Procrustes analyses       #####
##### -------------------------------------------------------------------- #####
################################################################################
####

# get mean shape and make mask
lanArray <- patternize::lanArray(landmarkList, adjustCoords = T,  imageList)
invisible(utils::capture.output(transformed <- Morpho::procSym(lanArray)))
refShape <- transformed$mshape
# write.table(refShape, "./images/mean_shape_reference.txt", quote = F, row.names = F, col.names = F)
refShape <- as.matrix(read.table("./images/mean_shape_reference.txt"))
sp.ref = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(refShape[-lndmks_2drop,])), 1)))
sp.ref <- spatialEco::rotate.polygon(sp.ref, angle=180)
sp.ref <- smoothr::smooth(sp.ref, method = "ksmooth", smooth = 1)
sp.ref$fish = 1
# sf::st_write(sf::st_as_sf(sp.ref), "./images/mean_shape_mask.geojson")
plot(sp.ref)

# get aspect ration of images
asp <- sapply(imageList, function(i) extent(i)[2]/extent(i)[4])

# transform images
imageTransformation(sampleList = imageList, landList = landmarkList, adjustCoords = TRUE, transformRef = refShape,
                         drop = lndmks_2drop,
                         crop = FALSE, cropOffset = c(0, 0, 0, 0), res = 300, keep.ASP  = T, aspect.ratio = mean(asp),
                         removebg.by = "landmarks", smooth = 1, rescale = FALSE,
                         transformType = "tps", focal = F, sigma = 3, interpolate =  5,
                         plot = FALSE, save = TRUE, dir = "images/new_Transformed_images", overwrite = FALSE)
plotRGB(stack("./images/new_Transformed_images/Abudefduf septemfasciatus_5.tif"))

img_filenames <- file.path("images/new_Transformed_images/", list.files("images/new_Transformed_images/"))
img_names <- sub("\\.[a-z]*$", "", basename(img_filenames))
length(img_names)

imgTrans.list_res300 <- makeList(img_names, 'image', prepath = "images/new_Transformed_images", extension = ".tif")

mask <- raster::raster("images/damsel_mask.tif")
crs(mask) = NA

# write as PNG
for (i in names(imgTrans.list_res300)) {
  img <- raster:mask(imgTransInt.list_res300[[i]], mask)
  raster::writeRaster(img, filename = file.path("images/new_Transformed_images", 
                                                           sub(".png", ".tif", i)), overwrite = TRUE)
  
  png(filename = paste0("images/TransIntImgPNG/", i, ".png"), width = 467, height = 300, bg = "transparent")
  plotRGB(img)
  dev.off()
}
