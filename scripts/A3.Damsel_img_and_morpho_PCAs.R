version
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
#load ColorAR package
#devtools::install_github("agarciaEE/ColorAR", auth_token = "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")
#install.packages("~/Unil/Research/Rpackages/ColAR", type = "source", repos = NULL)
library(ColorAR)
library(geomorph)

## set workfolder
path = "~/Unil/Research/5_Damselfish_evo/"
setwd(path) ###########

source("./Rscripts/custom_functions.R")

damsel_tree_subset <- read.tree("./data/damsel_subset.tre")

img_filenames <- file.path("images/TransIntImgPNGColorCorrected/", 
                           list.files("images/TransIntImgPNGColorCorrected/"))
img_names <- sub("\\.[a-z]*$", "", basename(img_filenames))

# number of images
length(img_names)

imgTransInt.list_res300 <- makeList(img_names, 'image', prepath = "images/TransIntImgPNGColorCorrected", extension = ".png")

# Remove 4 channel (alpha) if exists
if (!all(sapply(imgTransInt.list_res300, nlayers) == 3)) {
  imgTransInt.list_res300 <- sapply(imgTransInt.list_res300, function(i) i[[1:3]])
}

# flip y coordinates
imgTransInt.list_res300 <- sapply(imgTransInt.list_res300, raster::flip, "y")

# mask background
Rmask <- raster::stack("./images/mask2.png")[[4]]
Rmask[Rmask == 0] = NA
Rmask[Rmask > 0] = 1
Rmask <- raster::flip(Rmask, "y")

imgTransInt.list_res300 <- sapply(imgTransInt.list_res300, raster::mask, Rmask)

# load("./rdata/img1183_PCA1-8_whole.RData")

# Get species names from image ID
species_names <- sapply(img_names, function(i) {
  # Extract species name from the first part of the string (before the first underscore)
  name <- strsplit(i, "_")[[1]][1]
  
  # Capitalize the first letter and replace any hyphen with an underscore
  name <- paste0(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name)))
  name <- sub("-", "_", name)
  
  return(name)
})

# get unique species
species <- unique(species_names)

# get species in tree
tree_species <- damsel_tree_subset$tip.label

# check names in case misspeling or mismatch between image and tree species
check.names <- as.data.frame(t(sapply(species[!species %in% tree_species], function(i) check_spelling(i, tree_species))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree

# Correct misspelled species names
species_names[species_names == "Pycnochromis_lineata"] = "Pycnochromis_lineatus"
species_names[species_names == "Chrysiptera_maurinae"] = "Chrysiptera_maurineae"
species_names[species_names == "Plectroglyphidodon_leucozona"] = "Plectroglyphidodon_leucozonus"
species_names[species_names == "Pomacentrus_caeruleopuncatus"] = "Pomacentrus_caeruleopunctatus"
species_names[species_names == "Stegastes_flaviatus"] = "Stegastes_flavilatus"

# Correct image name IDs
img_names <- paste0(species_names, "-", sapply(strsplit(img_names, "_"), function(i) i[[2]]))
names(imgTransInt.list_res300) <- img_names
names(species_names) <- img_names

# obtain again unique species
species <- unique(species_names)

# check missing species 
species[!species %in% tree_species] # 13 missing species in the tree
tree_species[!tree_species %in% species] # 0 missing species from image data

# Reorder species based on the tree 
species <- species[match(tree_species, species)]
n_sp <- length(species)

# Assign color to species
genus <- sapply(species, function(i) strsplit(i, "_")[[1]][1])
names(genus) <- species
genus_cols_legend <- setNames(colorRampPalette(viridis(9))(n=length(unique(genus))),
                              unique(genus))
genus_cols <- setNames(rep(genus_cols_legend)[as.factor(genus)], names(genus))

saveRDS(imgTransInt.list_res300, "./rdata/imgTransInt.list_res300.rds")

######################################################################
#############                 Run image PCA             ##############
######################################################################
imgPCA_1183damsel <- imagePCA(imgTransInt.list_res300, tree = NULL, pcs = 1:8, impute.NA = TRUE, 
                              res = NULL, interpolate = FALSE, scale = FALSE, plot = FALSE) 

saveRDS(imgPCA_1183damsel, "./rdata/imgPCA_1183damsel.RDS")

imgPCA_1183damsel <- readRDS("./rdata/imgPCA_1183damsel.RDS")
imgTransInt.list_res300 <- readRDS("./rdata/imgTransInt.list_res300.rds")

###################################
##### species mean images
###################################
sp.mean_imgList <- list()
for (i in species){
  
  idx <- which(species_names == as.character(i))
  imgs <- imgTransInt.list_res300[idx]
  RR <- mean(stack(sapply(imgs, function(x) x[[1]]), na.rm =  T))
  GG <- mean(stack(sapply(imgs, function(x) x[[2]]), na.rm =  T))
  BB <- mean(stack(sapply(imgs, function(x) x[[3]]), na.rm =  T))
  sp.mean_imgList[[i]] <- stack(RR, GG, BB)
  names(sp.mean_imgList[[i]]) <- c("R", "G", "B")
}

###################################
##### plot image PCA
###################################
for (pc in c(1, 3, 5, 7)) {
  
  pdf(paste0("./1183imgs_PC", pc, pc +1, ".pdf"), width = 14, height = 11)
  plot.imgPCA(imgPCA_1183damsel, pcs = pc:(pc+1), group_vector = species_names[names(imgPCA_1183damsel$images)], group_imgList = sp.mean_imgList,
              show.sd = FALSE, display =  "images", pixel.contributions = TRUE,
              color.contributions = NULL, img.size = 0.1)
  dev.off()
}

###################################
##### Run image PCA on species mean images
###################################

imgPCA_1183damsel_spmean <- imagePCA(sp.mean_imgList, tree = NULL, pcs = 1:8, impute.NA = TRUE, 
                                     res = NULL, interpolate = FALSE, scale = FALSE, plot = FALSE) 

#############################################################nulldefault()#########
##### plot image PCA
###################################
for (pc in c(1, 3, 5, 7)) {
  
  pdf(paste0("./1183imgs_PC", pc, pc +1, "_spmean_scaled.pdf"), width = 14, height = 11)
  plot.imgPCA(imgPCA_1183damsel_spmean, pcs = pc:(pc+1), 
              show.sd = FALSE, display =  "images", pixel.contributions = TRUE,
              color.contributions = c("#F4A6A6", "#A6C8F4"), img.size = 0.1)
  dev.off()
}

###################################
##### Average PC scores by species
###################################

pcdata = as.data.frame(imgPCA_1183damsel$pca$x[,1:8])
pcdata$species <- species_names[names(imgPCA_1183damsel$images)]

# remove species not in tree
pcdata <- pcdata[pcdata$species %in% species,]

species_pcdata <- pcdata %>%
  group_by(species) %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))

write.csv(species_pcdata, "./data/1188img_colorPCA_spmean_scaled.csv", row.names = F)


######################################################################
############                 Run morpho PCA             ##############
######################################################################

imglist_nobg <- readRDS("./rdata/imglist_nobg.rds")

# load traits
eco_traits_subset <- read.csv("./data/eco_traits_subset_final.csv", row.names = 1)
morpho_gpa <- readRDS("./rdata/1166morho_gpagen.RDS")
morpho_data <- read.csv("./data/1166damsel_morpho.csv", row.names = 1)

sp_names <- setNames(morpho_data$species, rownames(morpho_data))

# Extract continuous trait columns
is_continuous <- sapply(morpho_data, is.numeric)
morpho_traits <- colnames(morpho_data)[is_continuous]

# Remove depth-related traits from morphology analysis
sp_morpho <- eco_traits_subset[tree_species,morpho_traits]

# Run phylogenetic PCA
morphoPCA_scaled <- phytools::phyl.pca(damsel_tree_subset, apply(sp_morpho, 2, scale), method = "lambda", mode = "corr")

# Convert to prcomp-like object (optional, if you need consistency with other workflows)
morphoPCA_scaled <- phytools::as.prcomp(morphoPCA_scaled)

saveRDS(morphoPCA_scaled, "./rdata/morpho_phyPCA_scaled.rds")

# Extract PC scores
morphoPCs <- morphoPCA$x

for (pc in c(1, 3, 5, 7)) {
  pcx <- pc
  pcy <- pc+1
  
  pdf(paste0("./morphoPCA_", pcx, pcy, "no_allardi.pdf"), width = 14, height = 11)
  plot_morphoPCA(morphoPCA, morpho_gpa, sp_names, genus_cols, imglist_nobg,
                 pcx = pcx, pcy = pcy, display = "images", cex.labels = 1.5, r = 0.5) 
  dev.off()
  
}

#