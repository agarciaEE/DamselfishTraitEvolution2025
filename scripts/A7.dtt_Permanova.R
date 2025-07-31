#####################################################################
############  Comparison PC scores vs ecological traits ##############
######################################################################
library(tidyverse)
library(vegan)      # for adonis2
library(ape)        # phylogenetics
library(phytools)
library(dispRity)
library(geomorph)
library(adephylo)
library(geiger)
library(caper)
library(dplyr)

## set workfolder
path = "~/Unil/Research/5_Damselfish_evo/"
setwd(path) ###########

# load input data and custom functions
source("./Rscripts/A0.Input_data.R")
tree_species <- tree$tip.label

# PERMANOVA to test variance partitioning
##############################################
adonis2_color <- adonis2(
  as.matrix(color_indPCA[, 1:8]) ~ species,
  data = color_indPCA,
  method = "euclidean",
  permutations = 999
)
print(adonis2_color)

adonis2_morpho <- adonis2(
  as.matrix(morpho_data[, 1:(ncol(morpho_data)-1)]) ~ species,
  data = morpho_data,
  method = "euclidean",
  permutations = 999
)

# Join species-level color data with ecology
color_eco <- data.frame(color_species_mean[tree_species,], eco_traits_subset[tree_species, c("DietEcotype", "Farming", "habitat1")])
colnames(color_eco)[1:8] <- paste0("color_", colnames(color_eco)[1:8])
colnames(color_eco)[9:11] <- c("diet", "symbiosis", "habitat")

adonis2_eco_color <- adonis2(
  as.matrix(dplyr::select(color_eco, starts_with("color"))) ~ diet + habitat + symbiosis,
  data = color_eco,
  method = "euclidean",
  permutations = 999
)
print(adonis2_eco_color)

# Join species-level morpho data with ecology
morpho_eco <- data.frame(morpho_species_mean[tree_species,], eco_traits_subset[tree_species, c("DietEcotype", "Farming", "habitat1")])
colnames(morpho_eco)[1:11] <- paste0("morpho_", colnames(morpho_eco)[1:11])
colnames(morpho_eco)[12:14] <- c("diet", "symbiosis", "habitat")

adonis2_eco_morpho <- adonis2(
  as.matrix(dplyr::select(morpho_eco, starts_with("morpho"))) ~ diet + habitat + symbiosis,
  data = morpho_eco,
  method = "euclidean",
  permutations = 999
)
print(adonis2_eco_morpho)

# phylogenetic disparity-through-time (DTT) 
##############################################
color_species_mean <- color_indPCA %>%
  dplyr::group_by(species) %>%
  summarise_all(mean, rm = TRUE) %>%
  column_to_rownames("species") %>%
  as.matrix()

# check data
color_species_mean <- color_species_mean[tree_species, ]

sum(is.na(color_species_mean))
all(rownames(color_species_mean) %in% tree_species)

# using traits
morpho_species_mean <- morpho_data %>%
  group_by(species) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("species") %>%
  as.matrix()

# using PCs 
morpho_species_mean <- morphoPCs[tree_species, ]

# check data
sum(is.na(morpho_species_mean))
all(rownames(morpho_species_mean) %in% tree_species)

# Calculate disparity-through-time
pdf("./new_figures/color_dtt_plots.pdf", height = 17, width = 12)
par(mfrow = c(4,2))
for (pc in 1:8){
  dtt_color <- dtt(phy = tree,
                   data = scale(color_species_mean[tree_species,pc]),
                   index = "avg.sq", 
                   nsim=100, plot=FALSE)
  dtt_color$MDIpVal <- compute.MDIpval(dtt_color)
  par(mar = c(7,7,4,7))
  plot.dtt(dtt_color, phy = tree, col = "grey20", ylab = "Relative disparity of coloration patterns", 
           cex = 1.5, cex.axis = 1.2, cex.lab = 1.5, main = paste("color PC", pc), cex.main = 2)
  
}
dev.off()

pdf("./new_figures/morpho_dtt_plots.pdf", height = 17, width = 12)
par(mfrow = c(4,2))
for (pc in 1:8){
dtt_morpho <- dtt(phy = tree,
                  data = morpho_species_mean[tree_species,pc],
                  index = "avg.sq", 
                  nsim=100, plot=FALSE)
dtt_morpho$MDIpVal <- compute.MDIpval(dtt_morpho)
par(mar = c(7,7,4,7))
plot.dtt(dtt_morpho, phy = tree, col = "grey20", ylab = "Relative disparity of morphological patterns", 
         cex = 1.5, cex.axis = 1.2, cex.lab = 1.5, main = paste("morpho PC", pc), cex.main = 2)

}
dev.off()

# Phylogenetic signal
##############################################
## color
image_data <- extractImageData(imgPCA_1183damsel$images, impute = TRUE)
image_species_data <- data.frame(species = color_indPCA$species, image_data)

# average by species
img_species_means <- image_species_data %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(across(everything(), mean)) %>%
  ungroup()

p <- ncol(image_data) 
n_species <- length(tree_species)  
img_species_means <- as.matrix(img_species_means[, -1])  

# Reshape this matrix into a 3D array (pixels × channels × species)
# thre RGB dimensions and each "landmark" is a single pixel in the image.
n_pixels <- ncol(img_species_means) / 3
img_species_means_array <- array(
  img_species_means, 
  dim = c(n_pixels, 3, n_species)  # n_pixels pixels (landmarks), 3 channels (R, G, B), n_species species
)

color_phylo_K <- physignal(A = img_species_means_array, phy = tree, iter = 999)
print(color_phylo_K)

## morpho
# Aggregate individual gpa coordinates by species
p <- dim(morpho_gpa$coords)[1]  # number of landmarks
k <- dim(morpho_gpa$coords)[2]  # number of dimensions
n_species <- length(tree_species)
species_mean_coords <- array(NA, dim = c(p, k, n_species), dimnames = list(NULL, NULL, tree_species))

# Compute species means
for (i in seq_along(tree_species)) {
  sp <- tree_species[i]
  sp_indices <- which(morpho_data$species == sp)
  species_mean_coords[,,i] <- apply(morpho_gpa$coords[,,sp_indices, drop = FALSE], c(1, 2), mean)
}

morpho_phylo_K <- physignal(A = species_mean_coords, phy = tree, iter = 999)
print(morpho_phylo_K)

# Compute distance matrices
# using traits
morpho_species_mean <- morpho_data %>%
  group_by(species) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("species") %>%
  as.matrix()

color_dist <- dist(color_species_mean)
morpho_dist <- dist(morpho_species_mean)
phylo_dist <- cophenetic(tree)

color_mantel_test <- mantel(morpho_dist, phylo_dist, method = "pearson", permutations = 999)
print(mantel_test)

morpho_mantel_test <- mantel(color_dist, phylo_dist, method = "pearson", permutations = 999)
print(mantel_test)

