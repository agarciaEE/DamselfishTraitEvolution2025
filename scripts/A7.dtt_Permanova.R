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
library(phyloch)

data(strat2012)

## set workfolder
path = "~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/"
setwd(path) ###########

# load input data and custom functions
source("./scripts/A1.Input_data.R")
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
print(adonis2_morpho)

color_ind_names <- rownames(color_indPCA)
morpho_ind_names <- sub(" ", "_", sub("_", "-", gsub("-", " ", rownames(morpho_data))))

ind_names <- intersect(color_ind_names, morpho_ind_names)

morpho_data <- morpho_data[!duplicated(morpho_ind_names),]
morpho_ind_names <- sub(" ", "_", sub("_", "-", gsub("-", " ", rownames(morpho_data))))
sum(duplicated(morpho_ind_names))
rownames(morpho_data) = morpho_ind_names

ind_names <- intersect(color_ind_names, morpho_ind_names)
color_morpho_ind_data <- cbind(color_indPCA[ind_names, 1:8], morpho_data[ind_names, 1:11])

adonis2_color_morpho <- adonis2(
  as.matrx(color_morpho_ind_data) ~ species,
  data = color_indPCA[ind_names,],
  method = "euclidean",
  permutations = 999
)
print(adonis2_color_morpho)

# Join species-level color data with ecology
color_species_mean <- color_indPCA %>%
  dplyr::group_by(species) %>%
  summarise_all(mean, rm = TRUE) %>%
  column_to_rownames("species") %>%
  as.matrix()

# check data
color_species_mean <- color_species_mean[tree_species, ]

sum(is.na(color_species_mean))
all(rownames(color_species_mean) %in% tree_species)

color_eco <- data.frame(color_species_mean[tree_species,], eco_traits_subset[tree_species, c("DietEcotype", "Symbiosis", "Habitat")])
colnames(color_eco)[1:8] <- paste0("color_", colnames(color_eco)[1:8])
colnames(color_eco)[9:11] <- c("diet", "symbiosis", "habitat")

adonis2_eco_color <- adonis2(
  as.matrix(dplyr::select(color_eco, starts_with("color"))) ~ diet + habitat + symbiosis,
  data = color_eco,
  method = "euclidean",
  permutations = 999, by = "margin"
)
print(adonis2_eco_color)

# Join species-level morpho data with ecology
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

morpho_eco <- data.frame(morpho_species_mean[tree_species,], eco_traits_subset[tree_species, c("DietEcotype", "Symbiosis", "Habitat")])
colnames(morpho_eco)[1:11] <- paste0("morpho_", colnames(morpho_eco)[1:11])
colnames(morpho_eco)[12:14] <- c("diet", "symbiosis", "habitat")

adonis2_eco_morpho <- adonis2(
  as.matrix(dplyr::select(morpho_eco, starts_with("morpho"))) ~ diet + habitat + symbiosis,
  data = morpho_eco,
  method = "euclidean",
  permutations = 999, by = "margin"
)
print(adonis2_eco_morpho)

color_morpho_species_mean <- color_morpho_ind_data %>%
  mutate(species = sapply(ind_names, \(x) strsplit(x, "-")[[1]][1])) %>%
  group_by(species) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  filter(species %in% tree_species) %>%
  column_to_rownames("species") %>%
  as.matrix() %>%
  scale() %>%
  na.exclude() %>%
  as.data.frame() 

# check data
sum(is.na(color_morpho_species_mean))
all(rownames(color_morpho_species_mean) %in% tree_species)

color_morpho_eco <- color_morpho_species_mean %>%
  rownames_to_column("species") %>%
  left_join(eco_traits_subset %>% 
              dplyr::select(DietEcotype, Symbiosis, Habitat) %>%
              rename(diet = DietEcotype, symbiosis = Symbiosis, habitat = Habitat) %>%
              rownames_to_column("species"), by = "species") %>%
  column_to_rownames("species")
  
colnames(color_morpho_eco)[1:8] <- paste0("color", colnames(color_morpho_eco)[1:8])

adonis2_eco_color_morpho <- adonis2(
  as.matrix(color_morpho_eco %>%
              dplyr::select(-c(diet, habitat, symbiosis))) ~ diet + habitat + symbiosis,
  data = color_morpho_eco,
  method = "euclidean",
  permutations = 999, by = "margin"
)
print(adonis2_eco_color_morpho)

save(adonis2_color, adonis2_morpho, adonis2_color_morpho, adonis2_eco_color, adonis2_eco_morpho, adonis2_eco_color_morpho, file = "./rdata/permanova_tests.rda")

# phylogenetic disparity-through-time (DTT) 
##############################################
dtt_colorList <- list()
# Calculate disparity-through-time
for (pc in 1:8){
  dtt_color <- dtt(phy = tree,
                   data = scale(color_species_mean[tree_species,pc]),
                   index = "avg.sq", 
                   nsim=1000, plot=FALSE)
  dtt_color$MDIpVal <- compute.MDIpval(dtt_color)
  dtt_colorList[[paste0("colorPC", pc)]] <- dtt_color
}

dtt_morphoList <- list()
for (pc in 1:8){
  dtt_morpho <- dtt(phy = tree,
                    data = morpho_species_mean[tree_species,pc],
                    index = "avg.sq", 
                    nsim=1000, plot=FALSE)
  dtt_morpho$MDIpVal <- compute.MDIpval(dtt_morpho)
  dtt_morphoList[[paste0("morphoPC", pc)]] <- dtt_morpho
}

pdf("./figures/DTT/color_dtt_plots.pdf", height = 17, width = 12)
par(mfrow = c(4,2))
lapply(1:length(dtt_colorList), function(pc) {
  if (pc > 4) {
    par(mar = c(4.5,5.5,4,4.5))
  } else {
    par(mar = c(4.5,4.5,4,5.5))
  }
  plot.dtt(dtt_colorList[[pc]], phy = tree, col = "grey20", ylab = "Relative disparity", cex.main = 2.5,
           cex = 2, cex.axis = 1.5, cex.lab = 2, main = paste("color PC", pc), cex.main = 2)
})
dev.off()

pdf("./figures/DTT/morpho_dtt_plots.pdf", height = 17, width = 12)
par(mfrow = c(4,2))
lapply(1:length(dtt_morphoList), function(pc) {
  if (pc > 4) {
    par(mar = c(4.5,5.5,4,4.5))
  } else {
    par(mar = c(4.5,4.5,4,5.5))
  }
  plot.dtt(dtt_morphoList[[pc]], phy = tree, col = "grey20", ylab = "Relative disparity", cex.main = 2.5,
           cex = 2, cex.axis = 1.5, cex.lab = 2, main = paste("morpho PC", pc), cex.main = 2)
})
dev.off()

cols <- c(
  "#332288", "#117733", "#44AA99", "#88CCEE",
  "#DDCC77", "#CC6677", "#AA4499", "#882255"
)

pdf("./figures/DTT/multidtt_plot.pdf", height = 15, width = 12)
par(mfrow = c(2,1), mar = c(5, 6, 4, 5) + 0.1)
# --- (a) Color DTT ---
plot.multidtt(
  dtt_colorList,
  phy = tree,
  col = cols,
  ylab = "Relative disparity",
  lwd = 3,
  cex.axis = 1.25,
  cex.legend = 1.3,
  cex.lab = 1.5,
  main = "Color disparity-through-time",
  cex.main = 1.5,
  time_axis = "absolute",
  xlab = "Time (Mya)",
  xaxt = "n"
)
axisGeo(GTS = strat2012, tip.time = 0, unit ="epoch", ages = 
          TRUE, cex = 1.2, col = RColorBrewer::brewer.pal(9, "YlOrBr")[1:5], 
        texcol = "black", offset = -0.22,
        gridcol = "grey60")
mtext("(a)", side = 3, adj = -0.1, padj = -0.3, line = 1.2, font = 2, cex = 1.5)
# --- (b) Morphological DTT ---
plot.multidtt(
  dtt_morphoList,
  phy = tree,
  col = cols,
  ylab = "Relative disparity",
  lwd = 3,
  cex.axis = 1.25,
  cex.legend = 1.3,
  cex.lab = 1.5,
  main = "Morphological disparity-through-time",
  cex.main = 1.5,
  time_axis = "absolute",
  xlab = "Time (Mya)",
  xaxt = "n"
)
axisGeo(GTS = strat2012, tip.time = 0, unit ="epoch", ages = 
          TRUE, cex = 1.2, col = RColorBrewer::brewer.pal(9, "YlOrBr")[1:5], 
        texcol = "black", offset = -0.45,
        gridcol = "grey60")
mtext("(b)", side = 3, adj = -0.1, padj = -0.3, line = 1.2, font = 2, cex = 1.5)
dev.off()

save(dtt_colorList, dtt_morphoList, file = "./rdata/dtt_analyses.rda")
load("./rdata/dtt_analyses.rda")

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
dimnames(img_species_means_array) <- list(1:n_pixels, c("R", "G", "B"), tree_species)

color_phylo_K <- geomorph::physignal(A = img_species_means_array, phy = tree, iter = 999)
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

morpho_phylo_K <- geomorph::physignal(A = species_mean_coords, phy = tree, iter = 999)
print(morpho_phylo_K)

save(color_phylo_K, morpho_phylo_K, file = "./rdata/phylo_signal.rda")

# Compute distance matrices
# using traits
color_dist <- dist(color_species_mean)
morpho_dist <- dist(morpho_species_mean)
phylo_dist <- cophenetic(tree)

color_mantel_test <- mantel(color_dist, phylo_dist, method = "pearson", permutations = 999)
print(color_mantel_test)

morpho_mantel_test <- mantel(morpho_dist, phylo_dist, method = "pearson", permutations = 999)
print(morpho_mantel_test)

