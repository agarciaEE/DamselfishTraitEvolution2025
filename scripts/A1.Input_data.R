#=============================#
# Libraries
#=============================#
lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", 
         "lme4", "mgcv", "phytools", "corHMM",
         "lmtest", "ColorAR", "castor", "tibble",
         "RPANDA","OUwie", "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

#=============================#
# Settings
#=============================#
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

# load custom functions
source("./scripts/custom_functions.R")

col.palette <- sample(unique(as.character(sapply(c("Zissou1", "Darjeeling1", "Darjeeling2", "FantasticFox1"), function(x) wes_palette(x, 5)))))
ncores <- 5
color_analysis = TRUE
morpho_analysis = TRUE

#=============================#
# Data Preparation
#=============================#

# tree
tree <- read.tree("./data/damsel_subset.tre")

tree$root.edge = NULL
nEdges <- Nedge(tree) # total number of edges
nTips <- Ntip(tree)

# load traits
eco_traits_subset <- read.csv("./data/eco_traits_subset_final.csv", row.names = 1)

# discrete traits as factors
eco_traits_subset$DietEcotype <- factor(eco_traits_subset$DietEcotype, levels = c("B", "I", "P"))
eco_traits_subset$Symbiosis <- factor(eco_traits_subset$Symbiosis, levels = c("No", "Farming", "Mutualism"))
eco_traits_subset$Habitat <- factor(eco_traits_subset$Habitat, levels = c("sea anemone", "freshwater", "seagrass-rubble-sand", 
                                                                            "rocky-reef", "coral-reef"))

levels(eco_traits_subset$DietEcotype) <- c("Benthic", "Intermediate", "Pelagic")
levels(eco_traits_subset$Symbiosis ) <- c("Free-living","Comensalistic","Mutualistic")
levels(eco_traits_subset$Habitat ) <- c("sea anemone", "freshwater", "non-reef", 
                                         "rocky-reef", "coral-reef")

discrete_traits <- c("DietEcotype", "Symbiosis", "Habitat")

# Biogeographic Reconstruction
bsm_pomacentridae <- readRDS("./Rdata/BSM_pomacentridae_DEC_J_CB2013_regs_strat_allwdareas_shortestdistmat_Descombes.rds")

# color PCA
imgPCA_1183damsel <- readRDS("./rdata/imgPCA_1183damsel_scaled.RDS") # image PCA
color_indPCA <- as.data.frame(imgPCA_1183damsel$pca$x[,1:8])
rownames(color_indPCA) <- names(imgPCA_1183damsel$images)
color_indPCA$species <- sapply(rownames(color_indPCA), function(i) strsplit(i, "-")[[1]][1])

colorPCs <- color_indPCA %>%
  dplyr::group_by(species) %>%
  dplyr::summarise_all(mean, rm = TRUE) %>%
  tibble::column_to_rownames("species") %>%
  as.matrix()

if (!all(tree$tip.label %in% rownames(colorPCs))) {
  cat("Not all species in tree are present in color PCA. Check possible misspelling:\n")
  check.names <- as.data.frame(t(sapply(rownames(colorPCs), function(i) check_spelling(i, tree$tip.label))))
  mispelled_tree <- check.names[check.names$score != 1,]
  print(mispelled_tree)
}
if (!all(rownames(colorPCs) %in% tree$tip.label)) {
  cat("Some species in color PCA are not present in the tree:\n")
  print(rownames(colorPCs)[!rownames(colorPCs) %in% tree$tip.label])
  colorPCs <- colorPCs[tree$tip.label,]
}

# independent PCs
for (i in 1:ncol(colorPCs)){
  assign(paste0("colorPC", i), setNames(colorPCs[,i], rownames(colorPCs)))
}

## Morpho PCA data 
morpho_gpa <- readRDS("./rdata/1166morho_gpagen.RDS") # Procrustes landmark coordinates
morpho_data <- read.csv("./data/1166damsel_morpho.csv", row.names = 1) # Morphometric dataset
morphoPCA <- readRDS("./rdata/morpho_phyPCA.rds") # morpho phylogenetic PCA
morphoPCs <- morphoPCA$x

if (!all(tree$tip.label %in% rownames(morphoPCs))) {
  cat("Not all species in tree are present in color PCA. Check possible misspelling:")
  check.names <- as.data.frame(t(sapply(rownames(morphoPCs), function(i) check_spelling(i, tree$tip.label))))
  mispelled_tree <- check.names[check.names$score != 1,]
  mispelled_tree
}
if (!all(rownames(morphoPCs) %in% tree$tip.label)) {
  cat("Some species in morpho PCA are not present in the tree:")
  print(rownames(morphoPCs)[!rownames(morphoPCs) %in% tree$tip.label])
  morphoPCs <- morphoPCs[tree$tip.label,]
}

# independent PCs
for (i in 1:ncol(morphoPCs)){
  assign(paste0("morphoPC", i), setNames(morphoPCs[,i], rownames(morphoPCs))[tree$tip.label])
}
