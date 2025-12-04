################################################################################
##### -------------------------------------------------------------------- ##### 
#####                           Data preparation                           ##### 
##### -------------------------------------------------------------------- ##### 
################################################################################
## set workfolder
path = "~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/"
setwd(path) ###########

# load custom functions
source("./scripts/custom_functions.R")

# Load the pre-trained UDPipe model for English
ud_model <- udpipe_download_model(language = "english")
ud_model <- udpipe_load_model(ud_model$file_model)

# The purpose of this script is to combine all data to:
# 1- Check species names and rename according to a standard taxonomic classification (Fishbase database)
# 2- Subset species from the tree of which we have images and ecological traits
# 3- Obtain extra traits from fish databases (fishbase, rls, FoA, Reeflex...)
# 4-  Combine and save datasets

# load transformed images
img_filenames <- list.files("./images/TransIntImgPNGColorCorrected/")
img_species_names <- sapply(img_filenames, function(i) sub("-", "_", strsplit(i, "_")[[1]][1]))
img_species_names <- sapply(img_species_names, function(i) paste0(toupper(substr(i, 1, 1)), substr(i, 2, nchar(i))))
species_names <- unique(img_species_names)
length(species_names)

### Load tree for visualizing on PCA plot
#damsel_tree <- ape::read.nexus("McCord2021_SupportingInfoFinal/S1_TreeFile_BestML.tre")
damsel_tree <- read.tree("./data/McCord2021_SupportingInfoFinal/S3_TreeFile_DamselOnlyTimeTree.phy")
damsel_tree$tip.label <- sub("[[:punct:]]$", "", damsel_tree$tip.label) # remove special characters at the end of the name
Ntip(damsel_tree)

## load eco traits
eco_traits <- read.csv("./data/McCord2021_SupportingInfoFinal/S3_Table_Traits.csv", row.names = 1)
rownames(eco_traits) <- sub("[[:punct:]]$", "", rownames(eco_traits)) # remove special characters at the end of the name
nrow(eco_traits)

# get all Pomacentridae genus of the tree
genus <- sapply(damsel_tree$tip.label, function(i) strsplit(i, "_")[[1]][1])
uniq_genus <- unique(genus)

rownames(eco_traits)[!(rownames(eco_traits) %in% damsel_tree$tip.label)] # mismatch between tree and ecological traits in McCord 2021
damsel_tree$tip.label[!damsel_tree$tip.label %in% rownames(eco_traits)]  # mismatch between tree and ecological traits in McCord 2021

##############
# Load taxons
##############
# # load all taxa from fishbase
# pomacentridae_fishbase <- read.csv("./data/pomacentridae_taxons.csv")
# # Search for the Pomacentridae family
# pomacentridae_itis <- taxize::children(uniq_genus, db = "itis")
# pomacentridae_itis <- do.call(rbind, pomacentridae_itis)
# pomacentridae_taxa <- unique(pomacentridae_itis$taxonname)
# # Get all pomacentridae species
# all_pomacentridae <- unique(c(pomacentridae_fishbase$Scientific_name, pomacentridae_taxa))
# #write.table(all_pomacentridae, "./data/all_pomacentridae_species_names.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
all_pomacentridae <- apply(read.table("./data/all_pomacentridae_species_names.txt"), 1, paste, collapse = "_")

##############
# Check species name matching between tree and all taxons db
##############
check.names <- as.data.frame(t(sapply(damsel_tree$tip.label, function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
damsel_tree$tip.label <- sub("Stegastes_rectifraeneum", "Stegastes_rectifraenum", damsel_tree$tip.label)
damsel_tree$tip.label <- sub("Pycnochromis_lineata", "Pycnochromis_lineatus", damsel_tree$tip.label)
damsel_tree$tip.label <- sub("Chromis_sanctahelenae", "Chromis_sanctaehelenae", damsel_tree$tip.label)
damsel_tree$tip.label <- sub("Chrysiptera_maurinae", "Chrysiptera_maurineae", damsel_tree$tip.label)
damsel_tree$tip.label <- sub("Amblypomacentrus_tricincta", "Amblypomacentrus_tricinctus", damsel_tree$tip.label)

##############
# Check species name matching between eco data and all taxons db
##############
check.names <- as.data.frame(t(sapply(rownames(eco_traits), function(i) check_spelling(i, all_pomacentridae))))
mispelled_eco <- check.names[check.names$score != 1,]
mispelled_eco
rownames(eco_traits) <- sub("Stegastes_rectifraeneum", "Stegastes_rectifraenum", rownames(eco_traits))
rownames(eco_traits) <- sub("Chromis_sanctahelenae", "Chromis_sanctaehelenae", rownames(eco_traits))
rownames(eco_traits) <- sub("Chrysiptera_maurinae", "Chrysiptera_maurineae", rownames(eco_traits))
rownames(eco_traits) <- sub("Pycnochromis_lineata", "Pycnochromis_lineatus", rownames(eco_traits))

##############
# Check species name matching between images and all taxons db
##############
check.names <- as.data.frame(t(sapply(species_names, function(i) check_spelling(i, all_pomacentridae))))
mispelled_img <- check.names[check.names$score != 1,]
mispelled_img
img_species_names <- sub("Stegastes_rectifraeneum", "Stegastes_rectifraenum", img_species_names)
img_species_names <- sub("Chrysiptera_maurinae", "Chrysiptera_maurineae", img_species_names)
img_species_names <- sub("Plectroglyphidodon_leucozona", "Plectroglyphidodon_leucozonus", img_species_names)
img_species_names<- sub("Pomacentrus_caeruleopuncatus", "Pomacentrus_caeruleopunctatus", img_species_names)
img_species_names <- sub("Pycnochromis_lineata", "Pycnochromis_lineatus", img_species_names)
img_species_names <- sub("Stegastes_flaviatus", "Stegastes_flavilatus", img_species_names)
species_names <- unique(img_species_names)

# get all species from all datasets
all_species_names <- unique(c(species_names, damsel_tree$tip.label, rownames(eco_traits)))

# missing data
missing_traits <- all_species_names[!all_species_names %in% rownames(eco_traits)]
missing_imgs <- all_species_names[!all_species_names %in% species_names]
missing_tree <- all_species_names[!all_species_names %in% damsel_tree$tip.label]

# select species with all data
selected_species <- all_species_names[!all_species_names %in% unique(c(missing_traits, missing_imgs, missing_tree))]
length(selected_species)

# write.table(selected_species, file = "~/Unil/Research/5_Damselfish_evo/data/listOfSpecies.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
selected_species <- read.table(file = "./data/listOfSpecies.txt")[,1]

##########################################
# subset tree with selected species
##########################################
damsel_tree_subset <- keep.tip(damsel_tree, selected_species)
plot(damsel_tree_subset)
write.tree(damsel_tree_subset, file = "./data/damsel_subset.tre")

##########################################
# subset eco data with selected species
##########################################
eco_traits_subset <- eco_traits[selected_species,c("DietEcotype", "Farming")]
# write.csv(eco_traits_subset, "./data/eco_traits_subset.csv", row.names = T)
# eco_traits_subset <- read.csv("~/Unil/Research/5_Damselfish_evo/data/eco_traits_subset.csv", row.names = 1)
nrow(eco_traits_subset)

##########################################
# Get Habitat information from Fishbase traits
##########################################
library(rfishbase)
# ecological traits
fb_data <- fb_tbl("species") %>% 
  dplyr::mutate(sci_name = paste(Genus, Species)) %>%
  dplyr::filter(sci_name %in% gsub("_", " ", selected_species)) 

unique(paste0(fb_data$Genus, "_", fb_data$Species))
colnames(fb_data)

###############################################################
# Extract and simplify habitat information from FishBase comments
###############################################################

habitat_levels <- c("sea anemone", "freshwater", "seagrass-rubble-sand", "coral-reef", "rocky-reef")
nhabitats <- length(habitat_levels)

# Recode habitats using dictionary
recode_dict <- c(
  # === sea anemone ===
  "anemone" = "sea anemone",
  "Heteractis" = "sea anemone",
  
  # === freshwater ===
  "freshwater" = "freshwater",
  "freshwater-stream" = "freshwater",
  "estuary" = "freshwater",
  
  # === seagrass-rubble-sand ===
  "mangrove" = "seagrass-rubble-sand",
  "algal-reefs" = "seagrass-rubble-sand",
  "rubble" = "seagrass-rubble-sand",
  "open-rubble" = "seagrass-rubble-sand",
  "coral-rubble" = "seagrass-rubble-sand",
  "coral-rubbles" = "seagrass-rubble-sand",
  "sand" = "seagrass-rubble-sand",
  "sandy" = "seagrass-rubble-sand",
  "sandy-area" = "seagrass-rubble-sand",
  "sandy-beach" = "seagrass-rubble-sand",
  "sandy-bottom" = "seagrass-rubble-sand",
  "sand-flat" = "seagrass-rubble-sand",
  "deep-sand" = "seagrass-rubble-sand",
  "seagrass" = "seagrass-rubble-sand",
  
  # === rocky-reef ===
  "rock" = "rocky-reef",
  "rocky" = "rocky-reef",
  "rocky-bottom" = "rocky-reef",
  "rocky-reef" = "rocky-reef",
  "rocky-reefs" = "rocky-reef",
  "rocky-inshore" = "rocky-reef",
  "rocky-outcrop" = "rocky-reef",
  "rocky-patch" = "rocky-reef",
  "rocky-area" = "rocky-reef",
  "rocky-substrate" = "rocky-reef",
  "rocky-coast" = "rocky-reef",
  "rocky-shore" = "rocky-reef",
  "rocky-outcrop" = "rocky-reef",
  "rocky-shore" = "rocky-reef",
  "rocky-substrate" = "rocky-reef",
  "rocky-patch" = "rocky-reef",
  "rocky-area" = "rocky-reef",
  "rocky-coast" = "rocky-reef",
  "rocky-shore" = "rocky-reef",
  "rocky-inshore" = "rocky-reef",
  "crevice" = "rocky-reef",
  "large-crevice" = "rocky-reef",
  "cliff" = "rocky-reef",
  "small-rock" = "rocky-reef",
  "dead-rock" = "rocky-reef",
  "inshore-boulder" = "rocky-reef",
  
  # === coral-reef ===
  "coral-reef" = "coral-reef",
  "coral-reefs" = "coral-reef",
  "coral-outcrop" = "coral-reef",
  "coral-outcropping" = "coral-reef",
  "coral-colony" = "coral-reef",
  "coral-area" = "coral-reef",
  "coral-bottom" = "coral-reef",
  "coral-head" = "coral-reef",
  "coral-thicket" = "coral-reef",
  "coral-basis" = "coral-reef",
  "coral-formation" = "coral-reef",
  "coral-growth" = "coral-reef",
  "coral-cover" = "coral-reef",
  "coral-patch" = "coral-reef",
  "coral-rock" = "coral-reef",
  "coral-habitat" = "coral-reef",
  "coral-reeg" = "coral-reef",
  "coral-rubbles" = "coral-reef",
  "live-coral" = "coral-reef",
  "soft-coral" = "coral-reef",
  "hard-coral" = "coral-reef",
  "mixed-coral" = "coral-reef",
  "rich-coral" = "coral-reef",
  "weedy-reef" = "coral-reef",
  "non-coral-reef" = "coral-reef",
  "black-coral" = "coral-reef",
  "dead-reefs" = "coral-reef",
  "dead-reef" = "coral-reef",
  "dead-coral" = "coral-reef",
  "reefs" = "coral-reef",
  
  # === non-specific (collapse to "") ===
  "lagoon" = NA,
  "deep-lagoon" = NA,
  "shallow-lagoon" = NA,
  "silty-lagoon" = NA,
  "clear-lagoon" = NA,
  "outer-lagoon" = NA,
  "reef" = NA,
  "reef-flat" = NA,
  "reef-flats" = NA,
  "reef-area" = NA,
  "coastal" = NA,
  "coastal-area" = NA,
  "coastal-reef" = NA,
  "coastal-reefs" = NA,
  "coastal-embayment" = NA,
  "coastal-water" = NA,
  "coastal-region" = NA,
  "coastal-reefs" = NA,
  "shore" = NA,
  "shoreline" = NA,
  "flat-shoreline" = NA,
  "outer-reef" = NA,
  "outer-reefs" = NA,
  "outer-slope" = NA,
  "outer-edge" = NA,
  "outer-rocky" = NA,
  "outer-channel" = NA,
  "outer-periphery" = NA,
  "outermost-reach" = NA,
  "outer-slope" = NA,
  "inshore" = NA,
  "inshore-area" = NA,
  "inshore-reef" = NA,
  "inshore-reefs" = NA,
  "inner-reef" = NA,
  "intertidal-reef" = NA,
  "offshore-reef" = NA,
  "offshore-reefs" = NA,
  "offshore-area" = NA,
  "offshore-platform" = NA,
  "offshore-trawl" = NA,
  "offshore-trawling" = NA,
  "deep-trawl" = NA,
  "deeper-trawl" = NA,
  "deep-trawling" = NA,
  "deep-reef" = NA,
  "deep-reefs" = NA,
  "deep-slope" = NA,
  "deep-patch" = NA,
  "deep-lagoon" = NA,
  "deeper-patch" = NA,
  "deeper-water" = NA,
  "deeper-surge" = NA,
  "shallow-reef" = NA,
  "shallow-reeg" = NA,
  "shallow-area" = NA,
  "shallow-water" = NA,
  "shallow-depth" = NA,
  "shallow-habitat" = NA,
  "shallow-rock" = NA,
  "shallower" = NA,
  "shelter" = NA,
  "shelter-water" = NA,
  "pool" = NA,
  "surge-pool" = NA,
  "surge-area" = NA,
  "silt" = NA,
  "silty-area" = NA,
  "silty-water" = NA,
  "silty-reef" = NA,
  "silty-region" = NA,
  "silt-affected-reef" = NA,
  "branch" = NA,
  "branching" = NA,
  "branching-coral" = "coral-reef",  # treat as coral-reef when specific
  "open-branching" = "coral-reef",
  "large-branching" = "coral-reef",
  "habitat" = NA,
  "habitats" = NA,
  "typical-habitat" = NA,
  "usual-habitat" = NA,
  "rich-habitat" = NA,
  "areas" = NA,
  "area" = NA,
  "-rich-area" = NA,
  "rich-area" = NA,
  "large-area" = NA,
  "bottom-area" = NA,
  "flat-area" = NA,
  "dead-area" = NA,
  "open-area" = NA,
  "calm-area" = NA,
  "live-patch" = NA,
  "small-patch" = NA,
  "patch" = NA,
  "outer-slope" = NA,
  "outer-channel" = NA,
  "outer-edge" = NA,
  "outer-periphery" = NA,
  "outermost-reach" = NA
)

# Extract nouns/adjectives from first sentence of comments
comments <- lapply(fb_data$Comments, function(txt) {
  clean <- gsub(" \\(Ref. [0-9]*\\)", "", txt)
  first_sentence <- strsplit(clean, "\\. ")[[1]][1]
  extract_nouns_with_adj(first_sentence)
})

# Habitat-related terms
hab_terms <- grep(
  "reef|lagoon|pool|outer|coastal|rubble|rock|shore|crevice|sand|coral|deep|shallow|mangrove|area|seagrass|freshwater|inshore|seaward|silt|patch|shelter|inlet|branch|estuary|habitat|anemone|Heteractis|cliff",
  unique(unlist(comments)), value = TRUE
)

# Keep only habitat-related terms and recode
hab_comments <- lapply(comments, function(i) {
  term = i[i %in% hab_terms]
  unique_recoded_terms = unique(recode_dict[term])
  na.exclude(unique_recoded_terms)
  })
names(hab_comments) <- fb_data$sci_name

# Convert to data.frame 
fb_hab_data <- do.call(rbind, lapply(names(hab_comments), function(sp) {
  vals <- hab_comments[[sp]]
  vals <- c(vals, rep(NA, nhabitats - length(vals)))[1:nhabitats]
  data.frame(t(vals), stringsAsFactors = FALSE, row.names = gsub(" ", "_", sp))
}))
colnames(fb_hab_data) <- paste0("habitat", 1:nhabitats)

habitat_specificity <- apply(fb_hab_data, 1, function(i) sum(!is.na(i)))
habitat_specificity <- habitat_specificity[habitat_specificity > 0]

fb_hab_data <- fb_hab_data[names(habitat_specificity),]

# Order habitats by specificity based on habitat levels
fb_hab_data <- t(apply(fb_hab_data, 1, function(row) {
  row <- factor(row, levels = habitat_levels)
  row <- sort(row, na.last = TRUE)
  row
}))
colnames(eco_traits_subset_old)
eco_traits_subset_old <- eco_traits_subset_old[, c("DietEcotype", "Farming", "habitat1", "source1", "source2", "source3")]
colnames(eco_traits_subset_old)[2:3] <- c("Symbiosis", "Habitat")

# Save simplified version
write.csv(fb_hab_data, "./data/fb_hab_comments_simplified.csv", row.names = TRUE)

# add specific habitats info
eco_traits_subset$Habitat <- NA
eco_traits_subset[rownames(fb_hab_data),"Habitat"] = fb_hab_data[,1]
eco_traits_subset$Habitat <- factor(eco_traits_subset$Habitat, levels = habitat_levels)

##############################################################################
##########################################
# Data set make-up, cleaning and ordering
##########################################

# include the use of sea anemones in Farming trait as a mutualism. Farming would be a form of commensalism
eco_traits_subset$Habitat == "sea anemone"
eco_traits_subset$Habitat[grep("Amphiprion", rownames(eco_traits_subset))] = "sea anemone"
eco_traits_subset[c("Amblypomacentrus_breviceps", "Amblypomacentrus_clarus", "Dascyllus_albisella","Dascyllus_trimaculatus"), "Habitat"] = "sea anemone"

eco_traits_subset$Symbiosis[eco_traits_subset$Habitat == "sea anemone"] = 2 # set reef fishes interacting with sea anemones with value 2 
# convert the numeric Farming categorical trait to binary trait symbiosis considering those with dependent ecological interactions and those that does not
eco_traits_subset$Symbiosis[eco_traits_subset$Symbiosis == 0] = "No" 
eco_traits_subset$Symbiosis[eco_traits_subset$Symbiosis == 1] = "Farming" 
eco_traits_subset$Symbiosis[eco_traits_subset$Symbiosis == 2] = "Mutualism" 

################################################################################
## SAVE ECO TRAITS DATASET
################################################################################
write.csv(eco_traits_subset, "./data/eco_traits_subset.csv", row.names = T)
################################################################################

#######
# !!!!!!
#######
# Further completion of missing data done manually on excel looking on reef fish databases 
# Saved as: "./data/eco_traits_subset_final.csv"
#######
# !!!!!!
########
