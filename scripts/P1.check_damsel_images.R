lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", 
         "lme4", "mgcv", "phytools", "corHMM",
         "lmtest", "ColorAR", "castor", "l1ou",
         "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

## Settings
###########
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

# load custom functions
source("./scripts/custom_functions.R")

getSpNames <- function(list){
  list <- sub(" ", "_", list)
  list <- sub("-", "_", list)
  unlist(lapply(strsplit(list, "_"), function(i) paste(i[1], i[2])))
}

# list of species
listOfSpecies <- read.table("./data/listOfSpecies.txt")[,1]

# all pomacentridae species names
all_pomacentridae <- apply(read.table("./data/all_pomacentridae_species_names.txt"), 1, paste, collapse = "_")

# read tree
damsel_tree <- read.tree(".data/McCord2021_SupportingInfoFinal/S3_TreeFile_DamselOnlyTimeTree.phy")
tree_sps <- sub("_", " ", damsel_tree$tip.label)

# images
all_damsel_pics <- list.files("./images/all_selected_images/")
all_damsel_landmarks <- list.files("./images/all_selected_landmarks/")

all_damsel_pics <- setNames(all_damsel_pics, getSpNames(all_damsel_pics))
all_damsel_landmarks <- setNames(all_damsel_landmarks, getSpNames(all_damsel_landmarks))

# check names on tree
check.names <- as.data.frame(t(sapply(damsel_tree$tip.label, function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
damsel_tree$tip.label <- sub("Stegastes_rectifraenum", "Stegastes_rectifraeneum", damsel_tree$tip.label)
damsel_tree$tip.label <- sub("Pycnochromis_lineata", "Pycnochromis_lineatus", damsel_tree$tip.label)
damsel_tree$tip.label <- sub("Chromis_sanctahelenae", "Chromis_sanctaehelenae", damsel_tree$tip.label)
damsel_tree$tip.label <- sub("Chrysiptera_maurinae", "Chrysiptera_maurineae", damsel_tree$tip.label)
damsel_tree$tip.label <- sub("Amblypomacentrus_tricincta", "Amblypomacentrus_tricinctus", damsel_tree$tip.label)

# check names selection pictures
check.names <- as.data.frame(t(sapply(names(all_damsel_pics), function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
names(all_damsel_pics)[names(all_damsel_pics) %in% "Stegastes flaviatus"] <- "Stegastes flavilatus"
names(all_damsel_pics)[names(all_damsel_pics) %in% "Plectroglyphidodon leucozona"] <- "Plectroglyphidodon leucozonus"
names(all_damsel_pics)[names(all_damsel_pics) %in% "pomacentrus maafu"] <- "Pomacentrus maafu"
names(all_damsel_pics)[names(all_damsel_pics) %in% "Chrysiptera maurinae"] <- "Chrysiptera maurineae"
names(all_damsel_pics)[names(all_damsel_pics) %in% "Pomacentrus caeruleopuncatus"] <- "Pomacentrus caeruleopunctatus"
names(all_damsel_pics)[names(all_damsel_pics) %in% "Pycnochromis lineata"] <- "Pycnochromis lineatus"

# check names new landmarked images
check.names <- as.data.frame(t(sapply(names(all_damsel_landmarks), function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
names(all_damsel_landmarks)[names(all_damsel_landmarks) %in% "Stegastes flaviatus"] <- "Stegastes flavilatus"
names(all_damsel_landmarks)[names(all_damsel_landmarks) %in% "Plectroglyphidodon leucozona"] <- "Plectroglyphidodon leucozonus"
names(all_damsel_landmarks)[names(all_damsel_landmarks) %in% "pomacentrus maafu"] <- "Pomacentrus maafu"
names(all_damsel_landmarks)[names(all_damsel_landmarks) %in% "Chrysiptera maurinae"] <- "Chrysiptera maurineae"
names(all_damsel_landmarks)[names(all_damsel_landmarks) %in% "Pomacentrus caeruleopuncatus"] <- "Pomacentrus caeruleopunctatus"
names(all_damsel_landmarks)[names(all_damsel_landmarks) %in% "Pycnochromis lineata"] <- "Pycnochromis lineatus"

# species with image not in tree
names(all_damsel_pics)[!names(all_damsel_pics) %in% sub("_", " ", damsel_tree$tip.label)]
# species in tree with no images
missing_sp_pic <- damsel_tree$tip.label[!sub("_", " ", damsel_tree$tip.label) %in% names(all_damsel_pics)]

all_damsel_pics_subset <- all_damsel_pics[names(all_damsel_pics) %in% sub("_", " ", damsel_tree$tip.label)]
all_damsel_pics_subset_tab <- table(names(all_damsel_pics_subset))
nrow(all_damsel_pics_subset_tab)

missing_sp_pic %in% names(all_damsel_pics) # images for two missing species not found

