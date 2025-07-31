################################################################################
##### -------------------------------------------------------------------- ##### 
#####                           Data preparation                           ##### 
##### -------------------------------------------------------------------- ##### 
################################################################################
## set workfolder
path = "~/Unil/Research/5_Damselfish_evo/"
setwd(path) ###########

# load custom functions
source("./Rscripts/custom_functions.R")

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

#imgTrans.list_res <- readRDS("./rdata/damselfish_transimages_res300.rds")

### Load tree for visualizing on PCA plot
#damsel_tree <- ape::read.nexus("McCord2021_SupportingInfoFinal/S1_TreeFile_BestML.tre")
damsel_tree <- read.tree("McCord2021_SupportingInfoFinal/S3_TreeFile_DamselOnlyTimeTree.phy")
damsel_tree$tip.label <- sub("[[:punct:]]$", "", damsel_tree$tip.label) # remove special characters at the end of the name
Ntip(damsel_tree)

## load eco traits
eco_traits <- read.csv("McCord2021_SupportingInfoFinal/S3_Table_Traits.csv", row.names = 1)
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
# load all taxa from fishbase
pomacentridae_fishbase <- read.csv("./data/pomacentridae_taxons.csv")
# Search for the Pomacentridae family
pomacentridae_itis <- taxize::children(uniq_genus, db = "itis")
pomacentridae_itis <- do.call(rbind, pomacentridae_itis)
pomacentridae_taxa <- unique(pomacentridae_itis$taxonname)
# Get all pomacentridae species
all_pomacentridae <- unique(c(pomacentridae_fishbase$Scientific_name, pomacentridae_taxa))
#write.table(all_pomacentridae, "./data/all_pomacentridae_species_names.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
all_pomacentridae <- apply(read.table("./data/all_pomacentridae_species_names.txt"), 1, paste, collapse = "_")

##############
# Check species name matching between tree and all taxons db
##############
check.names <- as.data.frame(t(sapply(damsel_tree$tip.label, function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
mispelled_tree <- mispelled_tree[c(1, 2, 4, 6, 7),]
damsel_tree$tip.label[damsel_tree$tip.label %in% rownames(mispelled_tree)] = unlist(mispelled_tree$suggested)

##############
# Check species name matching between eco data and all taxons db
##############
check.names <- as.data.frame(t(sapply(rownames(eco_traits), function(i) check_spelling(i, all_pomacentridae))))
mispelled_eco <- check.names[check.names$score != 1,]
mispelled_eco
mispelled_eco <- mispelled_eco[c(3, 4, 6),]
rownames(eco_traits)[rownames(eco_traits) %in% rownames(mispelled_eco)] = unlist(mispelled_eco$suggested)

##############
# Check species name matching between images and all taxons db
##############
check.names <- as.data.frame(t(sapply(species_names, function(i) check_spelling(i, all_pomacentridae))))
mispelled_img <- check.names[check.names$score != 1,]
mispelled_img
mispelled_img <- mispelled_img[c(3, 4, 5, 7, 8),]
for (i in 1:nrow(mispelled_img)){
  img_species_names[img_species_names %in% rownames(mispelled_img)[i]] = mispelled_img$suggested[[i]]
}
species_names <- unique(img_species_names)
#names(imgTrans.list_res)[names(imgTrans.list_res) %in% rownames(mispelled_img)] = unlist(mispelled_img$suggested)

# get all species from all datasets
all_species_names <- unique(c(names(imgTrans.list_res), damsel_tree$tip.label, rownames(eco_traits)))
all_species_names <- unique(c(species_names, damsel_tree$tip.label, rownames(eco_traits)))

# missing data
missing_traits <- all_species_names[!all_species_names %in% rownames(eco_traits)]
missing_imgs <- all_species_names[!all_species_names %in% species_names]
missing_tree <- all_species_names[!all_species_names %in% damsel_tree$tip.label]

# select species with all data
selected_species <- all_species_names[!all_species_names %in% unique(c(missing_traits, missing_imgs, missing_tree))]
length(selected_species)
#write.table(selected_species, file = "~/Unil/Research/5_Damselfish_evo/data/listOfSpecies_new.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
selected_species <- read.table(file = "~/Unil/Research/5_Damselfish_evo/data/listOfSpecies_new.txt")[,1]

##########################################
# subset tree with selected species
##########################################
damsel_tree_subset <- keep.tip(damsel_tree, selected_species)
plot(damsel_tree_subset)
#write.tree(damsel_tree_subset, file = "./data/damsel_subset.tre")
damsel_tree_subset <- read.tree("./data/damsel_subset.tre")

##########################################
# subset eco data with selected species
##########################################
eco_traits_subset <- eco_traits[selected_species,]
# write.csv(eco_traits_subset, "./data/eco_traits_subset.csv", row.names = T)
eco_traits_subset <- read.csv("~/Unil/Research/5_Damselfish_evo/data/eco_traits_subset.csv", row.names = 1)
rownames(eco_traits_subset) %in% selected_species
nrow(eco_traits_subset)
##########################################
# subset images with selected species
##########################################
imgTrans.list_res_subset <- imgTrans.list_res[selected_species]
res <- dim(imgTrans.list_res_subset[[1]])[1]
saveRDS(imgTrans.list_res_subset, file = paste0("./rdata/damselfish_transimages_res", res, "_subset.rds"))

##########################################
# Get Fishbase traits
##########################################
library(rfishbase)
# morphological tratis
fb_dammorpho <- rfishbase::morphometrics(gsub("_", " ", selected_species)) # fishbase morpho data
fb_dammorpho <- fb_dammorpho %>%
  group_by(Species) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
fb_dammorpho$genus <- sapply(fb_dammorpho$Species, function(i) strsplit(i, " ")[[1]][1])

# ecological traits
fb_data <- fb_tbl("species") %>% 
  mutate(sci_name = paste(Genus, Species)) %>%
  filter(sci_name %in% gsub("_", " ", selected_species)) 

unique(paste0(fb_data$Genus, "_", fb_data$Species))
colnames(fb_data)

##########################################
# complete missing species traits 
##########################################
# Body shape
BodyShape <- setNames(fb_data$BodyShapeI, paste0(fb_data$Genus, "_", fb_data$Species))
miss.BodyShape.sp <- selected_species[!selected_species %in% names(BodyShape)]
# complete missing data from similar categorization in other species and other references
miss.BodyShape.sp <- setNames(c("short and / or deep", "short and / or deep", "short and / or deep", 
                                "fusiform / normal", "fusiform / normal", "short and / or deep", 
                                "fusiform / normal", "fusiform / normal", "fusiform / normal", 
                                "fusiform / normal", "fusiform / normal") , miss.BodyShape.sp)
BodyShape <- c(BodyShape, miss.BodyShape.sp)
BodyShape <- factor(BodyShape[selected_species], levels = unique(BodyShape))

# habitat (not accurately depicting the different habitats)
# habitat <- setNames(fb_data$DemersPelag, paste0(fb_data$Genus, "_", fb_data$Species))
# miss.habitat.sp <- selected_species[!selected_species %in% names(habitat)]
# # complete missing data from similar categorization in other species and other references
# miss.habitat.sp <- setNames(c("reef-associated", "reef-associated", "reef-associated", 
#                               "reef-associated", "reef-associated", "reef-associated", 
#                               "reef-associated", "reef-associated", "reef-associated",
#                               "reef-associated", "reef-associated") , miss.habitat.sp)
# habitat <- c(habitat, miss.habitat.sp)
# habitat <- factor(habitat[selected_species], levels = unique(habitat))

# depth range
DepthRange <- data.frame(minDepth = fb_data$DepthRangeShallow,
                         maxDepth = fb_data$DepthRangeDeep,
                         row.names = sub(" ", "_", fb_data$sci_name))

for (i in 1:nrow(DepthRange)) {
  if (!is.na(DepthRange$maxDepth[i])){
    if (DepthRange$maxDepth[i] > 150){
      DepthRange[i, "DepthRange"] = "bathyal"
    } else if (DepthRange$maxDepth[i] > 60){
      DepthRange[i, "DepthRange"] = "deep"
    } else if (DepthRange$maxDepth[i] > 20){
      DepthRange[i, "DepthRange"] = "mid-deep"
    } else if (DepthRange$maxDepth[i] <= 20){
      DepthRange[i, "DepthRange"] = "shallow"
    }
  }
}
DepthRange <- DepthRange[selected_species,]
plot(DepthRange$minDepth, DepthRange$maxDepth, pch = 21, bg = rep(viridis(4))[as.factor(DepthRange$range)])

###############################################################
# get more specific habitat information from Fishbase comments
###############################################################
comments <- lapply(fb_data$Comments, function(i) extract_nouns_with_adj(strsplit(gsub(" \\(Ref. [0-9]*\\)", "", i), "\\. ")[[1]][1]))
hab_comments <- grep("reef|lagoon|pool|outer|coastal|rubble|rock|shore|crevice|sand|coral|deep|shallow|mangrove|area|seagrass|freshwater|inshore|seaward|silt|patch|shelter|inlet|branch|estuary|habitat|anemone|Heteractis|cliff", unique(unlist(comments)), value = T)

hab_comments <- lapply(comments, function(i) i[i %in% hab_comments])
names(hab_comments) <- fb_data$sci_name

# write file
lapply(1:length(hab_comments), function(i) {
  cat(paste(c(names(hab_comments)[i], hab_comments[[i]]), collapse = "\t"), file = "./data/fb_hab_comments.csv", append = T)
  cat("\n", file = "./data/fb_hab_comments.csv", append = T)
})

# simplify habitat traits to only 5 states
# states are order so highly specific habitat have priority over more general ones
habitat_levels <- c("sea anemone", "freshwater", "seagrass-rubble-sand", 
                    "rocky-reef", "coral-reef", "")

unique(unlist(as.vector(fb_hab_data[,1:4])))
for (i in 1:4){
  fb_hab_data[fb_hab_data[,i] == "coral-reefs",i] = "coral-reef"
  fb_hab_data[fb_hab_data[,i] == "branching coral",i] = "coral-reef"
  fb_hab_data[fb_hab_data[,i] == "reef-flats",i] = "" # non-specific
  fb_hab_data[fb_hab_data[,i] == "reef",i] = "" # non-specific
  fb_hab_data[fb_hab_data[,i] == "coral-outcrop",i] = "coral-reef"
  fb_hab_data[fb_hab_data[,i] == "lagoon",i] = "" # non-specific
  fb_hab_data[fb_hab_data[,i] == "crevice",i] = "rocky-reef"
  fb_hab_data[fb_hab_data[,i] == "branching-coral",i] = "coral-reef"
  fb_hab_data[fb_hab_data[,i] == "coral-reef",i] = "coral-reef"
  fb_hab_data[fb_hab_data[,i] == "estuary",i] = "freshwater"
  fb_hab_data[fb_hab_data[,i] == "mangrove",i] ="seagrass-rubble-sand"
  fb_hab_data[fb_hab_data[,i] == "algal-reefs",i] = "seagrass-rubble-sand"
  fb_hab_data[fb_hab_data[,i] == "rubble",i] = "seagrass-rubble-sand"
  fb_hab_data[fb_hab_data[,i] == "seagrass",i] = "seagrass-rubble-sand"
  fb_hab_data[fb_hab_data[,i] == "seagrass-rubble",i] = "seagrass-rubble-sand"
  fb_hab_data[fb_hab_data[,i] == "sandy-bottom",i] = "seagrass-rubble-sand"
  fb_hab_data[fb_hab_data[,i] == "tide-pool",i] = "" # non-specific
}
rownames(fb_hab_data) <- sub(" ", "_", rownames(fb_hab_data))

# order habitats in habitat columns according to the specific habitat levels 
for (i in 1:nrow(fb_hab_data[,c("habitat1", "habitat2", "habitat3", "habitat4")])) {
  row <- fb_hab_data[i,c("habitat1", "habitat2", "habitat3", "habitat4")]
  row[duplicated(as.vector(row))] = ""
  row <- factor(row, levels = habitat_levels)
  row <- sort(row)
  fb_hab_data[i,c("habitat1", "habitat2", "habitat3", "habitat4")] = row
}

# save
write.csv(fb_hab_data, "./data/fb_hab_comments_simplified.csv", row.names = T)

# other interesting traits obtained from comments section of Fishbase
##############################################################################
# nocturnal vs diurnal
DielActivity <- setNames(unlist(sapply(fb_data$Comments, function(i) {
  text <- grep("nocturnal|diurnal", tolower(strsplit(i, " ")[[1]]), value = T)
  ifelse(length(text) == 0, NA, text) })),
  sub(" ", "_", fb_data$sci_name))
unique(DielActivity)
sum(!is.na(DielActivity))
# aggregate vs solitary
Gregariousness <- setNames(unlist(sapply(fb_data$Comments, function(i) {
  text <- grep("solitary|gregarious|aggregate|aggregations|single|school|multiple", tolower(strsplit(i, " ")[[1]]), value = T)
  ifelse(length(text) == 0, NA, text) })),
  sub(" ", "_", fb_data$sci_name))
unique(Gregariousness)
sum(!is.na(Gregariousness))
# herbivors - carnivors  - omnivors
TrophicGroup <- setNames(unlist(sapply(fb_data$Comments, function(i) {
  text <- grep("herbivour|carnivour|algae|zooplankton|phytoplankton|omnivour", tolower(strsplit(i, " ")[[1]]), value = T)
  ifelse(length(text) == 0, NA, text) })),
  sub(" ", "_", fb_data$sci_name))
unique(TrophicGroup)
sum(!is.na(TrophicGroup))

##############################################################################
##########################################
# put all data together
##########################################
for (i in 1:nrow(eco_traits_subset)){
  if (rownames(eco_traits_subset)[i] %in% names(BodyShape)) {
    eco_traits_subset[i, "BodyShape"] = BodyShape[rownames(eco_traits_subset)[i]]
  }
  if (rownames(eco_traits_subset)[i] %in% names(DielActivity)) {
    eco_traits_subset[i, "DielActivity"] = DielActivity[rownames(eco_traits_subset)[i]]
  }
  if (rownames(eco_traits_subset)[i] %in% names(Gregariousness)) {
    eco_traits_subset[i, "Gregariousness"] = Gregariousness[rownames(eco_traits_subset)[i]]
  }
  if (rownames(eco_traits_subset)[i] %in% names(TrophicGroup)) {
    eco_traits_subset[i, "TrophicGroup"] = TrophicGroup[rownames(eco_traits_subset)[i]]
  }
}

# add Depth range info
eco_traits_subset <- cbind(eco_traits_subset, DepthRange[rownames(eco_traits_subset),])
# add specific habitats info
eco_traits_subset <- cbind(eco_traits_subset, fb_hab_data[rownames(eco_traits_subset),])

##############################################################################
##########################################
# Data set make-up, cleaning and ordering
##########################################

# remove numeric diet trait
eco_traits_subset <- eco_traits_subset[,-1]
eco_traits_subset[eco_traits_subset == ""] = NA # convert empty string to NAs

# include the use of sea anemones in Farming trait as a mutualism. Farming would be a form of comensalism
idx <- which(apply(eco_traits_subset[,c("habitat1", "habitat2", "habitat3", "habitat4")], 1, function(i) "sea anemone" %in% i))
eco_traits_subset[idx,c("habitat1", "habitat2", "habitat3", "habitat4")]
eco_traits_subset$Farming[idx] = 2 # set reef fishes interacting with sea anemones with value 2 
# convert the numeric Farming categorical trait to binary trait symbiosis considering those with dependent ecological interactions and those that does not
eco_traits_subset$Symbiosis <- eco_traits_subset$Farming
eco_traits_subset$Farming[eco_traits_subset$Farming == 0] = "No" 
eco_traits_subset$Farming[eco_traits_subset$Farming == 1] = "Farming" 
eco_traits_subset$Farming[eco_traits_subset$Farming == 2] = "Mutualism" 

eco_traits_subset$Symbiosis[eco_traits_subset$Symbiosis != "No"] = "Yes"
# diel activity
unique(eco_traits_subset$DielActivity)
eco_traits_subset$DielActivity[eco_traits_subset$DielActivity == "diurnal."] = "diurnal"

# gregariousness
unique(eco_traits_subset$Gregariousness)
eco_traits_subset$Gregariousness[grep("school|aggreg", eco_traits_subset$Gregariousness)] = "gregarious"
eco_traits_subset$Gregariousness[grep("solitary|single", eco_traits_subset$Gregariousness)] = "solitary"

# trophic group
unique(eco_traits_subset$TrophicGroup)
eco_traits_subset$TrophicGroup[grep("algae|phyto", eco_traits_subset$TrophicGroup)] = "planktivore"
eco_traits_subset$TrophicGroup[grep("zoo", eco_traits_subset$TrophicGroup)] = "invertivore"

################################################################################
## SAVE ECO TRAITS DATASET
################################################################################
write.csv(eco_traits_subset, "./data/eco_traits_subset.csv", row.names = T)
################################################################################

#######
# !!!!!!
#######
# Further completion of missing data done manually on excel looking on reef fish databases
#######
# !!!!!!
########

eco_traits_subset <- read.csv("./data/eco_traits_subset.csv", row.names = 1)

eco_traits_subset <- eco_traits_subset[,-grep("source", colnames(eco_traits_subset))] # remove source info columns

# encode habitat1 and habitat2 together and simplify by specialist (occupying one) and generalist (occupying more than 1) habitats
eco_traits_subset$habitat2[which(eco_traits_subset$habitat2 == eco_traits_subset$habitat1)] = NA
eco_traits_subset$Habitat <- paste0(eco_traits_subset$habitat1, "&", eco_traits_subset$habitat2)
eco_traits_subset$Habitat <- sub("&NA", "", eco_traits_subset$Habitat)
unique(eco_traits_subset$Habitat)
eco_traits_subset$Habitat[grep("sea anemone", eco_traits_subset$Habitat)] = "specialists"
eco_traits_subset$Habitat[grep("freshwater", eco_traits_subset$Habitat)] = "specialists"
eco_traits_subset$Habitat[grepl("sand", eco_traits_subset$Habitat) & grepl("reef", eco_traits_subset$Habitat)] = "generalists"
eco_traits_subset$Habitat[grep("&", eco_traits_subset$Habitat)] = "generalists"
eco_traits_subset$Habitat[!eco_traits_subset$Habitat %in% c("specialists", "generalists")] = "specialists"

# encode Symbiosis trait with mutualists and farmers as "Yes" and the rest as "No"
eco_traits_subset$Symbiosis <- eco_traits_subset$Farming
eco_traits_subset$Symbiosis[eco_traits_subset$Symbiosis != "No"] = "Yes"

# add reproductive mode (data collected by Maelys)
eco_traits_subset[, "ReproMode"] = "dioecism"

# Clownfishes are protrandrous
eco_traits_subset[grep("Amphiprion", rownames(eco_traits_subset)), "ReproMode"] = "protandry"

# Tang et al 2021
protogyny_species <- c("Dascyllus_aruanus", "Dascyllus_melanurus", "Dascyllus_carneus", "Dascyllus_flavicaudus", "Dascyllus_marginatus", 
                       "Dascyllus_reticulatus", "Dascyllus_albisella", "Dascyllus_trimaculatus", "Pomacentrus_amboinensis")
eco_traits_subset[protogyny_species, "ReproMode"] = "protogyny"

# subset only traits with values for more than 95% of the species
eco_traits_subset <- eco_traits_subset[, apply(is.na(eco_traits_subset), 2, sum) < nrow(eco_traits_subset) * 0.05]

# add column with species name
eco_traits_subset$species <- rownames(eco_traits_subset)

################################################################################
## SAVE FINAL ECO TRAITS DATASET
################################################################################
write.csv(eco_traits_subset, "./data/eco_traits_subset_final.csv", row.names = T)
################################################################################
