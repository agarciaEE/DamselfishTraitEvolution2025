################################################################################
####                                                                      ######
#### ------------      12. Biogeography and phylogenetic    ------------- ######
#### ------------                community ecology          ------------- ######
####                                                                      ######
################################################################################
library(dplyr)
library(sf)
library(ggplot2)
library(phytools)
library(geiger)
library(ape)
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
# seems to sometimes fail on simple problems (2-3 parameters)
library(rexpokit)
library(cladoRcpp)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)
library(rfishbase)
library(rgbif)
library(tidyr)
library(geojsonio)
library(sp)
library(broom)
library(ggplot2)
library(ape)
library(phytools)
library(geiger)
library(phangorn)
library(ColorAR)
library(treeio)
library(MultinomialCI) 

setwd("~/Unil/Research/5_Damselfish_evo/")

source("./Rscripts/custom_functions.R")

pomacentridae_species <- read.csv("./data/all_pomacentridae_species_names.txt", header = F)[,1]
pomacentridae_species <- sub(" ", "_", pomacentridae_species)

damsel_tree_subset <- read.tree("./data/damsel_subset.tre")
Ntip(damsel_tree_subset)

selected_species <- as.character(read.table("./data/listOfSpecies.txt")[,1])


################################################################################
# 12.1. Ancestral Area Reconstruction
################################################################################

#########################################
## --- Pomacentridae distributions ----##
#########################################
load("~/Unil/Research/5_Damselfish_evo/BioGeo/input/PA_acanth.rda")
names(PA_acanth)

selected_genus <- sapply(names(PA_acanth)[-c(1:2)], function(i) strsplit(i, "_")[[1]][1] %in% unique(genus))
PA_acanth <- PA_acanth[, c("Longitude", "Latitude", names(PA_acanth)[-c(1:2)][selected_genus])]

# Check species name matching between distributions and all taxons db
check.names <- as.data.frame(t(sapply(names(PA_acanth)[-c(1:2)], function(i) check_spelling(i, all_pomacentridae))))
mispelled_dis <- check.names[check.names$score != 1,]
mispelled_dis
mispelled_dis <- mispelled_dis[mispelled_dis$score > 0.5,]
names(PA_acanth)[-c(1:2)][names(PA_acanth)[-c(1:2)] %in% rownames(mispelled_dis)] = unlist(mispelled_dis$suggested)

names(PA_acanth)[-c(1:2)][!names(PA_acanth)[-c(1:2)] %in% all_species_names] # not in distribution data
all_species_names[!all_species_names %in% names(PA_acanth)[-c(1:2)]] # not in phylo data

PA_acanth_selection <- PA_acanth[, c("Longitude", "Latitude", names(PA_acanth)[-c(1:2)][names(PA_acanth)[-c(1:2)] %in% all_species_names] )]
#write.csv(PA_acanth_selection, file = "~/Unil/Research/5_Damselfish_evo/BioGeo/input/pomacentridae_distributions.csv", row.names = F)


# subset pomacentridae species
pomacentridae_genus <- unique(sapply(pomacentridae_species, function(i) strsplit(i, "_")[[1]][1]))
acanth_genus <- sapply(colnames(PA_acanth), function(i) strsplit(i, "_")[[1]][1])
pomacentridae_dis <- cbind(PA_acanth[,1:2], PA_acanth[, acanth_genus %in% pomacentridae_genus])

# check mispellings
check.names <- as.data.frame(t(sapply(colnames(pomacentridae_dis)[-c(1:2)], function(i) check_spelling(i, pomacentridae_species))))
mispelled_dis <- check.names[check.names$score != 1,]
colnames(pomacentridae_dis)[colnames(pomacentridae_dis) == "Chromis_brevirostris"] = "Azurina_brevirostris"
colnames(pomacentridae_dis)[colnames(pomacentridae_dis) == "Chromis_ovatiformes"] = "Pycnochromis_ovatiformis"
colnames(pomacentridae_dis)[colnames(pomacentridae_dis) == "Chromis_lepidolepis"] = "Azurina_lepidolepis"

write.csv(pomacentridae_dis, "./BioGeo/input/pomacentridae_distributions.csv", row.names = F)

# proportion of species distribution over total number of species
(ncol(pomacentridae_dis)-2) / length(pomacentridae_species) 

colnames(pomacentridae_dis)[colnames(pomacentridae_dis) %in% selected_species]
colnames(pomacentridae_dis)[-c(1:2)][!colnames(pomacentridae_dis)[-c(1:2)] %in% selected_species] # species distributions without tree or ecological data

# missing species distribution
missing_species_dis <- selected_species[!selected_species %in% colnames(pomacentridae_dis)] # missing species distributions

#########################################
## ---- Pomacentridae occurrences -----##
#########################################

# Search for occurrences of Pomacentridae
pomacentrida_occs <- rgbif::occ_search(
  scientificName = "Pomacentridae",
  limit = 10000,  # Adjust the limit as needed to retrieve more occurrences
  return = "data"
)
head(pomacentrida_occs$data)
pomacentrida_occ_data <- pomacentrida_occs$data[, c("species", "decimalLongitude", "decimalLatitude")]

# get occurrences from iNaturalist
iNat_obs <- read.csv("./BioGeo/input/observations-324213.csv")
iNat_obs <- iNat_obs[,c("scientific_name", "longitude", "latitude")]

iNat_obs2 <- read.csv("./BioGeo/input/observations-415875.csv")
iNat_obs2 <- iNat_obs2[,c("scientific_name", "longitude", "latitude")]

iNat_obs <- rbind(iNat_obs, iNat_obs2)
iNat_obs <- unique(iNat_obs)

# match column names
colnames(iNat_obs) <- c("species", "lon", "lat")
colnames(pomacentrida_occ_data) <- c("species", "lon", "lat")

# merge
pomacentrida_occ_data <- rbind(pomacentrida_occ_data, iNat_obs)

# remove occurences of only genus taxed
pomacentrida_occ_data <- pomacentrida_occ_data[!pomacentrida_occ_data$species %in% c(uniq_genus, "Pomacentridae", "Pomacentrinae"),]
pomacentrida_occ_data <- pomacentrida_occ_data[grep("_", sub(" ", "_", pomacentrida_occ_data$species)),]

# "Abudefduf caudobimaculatus" is synonymous of "Abudefduf vaigiensis"
x <- pomacentrida_occ_data[pomacentrida_occ_data$species == "Abudefduf vaigiensis",]
x$species = "Abudefduf caudobimaculatus"
pomacentrida_occ_data <- rbind(pomacentrida_occ_data, x)

pomacentrida_occ_data <- unique(pomacentrida_occ_data)
pomacentrida_occ_data <- na.exclude(pomacentrida_occ_data)

# Check species name matching between tree and all taxons db
check.names <- as.data.frame(t(sapply(unique(pomacentrida_occ_data$species), function(i) check_spelling(i, pomacentridae_species))))
mispelled_occs <- check.names[check.names$score != 1,]
mispelled_occs <- mispelled_occs[c(6,8),]
for (i in 1:nrow(mispelled_occs)){
  pomacentrida_occ_data$species[pomacentrida_occ_data$species == rownames(mispelled_occs)[i]] = sub("_", " ", unlist(mispelled_occs$suggested)[i])
}

# Save the occurrence data to a CSV file
write.csv(pomacentrida_occ_data, "./BioGeo/input/pomacentridae_occurrences.csv", row.names = FALSE)

#########################################
# missing species distributions recovered with occurrences
missing_species_dis[sub("_", " ", missing_species_dis) %in% unique(pomacentrida_occ_data$species)]

# still missing species
missing_species_occs <- missing_species_dis[!sub("_", " ", missing_species_dis) %in% unique(pomacentrida_occ_data$species)]

#########################################

###################################
## Using FAO areas from fishbase ##
###################################
# Function to get distribution information for a species
get_distribution <- function(speciesList) {
  
  info <- rfishbase::distribution(speciesList)
  
  # fix lon lat coordinates
  info$LatDeg <- sapply(1:length(info$LatDeg), function(i) if (!is.na(info$N_S[i]) & info$N_S[i] == "S") { -info$LatDeg[i] } else { info$LatDeg[i] })
  info$LongDeg <- sapply(1:length(info$LongDeg), function(i) if (!is.na(info$E_W[i]) & info$E_W[i] == "W") { -info$LongDeg[i] } else { info$LongDeg[i] })
  
  distribution <- data.frame(species = info$Species,
                             FAO_area = info$FAO,
                             Longitude = info$LongDeg,
                             Latitude = info$LatDeg)
  return(distribution)
}

# get all pomacentridae species (496) regional data
all_damsel_reg <- get_distribution(sub("_", " ", pomacentridae_species))

# download FAO json areas
download.file("http://www.fao.org/fishery/geoserver/fifao/ows?service=WFS&request=GetFeature&version=1.0.0&typeName=fifao:FAO_AREAS_CWP&outputFormat=json", dest="FAO.json")
spdf <- geojson_read("FAO.json",  what = "sp")

# Split FAO_Area column into Area_Number and Area_Name
fao_areas <- separate(fao_areas, FAO_Area, into = c("Area_Number", "Area_Name"), sep = " - ")

spdf_fortified <- tidy(spdf)
## make up some data to go with ...
fake_fish <- data.frame(id = as.character(1:324), value = rnorm(324))
spdf2 <- spdf_fortified %>% left_join(fake_fish, by = "id")
ggplot() +
  geom_polygon(data = spdf2, aes( x = long, y = lat, group = group,
                                  fill = value), color="grey") +
  scale_fill_viridis_c() +
  theme_void() +
  theme(plot.background = element_rect(fill = 'lightgray', colour = NA)) +
  coord_map() +
  coord_sf(crs = "+proj=cea +lon_0=0 +lat_ts=45") ## Gall projection


# get fishbase species regions
damsel_reg_FAO <- get_distribution(c(sub("_", " ", selected_species), "Abudefduf vaigiensis"))
damsel_reg_FAO$species[is.na(damsel_reg_FAO$FAO_area)] # missing species data
damsel_reg_FAO$species[damsel_reg_FAO$species == "Abudefduf vaigiensis"] = "Abudefduf caudobimaculatus" 
damsel_reg_FAO[damsel_reg_FAO$species == "Amblypomacentrus annulatus", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Amblypomacentrus tricinctus", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Chromis anadema", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Chromis enchrysura", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Chromis katoi", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Chromis kennensis", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Chrysiptera caesifrons", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Chrysiptera ellenae", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Chrysiptera leucopoma", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Chrysiptera maurineae", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Chrysiptera papuensis", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Dascyllus abudafur", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Pomacentrus flavioculus", 2:4] = c("Indian Ocean, Western", 55, -7.5)
damsel_reg_FAO[damsel_reg_FAO$species == "Stegastes sanctipauli", 2:4] = c("Indian Ocean, Western", 55, -7.5)

damsel_reg_tab_FAO <- t(table(damsel_reg_FAO$FAO_area, damsel_reg_FAO$species)) # convert to a table

species_nreg <- apply(damsel_reg_tab_FAO, 1, sum, na.rm = T) # number of regions occpied by each species
missing_species_reg <- names(species_nreg[species_nreg == 0]) # missing species regions data

FAO_regions <- colnames(damsel_reg_tab_FAO)
plot(damsel_sdm[[sub(" ", "_", missing_species_reg[sub(" ", "_", missing_species_reg) %in% colnames(damsel_dis)])]])

# Complete missing species regions
damsel_reg_tab_FAO["Chromis enchrysura", c("Atlantic, Northwest", "Atlantic, Southwest", "Atlantic, Western Central")] = 1
damsel_reg_tab_FAO["Stegastes sanctipauli", "Atlantic, Western Central"] = 1
damsel_reg_tab_FAO["Abudefduf caudobimaculatus", c("Indian Ocean, Eastern", "Indian Ocean, Western", 
                                               "Pacific, Northwest" , "Pacific, Western Central",  "Pacific, Southwest")] = 1
damsel_reg_tab_FAO["Amblypomacentrus annulatus", "Indian Ocean, Western"] = 1
damsel_reg_tab_FAO["Amblypomacentrus tricinctus", c("Pacific, Northwest" , "Pacific, Western Central",  "Pacific, Southwest")] = 1
damsel_reg_tab_FAO["Chromis anadema", "Pacific, Northwest"] = 1
damsel_reg_tab_FAO["Chromis katoi",  "Pacific, Northwest"] = 1
damsel_reg_tab_FAO["Chromis kennensis",  "Pacific, Southwest" ] = 1
damsel_reg_tab_FAO["Chrysiptera caesifrons", c("Pacific, Southwest", "Pacific, Western Central")] = 1
damsel_reg_tab_FAO["Chrysiptera ellenae", "Pacific, Western Central" ] = 1
damsel_reg_tab_FAO["Chrysiptera leucopoma", c("Pacific, Southwest", "Pacific, Western Central")] = 1
damsel_reg_tab_FAO["Chrysiptera maurineae", "Pacific, Western Central"] = 1
damsel_reg_tab_FAO["Chrysiptera papuensis", "Pacific, Western Central"] = 1
damsel_reg_tab_FAO["Dascyllus abudafur", c("Indian Ocean, Eastern", "Indian Ocean, Western")] = 1
damsel_reg_tab_FAO["Pomacentrus flavioculus", "Pacific, Western Central"] = 1

species_nreg <- apply(damsel_reg_tab_FAO, 1, sum, na.rm = T)# number of regions occpied by each species
missing_species <- names(species_nreg[species_nreg == 0]) # missing species regions data
missing_species # no missing data now

write.table(damsel_reg_tab_FAO, file = "./data/selected_species_FAOregions.txt", 
            quote = FALSE, sep = "\t")

# number of species endemic of only that one area
endemic_per_reg <- colSums(damsel_reg_tab[rowSums(damsel_reg_tab) == 1,])

# 1. remove areas with 0 endemism 
sel_regs <- names(endemic_per_reg)[endemic_per_reg > 0]
damsel_reg_tab.smpl <- damsel_reg_tab[, sel_regs]



# 2. simplify Atlantic, Indian and Pacific oceans
damsel_reg_tab_FAO.smpl <- damsel_reg_tab_FAO[,c("Asia - Inland waters", "Atlantic, Eastern Central", 
                                         "Indian Ocean, Eastern",  "Pacific, Eastern Central",
                                         "Oceania - Inland waters", "Mediterranean and Black Sea")]
colnames(damsel_reg_tab_FAO.smpl) <- c("Asia-inland", "Atlantic", "Indian", "Pacific", 
                                   "Oceania-inland", "Europe-inland")

damsel_reg_tab_FAO.smpl[,"Atlantic"] <- as.integer(rowSums(damsel_reg_tab_FAO[,grep("Atlantic", FAO_regions, value = TRUE)]) > 0)
damsel_reg_tab_FAO.smpl[,"Indian"] <- as.integer(rowSums(damsel_reg_tab_FAO[,grep("Indian", FAO_regions, value = TRUE)]) > 0)
damsel_reg_tab_FAO.smpl[,"Pacific"] <- as.integer(rowSums(damsel_reg_tab_FAO[,grep("Pacific", FAO_regions, value = TRUE)]) > 0)

write.table(damsel_reg_tab_FAO.smpl, file = "./data/selected_species_FAOregions_simlpified.txt", 
            quote = FALSE, sep = "\t")

#################################
# Using Cowman & Bellwood 2013 ##
#################################

# Cowman & Bellwood 2013 data
CB2013_data <- read.csv("./BioGeo/input/CB2013_pomacentridae_regions.csv", row.names = 1)
rownames(CB2013_data) <- sub(" ", "_", rownames(CB2013_data))
CB2013_species <- rownames(CB2013_data)
CB2013_species[!CB2013_species %in% selected_species] # in CB2013 not in ours
selected_species[!selected_species %in% CB2013_species] # in ours not in CB2013
colnames(CB2013_data) # CB2013 regions
# I decided to add Mediterranean and black sea as a new region for 2 reasons:
# FIRST: to see whether colonization of Mediterranean was posterior (through Atlantic Ocean) or early (through Indian Ocean)
# SECOND: to avoid Mediterranean species be classified as Atlantic in case some come from Indian Ocean or are actually not present in the Atlantic

CB2013_regions <- setNames(c("IO", "IAA", "CPO", "EPO", "AO", "MS"),
                    c("Indian", "Indo-Australian Archipelago", 
                      "Central Pacific", "East Pacific", "Atlantic", "Mediterranean"))

damsel_reg_tab_CB2013 <- data.frame(matrix(NA, nrow = length(pomacentridae_species),
                                           ncol = length(CB2013_regions), 
                                           dimnames = list(pomacentridae_species, 
                                                           CB2013_regions)))

# read meow ecoregions shape file to build a sf object with CB2013 regions
meow_regions <- sf::st_read("/Users/agarciaj/Unil/Research/1_NicheCompetition/snakemake_wd/data/MEOW_ECOS/", "meow_ecos")
plot(meow_regions[,"REALM"])
unique(meow_regions$REALM)
# subset targeted areas
meow_regions <- meow_regions[meow_regions$REALM %in% c("Western Indo-Pacific", "Central Indo-Pacific", "Eastern Indo-Pacific", 
                                       "Tropical Eastern Pacific", "Tropical Atlantic", "Temperate Southern Africa",
                                       "Temperate Australasia",  "Temperate South America" ,
                                       "Temperate Northern Pacific" , "Temperate Northern Atlantic"),]

# create CB2013 regions by renaming realms and provinces in a new column
meow_regions$CB2013[meow_regions$REALM == "Tropical Atlantic"] = "AO"
meow_regions$CB2013[meow_regions$REALM == "Western Indo-Pacific"] = "IO"
meow_regions$CB2013[meow_regions$REALM == "Central Indo-Pacific"] = "IAA"
meow_regions$CB2013[meow_regions$REALM == "Eastern Indo-Pacific"] = "CPO"
meow_regions$CB2013[meow_regions$REALM == "Tropical Eastern Pacific"] = "EPO"
plot(meow_regions[meow_regions$REALM == "Temperate Northern Atlantic","PROVINCE"])
unique(meow_regions$PROVINCE[meow_regions$REALM == "Temperate Northern Atlantic"])
meow_regions$CB2013[meow_regions$PROVINCE == "Lusitanian"] = "AO"
meow_regions$CB2013[meow_regions$PROVINCE == "Northern European Seas"] = "AO"
meow_regions$CB2013[meow_regions$PROVINCE == "Cold Temperate Northwest Atlantic"] = "AO"
meow_regions$CB2013[meow_regions$PROVINCE == "Warm Temperate Northwest Atlantic"] = "AO"
meow_regions$CB2013[meow_regions$PROVINCE == "Black Sea"] = "MS"
meow_regions$CB2013[meow_regions$PROVINCE == "Mediterranean Sea"] = "MS"
plot(meow_regions[meow_regions$REALM == "Temperate Northern Pacific","PROVINCE"])
unique(meow_regions$PROVINCE[meow_regions$REALM == "Temperate Northern Pacific"])
meow_regions$CB2013[meow_regions$PROVINCE == "Cold Temperate Northeast Pacific"] = "EPO"
meow_regions$CB2013[meow_regions$PROVINCE == "Warm Temperate Northeast Pacific"] = "EPO"
meow_regions$CB2013[meow_regions$PROVINCE == "Cold Temperate Northwest Pacific"] = "CPO"
meow_regions$CB2013[meow_regions$PROVINCE == "Warm Temperate Northwest Pacific"] = "CPO"
plot(meow_regions[meow_regions$REALM == "Temperate Southern Africa","PROVINCE"])
unique(meow_regions$PROVINCE[meow_regions$REALM == "Temperate Southern Africa"])
meow_regions$CB2013[meow_regions$PROVINCE == "Agulhas"] = "IO"
meow_regions$CB2013[meow_regions$PROVINCE == "Amsterdam-St Paul"] = "IO"
meow_regions$CB2013[meow_regions$PROVINCE == "Benguela"] = "AO"
plot(meow_regions[meow_regions$REALM == "Temperate South America","PROVINCE"])
unique(meow_regions$PROVINCE[meow_regions$REALM == "Temperate South America"])
meow_regions$CB2013[meow_regions$PROVINCE == "Warm Temperate Southeastern Pacific"] = "EPO"
meow_regions$CB2013[meow_regions$PROVINCE == "Juan Fernandez and Desventuradas"] = "EPO"
meow_regions$CB2013[meow_regions$PROVINCE == "Tristan Gough"] = "AO"
meow_regions$CB2013[meow_regions$PROVINCE == "Warm Temperate Southwestern Atlantic"] = "AO"
plot(meow_regions[meow_regions$PROVINCE == "Magellanic","ECOREGION"])
unique(meow_regions$ECOREGION[meow_regions$PROVINCE == "Magellanic"])
meow_regions$CB2013[meow_regions$ECOREGION %in% c("Chiloense", "Channels and Fjords of Southern Chile")] = "EPO"
meow_regions$CB2013[meow_regions$ECOREGION %in% c("North Patagonian Gulfs", "Malvinas/Falklands", "Patagonian Shelf")] = "AO"
plot(meow_regions[meow_regions$REALM == "Temperate Australasia","PROVINCE"])
unique(meow_regions$PROVINCE[meow_regions$REALM == "Temperate Australasia"])
meow_regions$CB2013[meow_regions$PROVINCE == "Southeast Australian Shelf" ] = "CPO"
meow_regions$CB2013[meow_regions$PROVINCE == "Southern New Zealand" ] = "CPO"
meow_regions$CB2013[meow_regions$PROVINCE == "Northern New Zealand" ] = "CPO"
meow_regions$CB2013[meow_regions$PROVINCE == "East Central Australian Shelf" ] = "CPO"
meow_regions$CB2013[meow_regions$PROVINCE == "West Central Australian Shelf" ] = "IO"
meow_regions$CB2013[meow_regions$PROVINCE == "Southwest Australian Shelf" ] = "IO"

# check result
plot(meow_regions[,"CB2013"])

# replace crs to lonlat 4326
meow_regions <- meow_regions[,"CB2013"] %>%
  st_transform(4326)
# convert regions to sp
meow_regions_sp <- as(meow_regions, "Spatial")
# save
sf::st_write(meow_regions, "./BioGeo/input/CB2013_regions_spdf.shp")
save(meow_regions_sp, file = "./BioGeo/input/CB2013_regions_spdf.rda")

### FILL CB2013 regions dataset 
#  1. From manually classified data
# load manually completed missing species regions data set
missing_species_occs <- read.csv("./BioGeo/input/missing_species_distribution_data.csv", row.names = 1)
missing_species_occs <- missing_species_occs[,1:5]
missing_species_occs$MS = 0
colnames(missing_species_occs) <- CB2013_regions

# transfer to main dataset
for (name in rownames(missing_species_occs)){
  damsel_reg_tab_CB2013[sub(" ", "_", name),] = missing_species_occs[name,]
}

# 2. From species occurrences
# load species occurrences
pomacentridae_occ <- read.csv("./BioGeo/input/pomacentridae_occurrences.csv")
sort(unique(pomacentridae_occ$species))
# convert occs to sp
pomacentridae_occ_sp <- pomacentridae_occ[,c(2:3, 1)]
coordinates(pomacentridae_occ_sp) <- ~lon+lat
crs(pomacentridae_occ_sp) <- crs(meow_regions_sp)

# check occurrences on CB2013 regions map
plot(meow_regions_sp, col = rep(viridis::viridis(6))[factor(meow_regions_sp$CB2013,
                                                           levels = CB2013_regions)], border = NA)
maps::map("world", add = T, col ="grey90", fill = T)
points(pomacentridae_occ_sp, pch = 21, bg = "red")

occs_CB2013_regions <- over(pomacentridae_occ_sp, meow_regions_sp)
pomacentridae_occ <- cbind(pomacentridae_occ, CB2013 = occs_CB2013_regions$CB2013)
pomacentridae_occ$CB2013 <- factor(pomacentridae_occ$CB2013, levels = CB2013_regions)

# check occurrences on CB2013 regions map
plot(meow_regions_sp, col = rep(viridis::viridis(6))[factor(meow_regions_sp$CB2013,
                                                            levels = CB2013_regions)], border = NA)
maps::map("world", add = T, col ="grey90", fill = T)
points(pomacentridae_occ[,2:3], pch = 21, 
       bg = rep(viridis::viridis(6))[pomacentridae_occ$CB2013])

# check if we may miss species occurrences with NA
unique(pomacentridae_occ$species[is.na(pomacentridae_occ$CB2013)]) %in% unique(pomacentridae_occ$species[!is.na(pomacentridae_occ$CB2013)])
# occs outside regions do not represent unique species distributions
pomacentridae_occ <- na.exclude(pomacentridae_occ) # remove NA regions

# check final result
plot(meow_regions_sp, col = rep(viridis::viridis(6))[factor(meow_regions_sp$CB2013,
                                                            levels = CB2013_regions)], border = NA)
maps::map("world", add = T, col ="grey90", fill = T)
points(pomacentridae_occ[,2:3], pch = 21, 
       bg = rep(viridis::viridis(6))[pomacentridae_occ$CB2013])

# get table of number of species observation per region
pomacentridae_nocc_tab <- table(pomacentridae_occ$species, pomacentridae_occ$CB2013)
# consider regions with less than 2.5% of occurrences per species as potential miss identifications or errors
pomacentridae_pocc_tab <- data.frame(matrix(as.integer(pomacentridae_nocc_tab/rowSums(pomacentridae_nocc_tab) > 0.025), 
                                 nrow = dim(pomacentridae_nocc_tab)[1],
                                 ncol = dim(pomacentridae_nocc_tab)[2], 
                                 dimnames = list(rownames(pomacentridae_nocc_tab), colnames(pomacentridae_nocc_tab))))

# CHECK MEDITERRANEAN SPECIES
# non-native (introduced) species in the Mediterranean: (Argyro Zenetos & Marika Galanidi 2020)
# 1. Abudefduf sexfasciatus - Giovos et al. 2018
# 2. Chrysiptera hemicyanea - Deidun et al. 2018
# 3. Abudefduf vaigiensis - Vella et al. 2016 / 2018-LY: Osca et al. 2020
pomacentridae_nocc_tab[pomacentridae_nocc_tab[,"MS"] > 1,]
pomacentridae_pocc_tab[pomacentridae_pocc_tab$MS == 1,]
# although Chrysiptera hemicyanea has been in the Mediterranean, does not occur in Mediterranean nor Atlantic based on fishbase 
# and it has been suggested to be an Aquarium release
pomacentridae_pocc_tab["Chrysiptera hemicyanea", c("AO", "MS")] = 0
# Chromis chromis seems the only species colonizing the Mediterranean naturally
####

# COMPARE with manual classification dataset
occ_miss_data <- pomacentridae_pocc_tab[rownames(pomacentridae_pocc_tab) %in% rownames(missing_species_occs),]
miss_occ_data <- missing_species_occs[rownames(missing_species_occs) %in% rownames(pomacentridae_pocc_tab),]
occ_miss_data[rowSums(occ_miss_data != miss_occ_data) > 0,]
miss_occ_data[rowSums(miss_occ_data != occ_miss_data) > 0,]

maps::map("world", col ="grey90", fill = T)
points(pomacentridae_occ[pomacentridae_occ$species == "Neopomacentrus taeniurus", 2:3], pch = 19)

# correct classifications
missing_species_occs["Azurina lepidolepis",] = pomacentridae_pocc_tab["Azurina lepidolepis",]
missing_species_occs["Neopomacentrus aktites",] = pomacentridae_pocc_tab["Neopomacentrus aktites",]
write.csv(missing_species_occs, "./data/missing_species_distribution_data_corrected.csv")

pomacentridae_pocc_tab["Neopomacentrus taeniurus",] = missing_species_occs["Neopomacentrus taeniurus",]
pomacentridae_pocc_tab["Pycnochromis acares",] = missing_species_occs["Pycnochromis acares",]
pomacentridae_pocc_tab["Pycnochromis alleni",] = missing_species_occs["Pycnochromis alleni",]
pomacentridae_pocc_tab["Pycnochromis iomelas",] = missing_species_occs["Pycnochromis iomelas",]

# LOOK FOR ODD ASSIGNMENTS
pomacentridae_pocc_tab[rowSums(pomacentridae_pocc_tab)> 2,]
# Neopomacentrus cyanomos does not occur in the Atlantic based on fishbase
pomacentridae_pocc_tab["Neopomacentrus cyanomos", "AO"] = 0
pomacentridae_pocc_tab[pomacentridae_pocc_tab$AO == 1,]

# transfer to main dataset
for (name in rownames(pomacentridae_pocc_tab)){
  damsel_reg_tab_CB2013[sub(" ", "_", name),] = pomacentridae_pocc_tab[name,]
}
rownames(na.exclude(damsel_reg_tab_CB2013))

# 3. From species distribution
# load species distributions
pomacentridae_dis <- read.csv("./BioGeo/input/pomacentridae_distributions.csv")
#pomacentridae_dis <- cbind(pomacentridae_dis[,1:2], pomacentridae_dis[,colnames(pomacentridae_dis) %in% selected_species])
pomacentridae_sdm <- NINA::raster_projection(pomacentridae_dis)
plot(sum(pomacentridae_sdm, na.rm = T))

pomacentridae_dis_locs <- tidyr::gather(pomacentridae_dis, species, presence, 3:ncol(pomacentridae_dis))
pomacentridae_dis_locs <- pomacentridae_dis_locs[pomacentridae_dis_locs$presence == 1,1:3]

# convert occs to sp
pomacentridae_dis_locs_sp <- pomacentridae_dis_locs
coordinates(pomacentridae_dis_locs_sp) <- ~Longitude+Latitude
crs(pomacentridae_dis_locs_sp) <- crs(meow_regions_sp)

occs_CB2013_regions <- over(pomacentridae_dis_locs_sp, meow_regions_sp)
pomacentridae_dis_locs <- cbind(pomacentridae_dis_locs, CB2013 = occs_CB2013_regions$CB2013)
pomacentridae_dis_locs$CB2013 <- factor(pomacentridae_dis_locs$CB2013, levels = CB2013_regions)

# check if we may miss species dis with NA
unique(pomacentridae_dis_locs$species[is.na(pomacentridae_dis_locs$CB2013)]) %in% unique(pomacentridae_dis_locs$species[!is.na(pomacentridae_dis_locs$CB2013)])
# occs outside regions do not represent unique species distributions
pomacentridae_dis_locs <- na.exclude(pomacentridae_dis_locs) # remove NA regions

# check final result
plot(meow_regions_sp, col = rep(viridis::viridis(6))[factor(meow_regions_sp$CB2013,
                                                            levels = CB2013_regions)], border = NA)
maps::map("world", add = T, col ="grey90", fill = T)
points(pomacentridae_dis_locs[,1:2], pch = 21, 
       bg = rep(viridis::viridis(6))[pomacentridae_dis_locs$CB2013])

# get table of number of species observation per region
pomacentridae_ndis_locs_tab <- table(pomacentridae_dis_locs$species, pomacentridae_dis_locs$CB2013)
# consider regions with less than 2.5% of occurrences per species as potential miss identifications or errors
pomacentridae_pdis_locs_tab <- data.frame(matrix(as.integer(pomacentridae_ndis_locs_tab/rowSums(pomacentridae_ndis_locs_tab) > 0.025), 
                                            nrow = dim(pomacentridae_ndis_locs_tab)[1],
                                            ncol = dim(pomacentridae_ndis_locs_tab)[2], 
                                            dimnames = list(rownames(pomacentridae_ndis_locs_tab), colnames(pomacentridae_ndis_locs_tab))))
# CHECK MEDITERRANEAN SPECIES
pomacentridae_ndis_locs_tab[pomacentridae_ndis_locs_tab[,"MS"] > 1,]
pomacentridae_pdis_locs_tab[pomacentridae_pdis_locs_tab$MS == 1,]
# Chromis pelloura is endemic from Red Sea (not in the Mediterranean)
pomacentridae_pdis_locs_tab["Chromis_pelloura", "MS"] = 0
###

colnames(pomacentridae_dis_locs)[1:2] <- c("lon", "lat")
pomacentridae_occ_dis_CB2013_regs <- rbind(pomacentridae_dis_locs, pomacentridae_occ[,c(2:3, 1, 4)])
write.csv(pomacentridae_occ_dis_CB2013_regs, "./BioGeo/input/pomacentridae_dis_occ_CB2013_regs.csv", row.names = F)

# COMPARE WITH OCCURRENCE DATASET
occ_dis_data <- pomacentridae_pocc_tab[sub(" ", "_", rownames(pomacentridae_pocc_tab)) %in% rownames(pomacentridae_ndis_locs_tab),]
dis_occ_data <- pomacentridae_pdis_locs_tab[rownames(pomacentridae_pdis_locs_tab) %in% sub(" ", "_", rownames(pomacentridae_pocc_tab)) ,]

# Conservative approach correcting with the data set with less number of regions in case of mismatch between datasets
occ_dis_data_mismatch <- occ_dis_data[rowSums(occ_dis_data != dis_occ_data) > 0,]
dis_occ_data_mismatch <- dis_occ_data[rowSums(dis_occ_data != occ_dis_data) > 0,]

# get differences in number of regions between occurrences and distributions data sets
occ_dis_nreg_diff <- rowSums(occ_dis_data_mismatch) - rowSums(dis_occ_data_mismatch)
pomacentridae_pocc_tab[rownames(occ_dis_data_mismatch[occ_dis_nreg_diff > 0,]),] = dis_occ_data_mismatch[occ_dis_nreg_diff > 0,]
pomacentridae_pdis_locs_tab[rownames(dis_occ_data_mismatch[occ_dis_nreg_diff < 0,]),] = occ_dis_data_mismatch[occ_dis_nreg_diff < 0,] 

# check same number of regions but different assignments between data sets
occ_dis_data_mismatch[occ_dis_nreg_diff == 0,]
dis_occ_data_mismatch[occ_dis_nreg_diff == 0,]

plot(pomacentridae_sdm$Chromis_xanthochira)
pomacentridae_pdis_locs_tab["Chromis_xanthochira",] = occ_dis_data_mismatch["Chromis xanthochira",]
plot(pomacentridae_sdm$Chromis_xouthos)
pomacentridae_pocc_tab["Chromis xouthos",] = dis_occ_data_mismatch["Chromis_xouthos",]

# LOOK FOR ODD ASSIGNMENTS
pomacentridae_pdis_locs_tab[rowSums(pomacentridae_pdis_locs_tab)> 2,]
rownames(pomacentridae_pdis_locs_tab[pomacentridae_pdis_locs_tab$AO == 1,])
rownames(pomacentridae_pdis_locs_tab[pomacentridae_pdis_locs_tab$EPO == 1,])
rownames(pomacentridae_pdis_locs_tab[pomacentridae_pdis_locs_tab$IO == 1,])
rownames(pomacentridae_pdis_locs_tab[pomacentridae_pdis_locs_tab$IAA == 1,])
rownames(pomacentridae_pdis_locs_tab[pomacentridae_pdis_locs_tab$CPO == 1,])
# Chromis enchrysura does not occur in the Eastern Pacific
pomacentridae_pdis_locs_tab["Chromis_enchrysura", "EPO"] = 0
# Chromis intercrusma does not occur in the Atlantic
pomacentridae_pdis_locs_tab["Chromis_intercrusma", "AO"] = 0

# transfer to main dataset
for (name in rownames(pomacentridae_pdis_locs_tab)){
  damsel_reg_tab_CB2013[name,] = pomacentridae_pdis_locs_tab[name,]
}

missing_reg_data <- rownames(damsel_reg_tab_CB2013[rowSums(is.na(damsel_reg_tab_CB2013)) > 0,])
missing_reg_data[missing_reg_data %in% selected_species] # missing data in the species selection
# Stegastes_baldwini occurs in the East Pacific (Clipperton Island)
damsel_reg_tab_CB2013["Stegastes_baldwini",] <- c(0,0,0,1,0,0)
damsel_reg_tab_CB2013["Altrichthys_alelia",] <- c(0,1,0,0,0,0)
damsel_reg_tab_CB2013["Amblypomacentrus_kuiteri",] <- c(0,1,0,0,0,0)
damsel_reg_tab_CB2013["Chromis_hangganan",] <- c(0,1,0,0,0,0)
damsel_reg_tab_CB2013["Chrysiptera_uswanasi",] <- c(0,1,0,0,0,0)
damsel_reg_tab_CB2013["Neopomacentrus_aquadulcis",] <- c(0,1,0,0,0,0)

#write.csv(damsel_reg_tab_CB2013, "./BioGeo/input/CB2013_regions_all_pomacentridae.csv")
damsel_reg_tab_CB2013 <- read.csv("./BioGeo/input/CB2013_regions_all_pomacentridae.csv", row.names = 1)
# remove species with missing data
damsel_reg_tab_CB2013 <- na.exclude(damsel_reg_tab_CB2013)
nrow(damsel_reg_tab_CB2013)

# check missing species
selected_species[!selected_species %in% rownames(damsel_reg_tab_CB2013)] 

# subset to species selection for analyses
damsel_reg_tab_CB2013_subset <- damsel_reg_tab_CB2013[selected_species,]
write.csv(damsel_reg_tab_CB2013_subset, "./BioGeo/input/CB2013_regions_subset_pomacentridae.csv")

# replace Mediterranean Region (only Chromis Chromis natively present) for Tethys for BioGeoBEARS modelling
colnames(damsel_reg_tab_CB2013_subset)[6] <- "TS"
CB2013_regions <- setNames(c("IO", "IAA", "CPO", "EPO", "AO", "TS"),
                           c("Indian", "Indo-Australian Archipelago", 
                             "Central Pacific", "East Pacific", "Atlantic", "Tethys"))
# remove occurrences in Tethys
damsel_reg_tab_CB2013_subset$TS = 0

################################

################################
## Preparing BioGeoBEARS data ##
################################

#x <- write_BioGeoBEARS_geodata(damsel_reg_tab, filename = "./BioGeo/input/pomacentridae_FAO")
#x <- write_BioGeoBEARS_geodata(damsel_reg_tab.smpl, filename = "./BioGeo/input/pomacentridae_FAOsimplified")
x <- write_BioGeoBEARS_geodata(damsel_reg_tab_CB2013_subset, filename = "./BioGeo/input/CB2021_pomacentridae")
x$description

damsel_data <- BioGeoBEARS::getranges_from_LagrangePHYLIP("./BioGeo/input/pomacentridae_FAO_geodata_input.txt")
damsel_data <- BioGeoBEARS::getranges_from_LagrangePHYLIP("./BioGeo/input/pomacentridae_FAO_simplified_geodata_input.txt")
damsel_data <- BioGeoBEARS::getranges_from_LagrangePHYLIP("./BioGeo/input/CB2021_pomacentridae_geodata_input.txt")

species_nreg <- apply(damsel_reg_tab, 1, sum, na.rm = T)
species_nreg <- apply(damsel_reg_tab.smpl, 1, sum, na.rm = T)
species_nreg <- apply(damsel_reg_tab_CB2013_subset, 1, sum, na.rm = T)

damsel_tree_subset <- read.tree("./data/damsel_subset.tre")
plot(damsel_tree_subset)
damsel_tree_subset_rootedge0 <- damsel_tree_subset
damsel_tree_subset_rootedge0$root.edge = NULL
plot(damsel_tree_subset_rootedge0)
write.tree(damsel_tree_subset_rootedge0, "./BioGeo/input/damsel_tree_subset_norootedge.tre")

tmp <- damsel_data@df
tmp[,1:ncol(tmp)] <- lapply(tmp[,1:ncol(tmp)], factor)
colnames(tmp) <- x$description$area
colors <- setNames(replicate(ncol(tmp), setNames(c("white", "darkgray"), 0:1), simplify = FALSE), 
                   colnames(tmp))
dev.off()
png("./BioGeo/DamselTree_CB2013_areas.png", width = 1200, height = 3200)
plotTree.datamatrix(damsel_tree_subset, tmp, fsize = 0.8, sep=0.2, srt=90, yexp=1, xexp=1.1,
                    space=0, header = TRUE, colors = colors)
legend("topright", c("absent", "present"), 
       pch = 22, pt.bg = c("white", "darkgray"), pt.cex = 2,
       cex = 1.8, bty = "n")
dev.off()

max_range_size = max(species_nreg)

bgb_run <- define_BioGeoBEARS_run(
  max_range_size = max_range_size,
  trfn = "/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/damsel_tree_subset_norootedge.tre",
  return_condlikes_table = TRUE)

bgb_run$geogfn <- "/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/CB2021_pomacentridae_geodata_input.txt"

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=bgb_run$geogfn)
tipranges
nrow(tipranges@df) # check
# check number of states to approximate runtime
# less than 1000 if you want the analysis to run in under a day
# less than 1500 if you want the analysis to run in under a week, and
# less than 2500 if you want it to run at all
nstates <- cladoRcpp::numstates_from_numareas(numareas=length(CB2013_regions), 
                                              maxareas=max_range_size, 
                                              include_null_range=TRUE)
print(nstates)

bgb_run$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
bgb_run$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Speed options and multicore processing if desired
bgb_run$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
bgb_run$speedup = TRUE          # shortcuts to speed ML search; use FALSE if worried (e.g. >3 params)
bgb_run$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();
# if "GenSA", use Generalized Simulated Annealing, which seems better on high-dimensional
# problems (5+ parameters), but seems to sometimes fail to optimize on simple problems
bgb_run$num_cores_to_use = 4 # (use more cores to speed it up; this requires library(parallel) and/or library(snow).

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
bgb_run$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# Set up a time-stratified analysis:
stratified = TRUE
select_states = FALSE
allowed_areas = TRUE
distance_matrix = TRUE
adjacency_matrix = FALSE

#bgb_run$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
if (allowed_areas){ bgb_run$areas_allowed_fn = "/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/damsel_biogeo_allwdmatrix_Descombes.txt" }
if (distance_matrix){ bgb_run$areas_adjacency_fn = "/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/damsel_biogeo_adjmatrix_Descombes.txt" }
if (adjacency_matrix){ bgb_run$distsfn = "/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/damsel_biogeo_distmatrix_Descombes.txt" }
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

if (stratified){
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  bgb_run$timesfn <- "/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/damsel_biogeo_timeperiods_Descombes.txt"
  bgb_run <-  readfiles_BioGeoBEARS_run(bgb_run)
  bgb_run = section_the_tree(inputs=bgb_run, make_master_table=TRUE, 
                             plot_pieces=FALSE, 
                             fossils_older_than=0.001, 
                             cut_fossils=FALSE)
  
  # The stratified tree is described in this table:
  bgb_run$master_table
  bgb_run$timeperiods

  # INPUT the NEW states list into the BioGeoBEARS_run_object
  if (select_states){
    selected_states <- readRDS("/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/CB2013_biogeo_selected_states_Descombes.rds")
    unique_states <- readRDS("./BioGeo/input/CB2013_biogeo_unique_states_Descombes.rds")
    if (is.list(last(selected_states))){
      for (i in 1:length(selected_states)){
        bgb_run$lists_of_states_lists_0based[[i]] = selected_states[[i]]
      }
      #bgb_run$states_list = unique_states
    } else {
      bgb_run$states_list = selected_states
    }
  }
} else {
  # INPUT the NEW states list into the BioGeoBEARS_run_object
  if (select_states){
    selected_states <- readRDS("/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/CB2013_biogeo_selected_states_Descombes.rds")
    bgb_run$states_list = selected_states
  }
}

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
bgb_run <-  readfiles_BioGeoBEARS_run(bgb_run)

# Good default settings to get ancestral states
bgb_run$return_condlikes_table = TRUE
bgb_run$calc_TTL_loglike_from_condlikes_table = TRUE
bgb_run$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# Look at the model object
bgb_run$BioGeoBEARS_model_object

# Run this to check inputs. Read the error messages if you get them!
bgb_run = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=bgb_run)
check_BioGeoBEARS_run(bgb_run)

model_specifications <- "_CB2013_regs_strat_allwdareas_shortestdistmat_Descombes"
models_path <- "./BioGeo/Result/models"

DEC.fit <- NULL
DEC_J.fit <- NULL
DIVAlike.fit <- NULL
DIVAlike_J.fit <- NULL
BAYAREAlike.fit <- NULL
BAYAREAlike_J.fit <- NULL

################################################################################
# 12.1.1. The Dispersal-Extinction-Cladogenesis model (DEC)
################################################################################
##################################
# --------- DEC model -----------#
##################################
DEC.fit <- bears_optim_run(bgb_run)

DEC.fit$optim_result

saveRDS(DEC.fit, file.path(models_path, paste0("DEC",  model_specifications, ".rds")))

layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
par(mar = c(4.1, 0.1, 3.1, 0.1))
plot_BioGeoBEARS_results(DEC.fit, analysis_titletxt = "DEC model", plotlegend = TRUE, 
                         addl_params=list("j"), plotwhat="text", label.offset=0.45, 
                         tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, 
                         plotsplits=TRUE, include_null_range=TRUE, 
                         tr=damsel_tree_subset_rootedge0, tipranges=tipranges)

# Pie chart
dev.off()
plot_BioGeoBEARS_results(DEC.fit, analysis_titletxt = "DEC model", addl_params=list("j"), 
                         plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, 
                         splitcex=0.6, titlecex=0.8, plotsplits=TRUE, 
                         include_null_range=TRUE, tr=damsel_tree_subset_rootedge0, tipranges=tipranges)

##################################
# ------- DEC+J model -----------#
##################################

dstart <- DEC.fit$outputs@params_table["d", "est"]
estart <- DEC.fit$outputs@params_table["e", "est"]
jstart <- 0.0001

bgb_run$BioGeoBEARS_model_object@params_table["d", "init"] <- dstart
bgb_run$BioGeoBEARS_model_object@params_table["d", "est"] <- dstart
bgb_run$BioGeoBEARS_model_object@params_table["e", "init"] <- estart
bgb_run$BioGeoBEARS_model_object@params_table["e", "est"] <- estart
bgb_run$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
bgb_run$BioGeoBEARS_model_object@params_table["j", "init"] <- jstart
bgb_run$BioGeoBEARS_model_object@params_table["j", "est"] <- jstart

check_BioGeoBEARS_run(bgb_run)

DEC_J.fit <- bears_optim_run(bgb_run)

DEC_J.fit$optim_result

saveRDS(DEC_J.fit, file.path(models_path, paste0("DECJ",  model_specifications, ".rds")))

layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
par(mar = c(4.1, 0.1, 3.1, 0.1))
plot_BioGeoBEARS_results(DEC_J.fit, analysis_titletxt = "DEC+J model",
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

################################################################################
# 12.1.2. The Dispersal-Vicariance model (DIVA-like)
################################################################################
##################################
# ------- DIVA-like model -------#
##################################
# Remove subset-sympatry
bgb_run$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
bgb_run$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
bgb_run$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

bgb_run$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
bgb_run$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
bgb_run$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
bgb_run$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
bgb_run$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
bgb_run$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
bgb_run$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

bgb_run = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=bgb_run)
check_BioGeoBEARS_run(bgb_run)

DIVAlike.fit <- bears_optim_run(bgb_run)

DIVAlike.fit$optim_result

saveRDS(DIVAlike.fit, file.path(models_path, paste0("DIVAlike",  model_specifications, ".rds")))

layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
par(mar = c(4.1, 0.1, 3.1, 0.1))
plot_BioGeoBEARS_results(DIVAlike.fit, analysis_titletxt = "DIVA-like model",
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

##################################
# ------ DIVA-like+J model ------#
##################################
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = DIVAlike.fit$outputs@params_table["d","est"]
estart = DIVAlike.fit$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
bgb_run$BioGeoBEARS_model_object@params_table["d","init"] = dstart
bgb_run$BioGeoBEARS_model_object@params_table["d","est"] = dstart
bgb_run$BioGeoBEARS_model_object@params_table["e","init"] = estart
bgb_run$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add jump dispersal/founder-event speciation
bgb_run$BioGeoBEARS_model_object@params_table["j","type"] = "free"
bgb_run$BioGeoBEARS_model_object@params_table["j","init"] = jstart
bgb_run$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
bgb_run$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
bgb_run$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

bgb_run = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=bgb_run)
check_BioGeoBEARS_run(bgb_run)

DIVAlike_J.fit <- bears_optim_run(bgb_run)

DIVAlike_J.fit$optim_result

saveRDS(DIVAlike_J.fit, file.path(models_path, paste0("DIVAlikeJ",  model_specifications, ".rds")))

layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
par(mar = c(4.1, 0.1, 3.1, 0.1))
plot_BioGeoBEARS_results(DIVAlike_J.fit, analysis_titletxt = "DIVA-like+J model",
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

################################################################################
# 12.1.2. The Bayesian Area model (BAYAREA-like)
################################################################################
##################################
# ----- BAYAREA-like model ------#
##################################
# No subset sympatry
bgb_run$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
bgb_run$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
bgb_run$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
bgb_run$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
bgb_run$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
bgb_run$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
bgb_run$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
bgb_run$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
bgb_run$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
bgb_run$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
bgb_run$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
bgb_run$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs; fixing any initial ("init") values outside min/max
bgb_run = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=bgb_run)
check_BioGeoBEARS_run(bgb_run)

BAYAREAlike.fit <- bears_optim_run(bgb_run)

BAYAREAlike.fit$optim_result

saveRDS(BAYAREAlike.fit, file.path(models_path, paste0("BAYAREAlike",  model_specifications, ".rds")))

layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
par(mar = c(4.1, 0.1, 3.1, 0.1))
plot_BioGeoBEARS_results(BAYAREAlike.fit, analysis_titletxt = "DIVA-like+J model",
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

##################################
# ---- BAYAREA-like+J model -----#
##################################

dstart = BAYAREAlike.fit$outputs@params_table["d","est"]
estart = BAYAREAlike.fit$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
bgb_run$BioGeoBEARS_model_object@params_table["d","init"] = dstart
bgb_run$BioGeoBEARS_model_object@params_table["d","est"] = dstart
bgb_run$BioGeoBEARS_model_object@params_table["e","init"] = estart
bgb_run$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
bgb_run$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
bgb_run$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
bgb_run$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
bgb_run$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
bgb_run$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
bgb_run$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
bgb_run$BioGeoBEARS_model_object@params_table["j","type"] = "free"
bgb_run$BioGeoBEARS_model_object@params_table["j","init"] = jstart
bgb_run$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
bgb_run$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
bgb_run$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
bgb_run$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
bgb_run$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
bgb_run$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
bgb_run$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
bgb_run$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
# machines. I can't replicate this on my Mac machines, but it is almost certainly
# just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to
# prevent this, but apparently optim/optimx sometimes go slightly beyond 
# these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
# slightly for each parameter:
#bgb_run$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
#bgb_run$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
#bgb_run$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
#bgb_run$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
#bgb_run$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
#bgb_run$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

bgb_run = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=bgb_run)
check_BioGeoBEARS_run(bgb_run)

BAYAREAlike_J.fit <- bears_optim_run(bgb_run)

BAYAREAlike_J.fit$optim_result

saveRDS(BAYAREAlike_J.fit, file.path(models_path, paste0("BAYAREAlikeJ",  model_specifications, ".rds")))

layout(matrix(1:2, 1, 2), widths = c(0.2, 0.8))
par(mar = c(4.1, 0.1, 3.1, 0.1))
plot_BioGeoBEARS_results(BAYAREAlike_J.fit, analysis_titletxt = "DIVA-like+J model",
                         plotlegend = TRUE, tipcex = 0.4, statecex = 0.4)

##################################
# -------- Check results --------#
##################################
null_models <- c("DEC", "DIVAlike", "BAYAREAlike")

models <- list(DEC = DEC.fit, 
               DEC_J = DEC_J.fit, 
               DIVAlike = DIVAlike.fit, 
               DIVAlike_J = DIVAlike_J.fit, 
               BAYAREAlike = BAYAREAlike.fit, 
               BAYAREAlike_J = BAYAREAlike_J.fit)
models <- models[sapply(models, length) > 0]

null_models <- null_models[null_models %in% names(models)]

logL_models <- sapply(models, get_LnL_from_BioGeoBEARS_results_object)
names(logL_models) <- names(models)
logL_models <- logL_models[colSums(sapply(models, sapply, length)) > 0]

AIC.table <- do.call(rbind, sapply(null_models, function(m) 
  conditional_format_table(AICstats_2models(logL_models[paste0(m, "_J")], 
                                            logL_models[m], 
                                            numparams1 = 3, 
                                            numparams2 = 2)), simplify = F))
AIC.table$alt = paste0(null_models, "_J")
AIC.table$null = null_models

print(AIC.table)

# DEC, null model for Likelihood Ratio Test (LRT)
res_table <- do.call(rbind, lapply(models, function(m) 
  extract_params_from_BioGeoBEARS_results_object(results_object=m, 
                                                 returnwhat="table", 
                                                 addl_params=c("j"), 
                                                 paramsstr_digits=4))) 

#######################################################
# Save the results tables for later -- check for e.g.
# convergence issues
#######################################################

# Loads to "restable"
saveRDS(res_table, file = paste0("./BioGeo/damsel_biogeo_restable",  model_specifications,".rds"))
res_table <- readRDS(file="./BioGeo/damsel_biogeo_restable_CB2013_regs_strat_allwdareas_shortestdistmat_Descombes.rds")

# Loads to "teststable"
saveRDS(AIC.table, file = paste0("./BioGeo/damsel_biogeo_testtable",  model_specifications,".rds"))
AIC.table <- readRDS(file="./BioGeo/damsel_biogeo_testable_CB2013_regs_strat_allwdareas_shortestdistmat_Descombes.rds")

# Also save to text files
write.table(unlist_df(AIC.table), file = paste0("./BioGeo/damsel_biogeo_testtable",  model_specifications,".txt"), quote=FALSE, sep="\t")

#######################################################
# Model weights of all six models
#######################################################

# With AICs:
AICtable = calc_AIC_column(LnL_vals=res_table$LnL, nparam_vals=res_table$numparams)
res_table = cbind(res_table, AICtable)
res_table = AkaikeWeights_on_summary_table(restable=res_table, colname_to_use="AIC")

# With AICcs -- factors in sample size
samplesize = length(damsel_tree_subset$tip.label)
AICtable = calc_AICc_column(LnL_vals=res_table$LnL, nparam_vals=res_table$numparams, samplesize=samplesize)
res_table = cbind(res_table, AICtable)
res_table = AkaikeWeights_on_summary_table(restable=res_table, colname_to_use="AICc")

# Save with nice conditional formatting
write.table(res_table, file = paste0("./BioGeo/damsel_biogeo_restable",  model_specifications,".txt"), quote=FALSE, sep="\t")

best_model_name <- rownames(res_table)[which.min(res_table$AICc)]
res <- models[[best_model_name]]

################################
## LOAD BEST BIOGEOGRAPHIC MODEL
################################
best_model_name <- "DEC_J"
res <- readRDS(file.path(models_path, paste0("DECJ",  model_specifications, ".rds")))

# Best model result
png(file.path(models_path, paste0("DECJ",  model_specifications, "_bestmodel.png")), width = 4000, height = 6000)
layout(matrix(c(3, 1, 2, 1, 3, 3, 1, 1, 1, 3, 3, 1, 1, 1, 3, 3, 1, 1, 1, 3), 
              nrow = 5, ncol = 4), heights = c(0.05, 0.1, 0.4, 0.4, 0.05))
par(mar = c(5.1, 4.1, 3.1, 0.1),  mgp=c(3, 5, 0))
my_plot_BioGeoBEARS_results(res, analysis_titletxt = paste(sub("_", "+", best_model_name), "model"),
                         plotlegend = FALSE, tipcex = 2, statecex = 2, titlecex = 10, 
                         cex.axis = 8, cex.lab = 8, axis.lwd = 5, title_offset = -15)
make_legend(res, states = 1:num_areas, names = names(CB2013_regions), cex = 6,  blank = T)

dev.off()


#######################################################
# Stochastic mapping on DEC non-stratified, non_adjacent and distance matrix
#######################################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
save(stochastic_mapping_inputs_list, file= paste0("./BioGeo/Result/Best_model_BSM/BSM_pomacentridae_", best_model_name, model_specifications, "_inputList.Rdata"))

# Run BSM
BSM_output = runBSM(res, 
                    stochastic_mapping_inputs_list = stochastic_mapping_inputs_list, 
                    maxnum_maps_to_try = 100, 
                    nummaps_goal = 50, 
                    maxtries_per_branch=40000, 
                    save_after_every_try=TRUE, 
                    savedir=getwd(), 
                    seedval=12345, 
                    wait_before_save=0.01, 
                    master_nodenum_toPrint=0)

saveRDS(BSM_output, file = paste0("./BioGeo/Result/Best_model_BSM/BSM_pomacentridae_", best_model_name, model_specifications, ".rds"))

best_model_name <- "DEC_J"
BSM_output <- readRDS(paste0("./BioGeo/Result/Best_model_BSM/BSM_pomacentridae_", best_model_name, model_specifications, ".rds"))

RES_clado_events_tables = BSM_output$RES_clado_events_tables
RES_ana_events_tables = BSM_output$RES_ana_events_tables

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

include_null_range = T
areanames = names(tipranges@df)
areas = areanames
max_range_size = max(species_nreg)

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
#states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

selected_states_0 <- list()
states <- NULL
for (i in 1:length(selected_states)){
  s <- sapply(selected_states[[1]], paste, collapse = "")
  selected_states_0 <- c(selected_states_0, selected_states[[1]][which(!s %in% states)])
  states <- c(states, s[which(!s %in% states)])
}

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, 
                                                           states_list_0based=selected_states_0, 
                                                           max_range_size=max_range_size, 
                                                           plot_null_range=TRUE)

############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

# cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
# colnums = match(cols_to_get, names(ana_events_table))
# ana_events_table_cols_to_add = ana_events_table[,colnums]
# anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
# ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
# rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
# master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

############################################
# Open a PDF
############################################
pdffn = paste0("./BioGeo/Result/", best_model_name, model_specifications, "_single_stochastic_map_n1.pdf")
pdf(file=pdffn, height=50, width=25)
stratified = TRUE
# Convert the BSM into a modified res object
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, 
                                        master_table_cladogenetic_events=master_table_cladogenetic_events, 
                                        stratified=stratified)

resmod$inputs$states_list <- selected_states_0

my_plot_BioGeoBEARS_results(results_object=resmod, 
                         analysis_titletxt="Stochastic map", 
                         addl_params=list("j"), 
                         label.offset=0.5, 
                         plotwhat="text", 
                         cornercoords_loc=scriptdir, 
                         root.edge=TRUE, 
                         colors_list_for_states=colors_list_for_states, 
                         skiptree=FALSE, titlecex = 2, cex.axis = 1, axis.lwd = 2,
                         show.tip.label=TRUE)

# Paint on the branch states
paint_stochastic_map_branches(res=resmod, 
                              master_table_cladogenetic_events=master_table_cladogenetic_events, 
                              colors_list_for_states=colors_list_for_states, 
                              lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)

my_plot_BioGeoBEARS_results(results_object=resmod, 
                         analysis_titletxt="Stochastic map", 
                         addl_params=list("j"), plotwhat="text", 
                         cornercoords_loc=scriptdir, root.edge=TRUE, 
                         colors_list_for_states=colors_list_for_states, titlecex = 2, cex.axis = 1, axis.lwd = 2,
                         skiptree=TRUE, show.tip.label=TRUE)

############################################
# Close PDF
############################################
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Plot all 50 stochastic maps to PDF
#######################################################
# Setup
include_null_range = include_null_range
areanames = areanames
areas = areanames
max_range_size = max_range_size
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=selected_states, 
                                                           max_range_size=max_range_size, plot_null_range=TRUE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# Loop through the maps and plot to PDF
pdffn = paste0("./BioGeo/Result/", best_model_name, "_", length(clado_events_tables), model_specifications, "BSMs_v1.pdf")
pdf(file=pdffn, height=45, width=15)

res$inputs$states_list <- selected_states_0

nummaps_goal = 50
for (i in 1:nummaps_goal) {
  clado_events_table = clado_events_tables[[i]]
  analysis_titletxt = paste0(best_model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  plot_BSM(results_object=res, clado_events_table=clado_events_table, 
           stratified=stratified, analysis_titletxt=analysis_titletxt, 
           addl_params=list("j"), label.offset=0.5, plotwhat="text", 
           cornercoords_loc=scriptdir, root.edge=TRUE, 
           colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, 
           include_null_range=include_null_range)
} # END for (i in 1:nummaps_goal)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Summarize stochastic map tables
#######################################################
dir_path <- paste0("./BioGeo/Result/Best_model_BSM/BSM_pomacentridae_" , best_model_name, model_specifications, "RES")
dir.create(dir_path)

length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = area_code[areanames]
actual_names

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

summary_counts_BSMs = counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Histogram of event counts
hist_event_counts(counts_list, pdffn=paste0(dir_path, "/", best_model_name, "_histograms_of_event_counts.pdf"))

#######################################################
# Print counts to files
#######################################################

tmpnames = names(counts_list)
cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
for (i in 1:length(tmpnames)){
  cmdtxt = paste0("item = counts_list$", tmpnames[i])
  eval(parse(text=cmdtxt))
  
  # Skip cubes
  if (length(dim(item)) != 2)
  {
    next()
  }
  
  outfn = paste0(tmpnames[i], ".txt")
  if (length(item) == 0)
  {
    cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
    cat("\n")
  } else {
    cat(outfn)
    cat("\n")
    write.table(conditional_format_table(item), file=file.path(dir_path, outfn), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  } # END if (length(item) == 0)
} # END for (i in 1:length(tmpnames))
cat("...done.\n")

#######################################################
# Check that ML ancestral state/range probabilities and
# the mean of the BSMs approximately line up
#######################################################
# For 95% CIs on BSM counts
check_ML_vs_BSM(res, clado_events_tables, best_model_name, tr=NULL, 
                plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)


results_path <- "./BioGeo/Result/Best_model_BSM/BSM_pomacentridae_DEC_J_CB2013_regs_strat_allwdareas_shortestdistmat_DescombesRES.rds"
tree_path = "/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/damsel_tree_subset_norootedge.tre"
geo_data_path <- "/Users/agarciaj/Unil/Research/5_Damselfish_evo/BioGeo/input/CB2021_pomacentridae_geodata_input.txt"
area_names <- c("IO", "IAA", "CPO", "EPO", "AO", "TS")

revbgb <- bgb_to_revgadgets(results_path, geo_data_path, tree_path, area_names = area_names) 
  
color_codes <- setNames(c("#D2B83F", "#7994CF", "#D5DDDC", "#332748", "#DE3922", "#4B6356"), area_names)
all_states <- unique(revbgb@data$end_state_1)

st.list <- sapply(names(color_codes), function(i) {
  grep(i,all_states)
})
arealist <- sapply(1:length(all_states), function(i) {
  names(st.list)[sapply(st.list, function(x) i %in% x)]
})
all_states <- setNames(sapply(arealist, paste, collapse = "\n"), all_states)
st.list <- sapply(1:length(all_states), function(i) {
  mix.colors(color_codes[names(st.list)[sapply(st.list, function(x) i %in% x)]])
})
names(st.list) <- names(all_states)

# simplify states for plotting
mix_states <- setdiff(names(all_states), area_names)
states <- setNames(revbgb@data$end_state_1[(Ntip(revbgb@phylo)+1):nrow(revbgb@data)], 
                   as.character(Ntip(revbgb@phylo)+1):nrow(revbgb@data))
transitions <- get_state_transitions(revbgb@phylo, states)
state_matrix <- to.matrix(states,
                          names(st.list))

simplified_state_matrix <- matrix(0, nrow = nrow(state_matrix), ncol = length(area_names), 
                                  dimnames = list(rownames(state_matrix), area_names))
for (i in 1:nrow(state_matrix)) {
  
  v <- state_matrix[i,]
  v_state <- names(v)[which(v == 1)]
  if(!v_state  %in% area_names) {
    v_main_states <- setNames(as.integer(sapply(area_names, function(a) grepl(a, v_state))), area_names)
  } else {
    v_main_states <- v[area_names]
    v_main_states[is.na(v_main_states)] = 0
    names(v_main_states) <- area_names
  }
  simplified_state_matrix[i, ] = v_main_states
}

# make legend
legend <- c("IO:     Indian Ocean", "IAA:   Indo-Australian Archipelago",
            "CPO: Central Pacific Ocean", "EPO: Eastern Pacific Ocean", 
            "AO:   Atlantic Ocean", "TS:    Tethys Sea")
legend <- sapply(area_names, function(i) grep(i, legend, value = TRUE))

# get mix states text
text.reg <- revbgb@data$end_state_1[(Ntip(revbgb@phylo)+1):nrow(revbgb@data)]
text.reg <- all_states[text.reg]
text.reg[text.reg %in% area_names] = ""

## PLOT
pdf("./new_figures/BSR/damsel_biogeographic_reconstruction.pdf", width = 12, height = 12)
phytools::plotTree(revbgb@phylo,show.tip.label = TRUE, type = "fan", part = 0.88,
         fsize = 0.5, offset = 12, lwd = 3)
xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(revbgb@phylo)]
yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(revbgb@phylo)]

# add ages lines
obj<-axis(1,pos=-2,at=seq(10,50,by=10),cex.axis=0.5,labels=FALSE)
text(obj,rep(-5,length(obj)), rev(obj),cex=1.1)
text(mean(obj),-8,"time (MYA)",cex=1.2)
for(i in 1:(length(obj)-1)){
  a1<-atan(-2/obj[i])
  a2<-0.88*2*pi
  plotrix::draw.arc(0,0,radius=obj[i],a1,a2,lwd=1.5,
                    col=make.transparent("grey40",0.25))
}
## add traits
points(xx*1.035, yy*1.035, pch = 21, cex = 1,
       bg = st.list[revbgb@data$end_state_1[1:Ntip(revbgb@phylo)]])
# add legend
legend(x="topleft", legend=legend,
       pt.cex=2, pch=16, col=color_codes, title = "", bty = "n")
# Add state-colored nodes
nodelabels(frame = "none", cex = ifelse(transitions, 0.35, 0),
  pie = simplified_state_matrix, piecol = color_codes)
dev.off()
