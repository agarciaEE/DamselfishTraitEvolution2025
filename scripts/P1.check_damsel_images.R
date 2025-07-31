## Settings
###########
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/")

# load custom functions
source("./Rscripts/custom_functions.R")

all_damsel_pics <- list.files("./images/1st_selection_images/")
actual_damsel_pics <- list.files("./images/3rd_selection_landmarks/")
actual_damsel_pics <- sub("-", " ",actual_damsel_pics)
new_damsel_pics <- list.files("./images/Mati_new_landmarks/")
extra_damsel_pics <- list.files("./images/extra_images/")

getSpNames <- function(list){
  list <- sub(" ", "_", list)
  unlist(lapply(strsplit(list, "_"), function(i) paste(i[1], i[2])))
}

all_damsel_pics <- setNames(all_damsel_pics, getSpNames(all_damsel_pics))
actual_damsel_pics <- setNames(actual_damsel_pics, getSpNames(actual_damsel_pics))
new_damsel_pics <- setNames(new_damsel_pics, getSpNames(new_damsel_pics))
extra_damsel_pics <- setNames(extra_damsel_pics, getSpNames(extra_damsel_pics))

listOfSpecies <- read.table("./data/listOfSpecies.txt")[,1]
damsel_tree <- read.tree("McCord2021_SupportingInfoFinal/S3_TreeFile_DamselOnlyTimeTree.phy")
tree_sps <- sub("_", " ", damsel_tree$tip.label)
all_pomacentridae <- apply(read.table("./data/all_pomacentridae_species_names.txt"), 1, paste, collapse = "_")

# check names on tree
check.names <- as.data.frame(t(sapply(damsel_tree$tip.label, function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
mispelled_tree <- mispelled_tree[c(1, 2, 4, 6, 7),]
damsel_tree$tip.label[damsel_tree$tip.label %in% rownames(mispelled_tree)] = unlist(mispelled_tree$suggested)

# check names 1st selection pictures
check.names <- as.data.frame(t(sapply(names(all_damsel_pics), function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
mispelled_tree <- mispelled_tree[grep("Stegastes|Plectroglyphidodon", rownames(mispelled_tree)),]
names(all_damsel_pics)[names(all_damsel_pics) %in% "Stegastes flaviatus"] <- "Stegastes flavilatus"
names(all_damsel_pics)[names(all_damsel_pics) %in% "Plectroglyphidodon leucozona"] <- "Plectroglyphidodon leucozonus"
names(all_damsel_pics)[names(all_damsel_pics) %in% "pomacentrus maafu"] <- "Pomacentrus maafu"

# check names extra images
check.names <- as.data.frame(t(sapply(names(extra_damsel_pics), function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
names(extra_damsel_pics)[names(extra_damsel_pics) %in% "Chromis sanctahelenae"] <- "Chromis sanctaehelenae"
names(extra_damsel_pics)[names(extra_damsel_pics) %in% "Pomacentrus caeruleopuncatus"] <- "Pomacentrus caeruleopunctatus"
names(extra_damsel_pics)[names(extra_damsel_pics) %in% "Pycnochromis lineata"] <- "Pycnochromis lineatus"
names(extra_damsel_pics)[names(extra_damsel_pics) %in% "pomacentrus maafu"] <- "Pomacentrus maafu"

# check names new landmarked images
check.names <- as.data.frame(t(sapply(names(new_damsel_pics), function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
names(new_damsel_pics)[names(new_damsel_pics) %in% "Stegastes flaviatus"] <- "Stegastes flavilatus"
names(new_damsel_pics)[names(new_damsel_pics) %in% "Plectroglyphidodon leucozona"] <- "Plectroglyphidodon leucozonus"
names(new_damsel_pics)[names(new_damsel_pics) %in% "pomacentrus maafu"] <- "Pomacentrus maafu"

# check names current landmarked images
check.names <- as.data.frame(t(sapply(names(actual_damsel_pics), function(i) check_spelling(i, all_pomacentridae))))
mispelled_tree <- check.names[check.names$score != 1,]
mispelled_tree
names(actual_damsel_pics)[names(actual_damsel_pics) %in% "Chrysiptera maurinae"] <- "Chrysiptera maurineae"
names(actual_damsel_pics)[names(actual_damsel_pics) %in% "Plectroglyphidodon leucozona"] <- "Plectroglyphidodon leucozonus"
names(actual_damsel_pics)[names(actual_damsel_pics) %in% "Pomacentrus caeruleopuncatus"] <- "Pomacentrus caeruleopunctatus"
names(actual_damsel_pics)[names(actual_damsel_pics) %in% "Pycnochromis lineata"] <- "Pycnochromis lineatus"
names(actual_damsel_pics)[names(actual_damsel_pics) %in% "pomacentrus maafu"] <- "Pomacentrus maafu"


current_damsel_pics <- c(actual_damsel_pics, new_damsel_pics)
current_pics_tab <- table(names(current_damsel_pics))
nrow(current_pics_tab)

# species with image not in tree
names(current_pics_tab)[!names(current_pics_tab) %in% sub("_", " ", damsel_tree$tip.label)]
# species in tree with no images
missing_sp_pic <- damsel_tree$tip.label[!sub("_", " ", damsel_tree$tip.label) %in% names(current_pics_tab)]


all_damsel_pics <- c(all_damsel_pics, extra_damsel_pics)
all_damsel_pics <- all_damsel_pics[names(all_damsel_pics) %in% sub("_", " ", damsel_tree$tip.label)]
all_pics_tab <- table(names(all_damsel_pics))
nrow(all_pics_tab)

missing_sp_pic %in% names(all_damsel_pics) # images for two missing species not found

current_damsel_pics <- current_damsel_pics[names(current_damsel_pics) %in% sub("_", " ", damsel_tree$tip.label)]
current_pics_tab <- table(names(current_damsel_pics))
nrow(current_pics_tab)
 
list_todo <- NULL
for (i in 1:nrow(current_pics_tab)) {
  sp <- names(current_pics_tab)[i]
  n <- current_pics_tab[i]
  
  if (n < 5 & sp %in% names(all_pics_tab)) {
    candidate_pics <- all_damsel_pics[names(all_damsel_pics) == sp]
    # remove already landmarked images
    candidate_pics <- candidate_pics[!sub("\\.[a-z]*$", "", candidate_pics) %in% sub("\\.[a-z]*$", "", current_damsel_pics)]
  
    # sample 5 - n images
    m <- 5 - n
    if (length(candidate_pics) >= m) {
      list_todo <- c(list_todo, sample(candidate_pics, m))
    } else {
      list_todo <- c(list_todo, candidate_pics)
    }
  }
}

length(list_todo)

write.csv(list_todo, file = "./images_2_landmark.txt")

dir.create("./images/images2landmark")
for (i in list_todo){
  if (i %in% list.files("./images/1st_selection_images/")) {
    #system(paste0("cp ", "./images/1st_selection_images/'", i, "' ./images/images2landmark"))
  }
  else if (i %in% list.files("./images/extra_images/"))  {
    #system(paste0("cp ", "./images/extra_images/'", i, "' ./images/images2landmark"))
  }
  else {
    print(i)
  }
}
