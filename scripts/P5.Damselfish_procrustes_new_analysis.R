### Open relevant packages ###
library(FactoMineR)
library(factoextra)
library(ggplot2) # To draw
library(ggmap) # To draw maps (e.g. from Google Maps)
library(plyr) # To transform data set
library(geomorph) # To perform morphometric analyses
library(scatterplot3d) # To draw 3D-plots
library(abind) # To combine muti-dimesionnal arrays
library(ade4) # To perform multivariate analyses
library(jpeg)
library(readr)
library(gtools)
library(ggpubr)
library(gridExtra)
library(sp)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggmap)
require(cowplot)
library(grid)
library(PCAmixdata)
library(psych)
library(raster)
library(sf)
library(patternize)

setwd("~/Unil/Research/5_Damselfish_evo/")
source("./Rscripts/Damsel_functions.R")

###############################
### GEOMORPHOMETRY ANALYSES ###
###############################

## Load images filenames
img_filenames <- file.path("images/all_selected_images", list.files("images/all_selected_images/"))
img_names <- sub("\\.[a-z]*$", "", basename(img_filenames))

# remove duplicates
img_filenames <- img_filenames[!duplicated(img_names)]
img_names <- img_names[!duplicated(img_names)]
length(img_names)

## load landmarks filenames
lmk_filenames <- file.path("images/all_selected_landmarks", list.files("images/all_selected_landmarks/"))
lmk_names <- sub("\\.[a-z]*$", "", basename(lmk_filenames))

# remove duplicates
lmk_filenames <- lmk_filenames[!duplicated(lmk_names)]
lmk_names <- lmk_names[!duplicated(lmk_names)]
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

# set store raster on memory to avoid issues
rasterOptions(todisk = F)

# get species names and genera
sp_names <- sapply(img_names, function(i) strsplit(i, "_")[[1]][1])
sp_names <- sub("-", "_", sp_names)
sp_names <- sub(" ", "_", sp_names)
table(sp_names)

sp_genera <- sapply(sp_names, function(i) strsplit(i, "_")[[1]][1])
table(sp_genera)

selected_species <- as.character(read.table("./data/listOfSpecies.txt")[,1])

# Check species name matching between images and all taxons db
check.names <- as.data.frame(t(sapply(unique(sp_names), function(i) check_spelling(i, selected_species))))
mispelled_lmk <- check.names[check.names$score != 1,]
mispelled_lmk
sp_names[sp_names == "Stegastes_flaviatus"] = "Stegastes_flavilatus"
mispelled_lmk <- mispelled_lmk[mispelled_lmk$score > 0.6,]
for (i in 1:length(rownames(mispelled_lmk))) {
  sp_names[sp_names == rownames(mispelled_lmk)[i]] = unlist(mispelled_lmk$suggested)[i]
}

landmarkList <- landmarkList[sp_names %in% selected_species]
sp_names <- sp_names[sp_names %in% selected_species]

# check number of landmarks are the correct (46)
names(landmarkList)[sapply(landmarkList, nrow) != 46]

# drop fin-related landmarks to keep only body landmarks and avoid weird transformations
drop_lnmks <- c(9:10, 16:22, 25:29, 32:34, 37:39)
landmarkList_dropped <- lapply(landmarkList, function(x) x[-drop_lnmks, ])

## CREATING MATRIX FOR PROCRUSTES ANALYSES 
pic_names <- names(landmarkList_dropped)
Ldmks_array <- array(unlist(landmarkList_dropped), dim = c(nrow(landmarkList_dropped[[1]]), ncol(landmarkList_dropped[[1]]), length(landmarkList_dropped)))
dim(Ldmks_array)
dimnames(Ldmks_array)[[3]] <- pic_names

all_temp <- gpagen(Ldmks_array[,,])

all_m.temp <- matrix(NA,dim(Ldmks_array)[3],11) # 11 corresponds to the number of traits to be recorded...

## FUNCTIONS 
get_point <- function(pt, k) c(all_temp$coords[pt, 1, k], all_temp$coords[pt, 2, k])
all_dist.p <- function(pt1, pt2, k) dist(rbind(if (length(pt1) == 2) pt1 else get_point(pt1, k), 
                                               if (length(pt2) == 2) pt2 else get_point(pt2, k)))
all_mean.p <- function(pt1, pt2, k) colMeans(rbind(get_point(pt1, k), get_point(pt2, k)))

for(j in 1:dim(Ldmks_array)[3]){
  # 1 head length = 11 -> c(5,6)
  all_m.temp[j,1]<- all_dist.p(9,all_mean.p(5,6,j),j) # modified  
  
  # 2 snout length = 11 -> (1x,11y)
  all_m.temp[j,2]<- all_dist.p(get_point(9,j), c(get_point(1,j)[1], all_mean.p(5,6,j)[2]),j)  

  # 3 eye greatest diameter
  all_m.temp[j,3]<-max(all_dist.p(1,2,j),all_dist.p(3,4,j))
  
  # 4 standard length = 1 -> (24,30)
  all_m.temp[j,4]<-stl<-all_dist.p(9,all_mean.p(15,16,j),j)
  
  # 5 anal fin length = 31 -> 35
  all_m.temp[j,5]<-all_dist.p(17,18,j)
  
  # 6 body depth = mean (15,40 & 45,35)
  all_m.temp[j,6]<-mean(all_dist.p(13,20,j), all_dist.p(25,18,j))
  
  # 7 angle of snout profile = angle between standard length and 11 -> 15
  A<-rbind(get_point(9,j),all_mean.p(15,16,j))
  B<-rbind(get_point(9,j),get_point(13,j))
  all_m.temp[j,7]<-acos(sum(A*B)/(sqrt(sum(A*A))*sqrt(sum(B*B))))
  
  # 8 eye vertical position relative to head
  all_eye.mean<- apply(rbind(get_point(1,j),
                            get_point(2,j),
                            get_point(3,j),
                            get_point(4,j)),2,mean)
  # from horizontal line that crosses 24 and 30 to eye.mean
  all_eye.height<- all_eye.mean[2] - all_mean.p(15,16,j)[2]
  all_m.temp[j,8]<-all_m.temp[j,6]/all_eye.height

  # 9 peduncle area
  all_poly <- data.frame(all_temp$coords[c(14,15,16,17,14),1:2,j])
  polygon <- sf::st_polygon(list(as.matrix(all_poly)))  
  all_m.temp[j,9]<- sf::st_area(polygon)

  # 10 pelvic fin 36->40
  all_m.temp[j,10]<- all_dist.p(19,20,j)
  
  # 11 peduncle_depth 23->31
  all_m.temp[j,11]<- all_dist.p(14,17,j)
  
}
all_m <- all_m.temp

### To extract trait values from landmarks ###
all_morpho_ratio <- data.frame(body_ratio=all_m[,4]/all_m[,6], head_ratio=all_m[,4]/all_m[,1], snout_ratio=all_m[,1]/all_m[,2], eye_ratio=all_m[,1]/all_m[,3], peduncle_ratio=all_m[,1]/all_m[,9],
                               anal_fin_ratio=all_m[,4]/all_m[,5], pelvic_fin_ratio=all_m[,4]/all_m[,10], eye_height=all_m[,8], snout_angle=all_m[,7], peduncle_depth_factor=all_m[,6]/all_m[,11], relative_body_length=all_m[,4])

rownames(all_morpho_ratio) <- dimnames(Ldmks_array)[[3]]
all_morpho_ratio$species <- sp_names[rownames(all_morpho_ratio)]

saveRDS(all_temp, file = "./rdata/1166morho_gpagen.RDS")
write.csv(all_morpho_ratio, file = "./data/1166damsel_morpho.csv", row.names = TRUE)

eco_traits_subset <- read.csv("./data/eco_traits_subset_final.csv", row.names = 1)

all_morpho_ratio_sp <- as.data.frame(all_morpho_ratio %>%
  group_by(species) %>%
  dplyr::summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))))

rownames(all_morpho_ratio_sp) <- all_morpho_ratio_sp$species
all_morpho_ratio_sp <- all_morpho_ratio_sp[,-1]

all(rownames(all_morpho_ratio_sp) %in% rownames(eco_traits_subset))
eco_traits_subset <- cbind(eco_traits_subset[, !colnames(eco_traits_subset) %in% colnames(all_morpho_ratio_sp)],
                           all_morpho_ratio_sp[rownames(eco_traits_subset),])

write.csv(eco_traits_subset, file = "./data/eco_traits_subset_final.csv", row.names = T)

