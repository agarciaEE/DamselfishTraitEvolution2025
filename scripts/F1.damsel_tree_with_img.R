library(ColorAR)
library(wesanderson)
library(ape)
library(raster)
library(phytools)
library(phyloch)

data(strat2012)

setwd("~/Unil/Research/5_Damselfish_evo/")

source("./Rscripts/custom_functions.R")

# Plot circular tree with fish examples of each genus
#####################################################

# load tree
tree <- read.tree("~/Unil/Research/5_Damselfish_evo//data/damsel_subset.tre")


## Load images filenames
img_filenames <- list.files("images/3rd_selection_images/", pattern = "\\.jpg$")
img_filenames <- stringr::str_replace_all(img_filenames, ".jpg", "")
length(img_filenames)

## load landmarks filenames
#lmk_filenames <- list.files("Mati_landmarks/", pattern = "\\.txt$")
lmk_filenames <- list.files("images/3rd_selection_landmarks/", pattern = "\\.txt$")
lmk_filenames <- stringr::str_replace_all(lmk_filenames, ".txt", "")
length(lmk_filenames)

# check match between landmarks and images
img_filenames[!img_filenames %in% lmk_filenames]
lmk_filenames[!lmk_filenames %in% img_filenames]

# subset images with landmarks
img_filenames <- img_filenames[img_filenames %in% lmk_filenames]

## Load images 
prepath <- 'images/3rd_selection_images'
extension <- '.jpg'
imageList = makeList(img_filenames, 'image', prepath, extension)

## load landmarks
prepath <- 'images/3rd_selection_landmarks'
extension <- '.txt'
landmarkList = makeList(lmk_filenames, 'landmark', prepath, extension)

# check number of landmarks are the correct (46)
names(landmarkList)[sapply(landmarkList, nrow) != 46]

# set internal landmarks to drop when building the mask
lndmks_2drop <- c(1:10, 45, 46)

# set store raster on memory to avoid issues
rasterOptions(todisk = F)

# name image and landmark list objects according to species names
names(imageList) <- gsub("-", "_", sub("_[0-9]*", "", names(imageList)))
names(landmarkList) <- gsub("-", "_", sub("_[0-9]*", "", names(landmarkList)))

# order landmarkList to match imageList
landmarkList <- landmarkList[names(imageList)]

selected_species <- tree$tip.label

check.names <- as.data.frame(t(sapply(names(imageList), function(i) check_spelling(i, selected_species))))
mispelled_img <- check.names[check.names$score != 1,]
mispelled_img
mispelled_img <- mispelled_img[c(6, 9, 11, 15),]
names(imageList)[names(imageList) %in% rownames(mispelled_img)] = unlist(mispelled_img$suggested)
names(landmarkList)[names(landmarkList) %in% rownames(mispelled_img)] = unlist(mispelled_img$suggested)

imageList <- imageList[selected_species]
landmarkList <- landmarkList[selected_species]

## removing bg of original images
imglist_nobg = sapply(selected_species, function(sp) {
  e <- extent(imageList[[sp]])
  lmk <- landmarkList[[sp]][-lndmks_2drop,]
  lmk[,2] <- e[4] - lmk[,2]
  mask <- SpatialPolygons( list(  Polygons(list(Polygon(lmk)), 1)))
  mask <- smoothr::smooth(mask, method = "ksmooth", smooth = 1)
  raster::mask(imageList[[sp]], mask)
})

par(mfrow = c(2,5))
sapply(imglist_nobg[sample(selected_species,10)], plotRGB)
## get color palette
palette <- unique(unlist(lapply(c("Zissou1", "Darjeeling1", "Royal1", "BottleRocket2", "Moonrise1", "AsteroidCity3", "GrandBudapest2",
                                  "Darjeeling2", "FantasticFox1"), function(i) as.character(wes_palette(i)))))

## Select species examples
genus <- unique(sapply(tree$tip.label, function(i) strsplit(i, "_")[[1]][1]))
species_bygenus.List <- sapply(genus, function(i) grep(i, tree$tip.label, value = T))

fish_examples <- tree$tip.label[seq(5, Ntip(tree)-5, 10)]

# get colors
palette <- sample(palette, length(genus))
# get images
img_examples <- imglist_nobg[fish_examples] # original images

tree2 <- tree
tree2$tip.label <- gsub("^([A-Za-z])[A-Za-z]*_([A-Za-z]+)", "\\1. \\2", tree2$tip.label)
# paint tree
for(i in 1:length(genus)){
  ii <- grep(genus[i], tree$tip.label)
  ca <- if(length(ii) > 1) {
    phytools::findMRCA(tree, tree$tip.label[ii])
  }  else { ii }
  tree <- phytools::paintSubTree(tree, ca, state = as.character(i),
                             anc.state= "0", stem=TRUE)
}
tol<-max(phytools::nodeHeights(tree))*1e-12
tree$maps <- lapply(tree$maps, function(x,tol) 
  if(length(x)>1) x[-which(x<tol)] else x,tol=tol)
cols <- setNames(c("grey", palette),0:length(genus))

### PLOT CIRCULAR
# plot trait on tips
plot(tree, color = cols, tips = FALSE, ftype="off", type = "fan", part = 0.99,
         fsize = 0.5, offset = 10, lwd = 3)
xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(tree)]
yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(tree)]
centr <- c(mean(xx), mean(yy))
xlim = get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim
ylim = get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim
plotSimmap(tree, color = cols, tips = FALSE, ftype="off", type = "fan", part = 0.88,
          lwd = 3, xlim = xlim * 1.5, ylim = ylim * 1.5)
xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(tree)]
yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(tree)]
# add ages lines
obj<-axis(1,pos=-2,at=seq(10,50,by=10),cex.axis=0.5,labels=FALSE)
text(obj,rep(-5,length(obj)), rev(obj),cex=0.6)
text(mean(obj),-8,"Million Yeas Ago",cex=0.8)
for(i in 1:(length(obj)-1)){
  a1<-atan(-2/obj[i])
  a2<-0.88*2*pi
  plotrix::draw.arc(0,0,radius=obj[i],a1,a2,lwd=1.5,
                    col=make.transparent("grey40",0.25))
}
par(bg = "transparent")
plotSimmap(tree2, color = cols, tips = FALSE, ftype="i", fsize = 0.4, offset = 5, type = "fan", part = 0.88,
           lwd = 3, xlim = xlim * 1.5, ylim = ylim * 1.5, add = T )
par(bg = "white")

## add fishes
xx <- setNames(xx[sapply(fish_examples, function(i) which(tree$tip.label == i))], tree$tip.label[sapply(fish_examples, function(i) which(tree$tip.label == i))])
yy <- setNames(yy[sapply(fish_examples, function(i) which(tree$tip.label == i))], tree$tip.label[sapply(fish_examples, function(i) which(tree$tip.label == i))])
angles <- sapply(1:length(img_examples), function(i) get_angle(centr[1], centr[2], xx[i], yy[i]))
plotImages((xx) * 1.35, (yy) * 1.35, img_examples, 
           width = 0.1, angle = angles)
text((xx) * 1.55, (yy) * 1.55, font = 3,
     gsub("^([A-Za-z])[A-Za-z]*_([A-Za-z]{3})[A-Za-z]*", "\\1.\\2", fish_examples))


### PLOT VERTICAL

# plot trait on tips
phytools::plotTree(tree2, xlim = c(0,75), mar = c(1.5,1,1,1),
     fsize = 0.5, offset = 0, lwd = 2, ftype = "off")
xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(tree)]
yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(tree)]
axisGeo(GTS = strat2012, tip.time = 0, unit ="epoch", ages = 
          TRUE, cex = 0.8, col = RColorBrewer::brewer.pal(9, "YlOrBr")[1:5], 
        texcol = "black", gridty = "longdash", offset = -2.5,
        gridcol = "grey60")
par(bg = "white")

## add fishes
xx <- setNames(xx[sapply(fish_examples, function(i) which(tree$tip.label == i))], tree$tip.label[sapply(fish_examples, function(i) which(tree$tip.label == i))])
yy <- setNames(yy[sapply(fish_examples, function(i) which(tree$tip.label == i))], tree$tip.label[sapply(fish_examples, function(i) which(tree$tip.label == i))])
yy0 <- yy
yy <- yy + sample(c(1, 2), length(yy), replace = T) * sample(c(-1, 1), length(yy), replace = T)

text((xx) * 1.2, (yy) , font = 3, adj = 0, cex = 1,
     gsub("^([A-Za-z])[A-Za-z]*_([A-Za-z]*)[A-Za-z]*", "\\1. \\2", fish_examples))

sapply(1:length(img_examples), function(i) 
  lines(x = c(xx[i],(xx)[i] * 1.05), y = c(yy0[i],yy[i]), lty = 2, col = "grey60"))


sapply(1:length(img_examples), function(i) plotImages(((xx) * 1.1)[i], yy[i], img_examples[i], 
                                                      width = 10, height = 12.5))

