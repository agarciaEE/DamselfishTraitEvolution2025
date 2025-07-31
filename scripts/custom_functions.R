## Custom Functions
####################
library(udpipe)
library(rvest)
library(raster)
library(sp)
library(foreach)
library(doParallel)
library(vegan)
library(phytools)
library(expm)
library(ape)

# copy of Patternize function of the same name but removing the crs assignment
makeList <- function (IDlist, type, prepath = NULL, extension = NULL, format = "imageJ", 
                      tpsFile = NULL, skipLandmark = NULL) 
{
  objectList <- list()
  if (!is.null(skipLandmark)) {
    skipLandmark <- -1 * skipLandmark
  }
  for (n in 1:length(IDlist)) {
    if (format == "imageJ") {
      print(paste("sample", n, IDlist[n], "added to list", 
                  sep = " "))
      if (type == "landmark") {
        if (is.null(prepath)) {
          landmarks <- read.table(paste(IDlist[n], extension, 
                                        sep = ""), header = FALSE, stringsAsFactors = FALSE, 
                                  colClasses = c("numeric", "numeric"))
        }
        else {
          landmarks <- read.table(paste(prepath, "/", 
                                        IDlist[n], extension, sep = ""), header = FALSE, 
                                  stringsAsFactors = FALSE, colClasses = c("numeric", 
                                                                           "numeric"))
        }
        landmarks <- as.matrix(landmarks)
        colnames(landmarks) <- NULL
        if (!is.null(skipLandmark)) {
          landmarks <- landmarks[skipLandmark, ]
        }
        objectList[[IDlist[n]]] <- landmarks
      }
    }
    if (type == "image") {
      if (is.null(prepath)) {
        suppressWarnings(image <- raster::stack(paste(IDlist[n], 
                                                      extension, sep = "")))
        #crs(image) <- sp::CRS("+init=EPSG:4326")
      }
      else {
        suppressWarnings(image <- raster::stack(paste(prepath, 
                                                      "/", IDlist[n], extension, sep = "")))
        #crs(image) <- sp::CRS("+init=EPSG:4326")
      }
      objectList[[IDlist[n]]] <- image
    }
  }
  if (all(c(type == "landmark", format == "tps"))) {
    objectListX <- readland.tps(tpsFile, specID = "imageID", 
                                warnmsg = FALSE)
    objectList <- lapply(1:dim(objectListX)[3], function(i) objectListX[, 
                                                                        , i])
    names(objectList) <- IDlist
  }
  return(objectList)
}

# add transparency to color
alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

# function to check color luminosity
check.luma <- function(col){
  #col = as.numeric(grDevices::col2rgb(col, alpha = FALSE))
  luma = 0.2126 * col[1] + 0.7152 * col[2] + 0.0722 * col[3]
  return(luma)
}

# function to build mean shape mask from a list of landmarks
build_mask <- function(landmaks, images, adjustCoords = TRUE, lndmks_2drop = NULL, smooth = 1){
  
  
  lanArray <- patternize::lanArray(landmarkList, adjustCoords = T, 
                                   imageList)
  
  invisible(utils::capture.output(transformed <- Morpho::procSym(lanArray)))
  refShape <- transformed$mshape
  if (!is.null(lndmks_2drop)){ refShape <- refShape[-lndmks_2drop, ] }
  
  sp.ref = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(refShape)), 1)))
  sp.ref <- smoothr::smooth(sp.ref, method = "ksmooth",  smooth = 1)
  sp.ref <- spatialEco::rotate.polygon(sp.ref, angle = 180)
  
  return(sp.ref)
}

# function to homogenize images and reduce NAs by interpolation
homogenize <- function(rasList, th = NULL, rasTemp = NULL, fact = 5, method = "bilinear", plot = TRUE, parallel = FALSE, cores = 2) {
  
  if (is.null(rasTemp)) {
    ext <- extent(rasList[[1]])
    ras_dim <- dim(rasList[[1]])
    rasTemp <- raster(nrows = ras_dim[1], ncol = ras_dim[2], ext = ext)
  }
  
  NAcounts <- sapply(rasList, function(i) sum(is.na(i[[1]][])))
  
  if (plot) {
    plotRGB(rasList[[names(NAcounts)[which.max(NAcounts)]]])
    text(as.numeric(dim(rasList[[1]]) / 2)[2], 0,
         paste0("Image with higher number of NAs (", NAcounts[which.max(NAcounts)], ")"))
  }
  
  if (is.null(th)) {
    NA_outliers <- boxplot.stats(NAcounts)$out
  } else {
    NA_outliers <- NAcounts[NAcounts > th]
  }
  
  message(length(NA_outliers), " images to homogenize...\n")
  
  # Register parallel backend if enabled
  if (parallel) {
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  }
  
  rasList_updated <- foreach(i = names(NA_outliers), .packages = c("raster", "sp"), .combine = "c") %dopar% {
    message("Interpolating ", i, " sample...")
    ras <- sp::disaggregate(rasList[[i]], fact = fact, method = method)
    ras <- raster::resample(ras, rasTemp, method = "ngb")
    list(ras) # Return as list for combining
  }
  
  if (parallel) {
    stopCluster(cl)
  }
  
  # Update the original rasList with the modified rasters
  for (i in seq_along(names(NA_outliers))) {
    rasList[[names(NA_outliers)[i]]] <- rasList_updated[[i]]
  }
  
  NAcounts_new <- sapply(rasList, function(i) sum(is.na(i[[1]][])))
  
  if (plot) {
    par(mfrow = c(1, 2))
    hist(NAcounts, main = "Distribution of NA counts")
    hist(NAcounts_new, main = "New distribution of NA counts")
    par(mfrow = c(1, 1))
  }
  
  return(rasList)
}

# check image landmakrs visually
check_landmarks <- function(imageList, landmarkList, imgTransList, txt_col = "yellow"){
  
  rep_lndmk <- list()
  stop = FALSE
  
  for (i in names(imageList)){
    
    img <- imageList[[i]]
    ext <- extent(img)
    lmk <- landmarkList[[i]]
    lmk[,2] <- ext[4] - lmk[,2]
    lmk_adj <- lmk
    
    par(mfrow = c(1,2))
    plotRGB(img)
    points(lmk, pch = 19, col = txt_col)
    text(lmk[,1], lmk[,2], 1:nrow(lmk), col = txt_col, adj = 2)
    plotRGB(imgTransList[[i]])
    par(mfrow = c(1,1))
    
    rep_lndmk_i <- NULL
    while (TRUE) {
      input <- readline(prompt="Introduce landmark number to repeat:")
      input <- as.integer(input)
      if (is.na(input)){
        break
      }
      else if (is.integer(input)){
        rep_lndmk_i <- c(rep_lndmk_i, input)
      }
      else if (input == "stop"){
        stop = TRUE
        break
      }
    }
    rep_lndmk[[i]] <- rep_lndmk_i
    if (stop) { break }
  }
  return(rep_lndmk)
}

# check species name misspelling
check_spelling <- function(spname, namesList, gs_sep = "_") {
  
  # Create a file with on each line one of the words we want to check
  #"[[:punct:]]"
  spname <- gsub(gs_sep," ", spname) # remove special characters
  spname <- strsplit(spname, " ", fixed=TRUE)[[1]] # split genus and species
  
  g <- strsplit(spname[1], "")[[1]]
  s <- strsplit(spname[2], "")[[1]]
  
  dict <- sapply(gsub("[[:punct:]]"," ", namesList) , function(i) strsplit(i, " ", fixed=TRUE)[[1]])
  
  match.score_genus <- apply(dict, 2, function(i) suppressWarnings(sum(g == strsplit(i[1], "")[[1]])/length(g)))
  match.length_genus <- scales::rescale(apply(dict, 2, function(i) 1 - abs(length(g) - nchar(i[1]))), c(0,1))
  score_genus <- match.score_genus + match.length_genus
  
  genus_match <- score_genus[which(score_genus == max(score_genus))]
  
  if (length(genus_match) == 1) {
    match.score_species <- suppressWarnings(sum(s == strsplit( dict[2,names(genus_match)], "")[[1]])/length(s))
    match.length_species <- 1
    
  } else {
    match.score_species <- apply(dict[,names(genus_match)], 2, function(i) suppressWarnings(sum(s == strsplit(i[2], "")[[1]])/length(s)))
    match.length_species <- scales::rescale(apply(dict[,names(genus_match)], 2, function(i) 1 - abs(length(s) - nchar(i[2]))), c(0.5,1))
  }
  
  score_species <- match.score_species * match.length_species
  
  res <- score_species[which(score_species == max(score_species))]
  
  if (length(res) > 1) {
    res <- res[which.max(match.score_species[names(res)])]
  }
  names(res) <- sub(" ", gs_sep, names(res))
  return(data.frame(score = res, suggested = names(res), 
                    row.names = paste(spname, collapse = gs_sep)))
}

# Function to extract nouns from a sentence
extract_nouns <- function(sentence) {
  # Tokenize the sentence
  tokens <- udpipe_annotate(ud_model, x = sentence)
  
  tokens <- as.data.frame(tokens)
  # Extract nouns
  nouns <- tokens$lemma[tokens$upos == "NOUN"]
  
  return(nouns)
}

# Function to extract nouns along with preceding adjectives from a sentence
extract_nouns_with_adj <- function(sentence) {
  # Tokenize the sentence
  tokens <- udpipe_annotate(ud_model, x = sentence)
  
  tokens <- as.data.frame(tokens)
  
  # Extract nouns along with preceding adjectives
  nouns_with_adj <- c()
  i <- 1
  while (i <= length(tokens$lemma)) {
    # Check if the current token is a noun
    if (tokens$upos[i] == "NOUN") {
      # Look for preceding adjectives
      adjectives <- c()
      j <- i - 1
      while (j >= 1 && tokens$upos[j] %in% c("ADJ", "CCONJ")) {
        if (tokens$upos[j] == "ADJ"){
          adjectives <- c(tokens$lemma[j], adjectives)
        }
        j <- j - 1
      }
      # Combine adjectives and noun with hyphen
      if (length(adjectives) > 0) {
        nouns_with_adj <- c(nouns_with_adj, paste0(rev(adjectives), "-", tokens$lemma[i]))
      } else {
        nouns_with_adj <- c(nouns_with_adj, tokens$lemma[i])
      }
    }
    i <- i + 1
  }
  
  return(nouns_with_adj)
}

# function to extract habitat preference and source link from Reef Life Survey website
extract_rls_info <- function(species_name) {
  # Construct the URL for the species on the Reef Life Survey website
  url <- paste0("https://reeflifesurvey.com/species/", gsub("_| ", "-", species_name), "/")
  #url <- paste0("https://fishbase.mnhn.fr/summary/", tolower(gsub("_| ", "-", species_name)), "/")
  
  # Read the HTML content of the webpage
  webpage <- tryCatch(read_html(url), error = function(e) NA)
  
  if (!is.na(webpage)){
    # Extract habitat preference
    txt <- webpage %>%
      html_text() %>%
      trimws()
    
    txt <- strsplit(txt, "\\{")[[1]][166]
    txt <- strsplit(txt, "Distribution")[[1]][2]
    
    dis <- strsplit(txt, "Description")[[1]][1]
    txt <- strsplit(txt, "Description")[[1]][2]
    
    des <- strsplit(txt, "Information")[[1]][1]
    txt <- strsplit(txt, "Information")[[1]][2]
    
    txt <- strsplit(txt, "Max Size: ")[[1]][2]
    
    MaxSize <- strsplit(txt, "Sea Temperature Range: ")[[1]][1]
    txt <- strsplit(txt, "Sea Temperature Range: ")[[1]][2]
    
    SeaTemp <- strsplit(txt, "Depth: ")[[1]][1]
    txt <- strsplit(txt, "Depth: ")[[1]][2]
    
    Depth <- strsplit(txt, "Habitat Generalization Index: ")[[1]][1]
    txt <- strsplit(txt, "Habitat Generalization Index: ")[[1]][2]
    
    HGI <- as.numeric(substr(txt, 1, 5))
    txt <- strsplit(txt, "Status: ")[[1]][2]
    
    IUCN <- strsplit(txt, "Occurrence: ")[[1]][1]
    txt <- strsplit(txt, "Occurrence: ")[[1]][2]
    
    occ <- strsplit(strsplit(txt, "Abundance: ")[[1]][1], " ")[[1]][1]
    txt <- strsplit(txt, "Abundance: ")[[1]][2]
    
    abu <- strsplit(txt, " ")[[1]][1]
    
    # Extract source link
    source_link <- url
    
    res <- data.frame(species = species_name,
                      distribution = dis,
                      maxsize = MaxSize,
                      seatemp = SeaTemp,
                      depth = Depth,
                      HGI = HGI,
                      IUCN = IUCN,
                      occurrence = occ,
                      abundance = abu,
                      source = url)
    
    return(res)
  } 
}

# function to rescale a tree
rescaleTree<-function(tree, scale){
  tree$edge.length<- tree$edge.length / max(nodeHeights(tree)) * scale
  tree
}

# Adjusted densityTree function
densityTree<-function(trees,colors="blue",alpha=NULL,method="plotTree",
                      fix.depth=FALSE,use.edge.length=TRUE,compute.consensus=TRUE,
                      use.gradient=FALSE,show.axis=TRUE,...){
  N<-length(trees)
  if(any(sapply(trees,function(x) is.null(x$edge.length)))) 
    use.edge.length<-FALSE
  if(!use.edge.length) trees<-lapply(trees,compute.brlen)
  h<-sapply(trees,function(x) max(nodeHeights(x)))
  if(fix.depth){
    if(method=="plotTree"){
      trees<-lapply(trees,rescaleTree,mean(h))
      class(trees)<-"multiPhylo"
    } else if(method=="plotSimmap"){ 
      trees<-rescaleSimmap(trees,depth=mean(h))
      print(class(trees))
      print(class(trees[[1]]))
    }
    h<-sapply(trees,function(x) max(nodeHeights(x)))
  }
  tips<-setNames(1:Ntip(trees[[1]]), 
                 if(compute.consensus) untangle(consensus(trees),
                                                "read.tree")$tip.label 
                 else trees[[1]]$tip.label)
  if(is.null(alpha)) alpha<-max(c(1/N,0.01))
  args<-list(...)
  args$direction<-"leftwards"
  args$tips<-tips
  args$add<-FALSE
  if(is.null(args$nodes)) args$nodes<-"inner"
  if(is.null(args$mar)) args$mar<-if(show.axis) c(4.1,1.1,1.1,1.1) else rep(1.1,4)
  if(is.null(args$ftype)) args$ftype<-"i"
  if(!use.gradient){
    plotTree(trees[[which(h==max(h))[1]]],direction="leftwards",mar=args$mar,
             plot=FALSE)
    args$xlim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim[2:1]
    if(method=="plotTree"){
      args$color<-make.transparent(colors[1],alpha)
      for(i in 1:length(trees)){
        args$tree<-trees[[i]]
        do.call(plotTree,args)
        if(i==1){ 
          if(show.axis) axis(1)
          args$ftype<-"off"
          args$add<-TRUE
        }
      }
    } else if(method=="plotSimmap"){
      states<-sort(unique(as.vector(mapped.states(trees))))
      print(states)
      print(colors)
      if(length(colors)!=length(states)){
        colors<-setNames(c("grey",palette()[2:length(states)]),
                         states)
      }
      colors<-sapply(colors,make.transparent,alpha=alpha)
      args$colors<-colors
      for(i in 1:length(trees)){
        args$tree<-trees[[i]]
        do.call(plotSimmap,args)
        if(i==1){ 
          if(show.axis) axis(1)
          args$ftype<-"off"
          args$add<-TRUE
        }
      }
    }
  } else if(use.gradient){
    nulo<-capture.output(rf<-multiRF(trees),type="message")
    mds<-cmdscale(rf,k=1)[,1]
    trees<-trees[order(mds)]
    args$ylim<-c(0,Ntip(trees[[1]])+1)
    plotTree(trees[[which(h==max(h))[1]]],direction="leftwards",mar=args$mar,
             ylim=args$ylim,plot=FALSE)
    args$xlim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim[2:1]
    colors<-sapply(rainbow(n=length(trees)),make.transparent,alpha=alpha)
    ftype<-args$ftype
    for(i in 1:length(trees)){
      y.shift<-(i-median(1:length(trees)))/length(trees)/2
      args$tree<-trees[[i]]
      args$tips<-tips+y.shift
      args$color<-colors[i]
      args$ftype<-if(i==floor(median(1:length(trees)))) ftype else "off"
      do.call(plotTree,args)
      if(i==1){ 
        if(show.axis) axis(1)
        args$ftype<-"off"
        args$add<-TRUE
      }
    }
  }
}

# function to make a 3-state continuous legend
make.3state.legend <- function(cols, names = NULL, fill.center = NULL, add = TRUE,
                               size = 1.5, border.col = c("transparent", "color"), 
                               xlims = c(-1, 1), ylims = c(-1, 1), pc = 1, b = 100,
                               cex = 1, offset = 0.5, ...){
  
  if (length(cols) != 3) {
    stop("Number of colors is not 3.")
  } else {
    col1 <- cols[1]
    col2 <- cols[2]
    col3 <- cols[3]
  }
  
  if (!add){
    plot(NULL, xaxt="n", yaxt="n", type="n", xlab="", ylab="",  bty="n",
         xlim = xlims + diff(xlims)*0.1 * c(-1,1),
         ylim = ylims + diff(ylims)*0.1 * c(-1,1))
  }
  
  border.col <- border.col[1]
  plot.center = FALSE
  ## make a triangle of cell coordinates
  cells.x<-cells.y<-vector()
  k<-1
  x0 <- 100 # start x coordinates
  asp<-abs(diff(ylims)/diff(xlims))*par()$pin[2]/par()$pin[1]
  for(i in 1:(b/2)){
    for(j in i:(b-i+1)){
      cells.x[k]<-x0+(j-1)/100
      cells.y[k]<-b+(i-1)*2*asp/100
      k<-k+1
    }
  }
  pch <- rep(c(24, 25), length(cells.y)/2)
  size.correction <- (par("cin")[2]/par("pin")[2])*(ylims[2] - ylims[1]) * 0.25 * size * pc
  
  # adjust points for pch 25
  cells.y[1:length(cells.y)%%2==0] = cells.y[1:length(cells.y)%%2==0] + abs(size.correction)
  
  # remove extra diagonal
  res <- cumsum(c(b, b - (1:(b))[1:(b)%%2 == 0]))
  #res <- cumsum(rev((1:(b/2))[1:(b/2)%%2==0]))
  cells.x <- cells.x[-res]
  cells.y <- cells.y[-res]
  pch <- pch[-res]
  
  # get weigths for each color:
  # color 1 (left)
  hm <- cells.y[which.min(abs(cells.y - (min(cells.y) + diff(range(cells.y))*3/5)))]
  xm1 <- max(cells.x[cells.y == hm])
  w1 <- cells.x
  w1[cells.y > hm] = 0
  w1[cells.y <= hm] = 1 - scales::rescale(cells.y[cells.y <= hm], to = c(0,1))
  w1[cells.x > xm1] = 0
  w1[cells.x <= xm1] = w1[cells.x <= xm1] * (1 - scales::rescale(cells.x[cells.x <= xm1], to = c(0,1)))
  # color 2 (rigth)
  xm2 <- min(cells.x[cells.y == hm])
  w2 <- cells.x
  w2[cells.y > hm] = 0
  w2[cells.y <= hm] = 1 - scales::rescale(cells.y[cells.y <= hm], to = c(0,1))
  w2[cells.x < xm2] = 0
  w2[cells.x >= xm2] = w2[cells.x >= xm2] * (scales::rescale(cells.x[cells.x >= xm2], to = c(0,1)))
  # color 3 (top)
  hm <- cells.y[which.min(abs(cells.y - (min(cells.y) + diff(range(cells.y))*2/5)))]
  w3 <- cells.x
  w3[cells.y < hm] = 0
  w3[cells.y >= hm] = scales::rescale(cells.y[cells.y >= hm], to = c(0,1))
  
  if (!is.null(fill.center)){
    # get center
    cells.x <- scales::rescale(cells.x, to = c(0, 100))
    cells.y <- scales::rescale(cells.y, to = c(0, 100))
    center_x <- mean(cells.x)
    center_y <- mean(cells.y)
    wc <- log(sqrt((cells.x - center_x)^2 + (cells.y - center_y)^2))
    wc[wc<0] = 0
    wc <- scales::rescale(wc, to = c(0,1))
    
    w1 <- scales::rescale(w1 * wc, to = c(0,1))
    w2 <- scales::rescale(w2 * wc, to = c(0,1))
    w3 <- scales::rescale(w3 * wc, to = c(0,1))
    
    colsc <- sapply(wc, function(i) alpha(fill.center, 1-i))
    
    plot.center = TRUE
    
  }
  
  cols1 <- sapply(w1, function(i) alpha(col1, i))
  cols2 <- sapply(w2, function(i) alpha(col2, i))
  cols3 <- sapply(w3, function(i) alpha(col3, i))
  
  if (border.col == "color") {
    bcols1 <- cols1
    bcols2 <- cols2
    bcols3 <- cols3
  } else {
    bcols1 <- bcols2 <- bcols3 <- border.col
  }
  
  cells.x <- scales::rescale(cells.x, to = xlims)
  cells.y <- scales::rescale(cells.y, to = ylims)
  
  points(cells.x, cells.y, pch = pch, cex = size, bg = cols1, col = bcols1)
  points(cells.x, cells.y, pch = pch, cex = size, bg = cols2, col = bcols2)
  points(cells.x, cells.y, pch = pch, cex = size, bg = cols3, col = bcols3)
  if (plot.center) points(cells.x, cells.y, pch = pch, cex = size, bg = colsc, col = colsc)
  
  if (!is.null(names)) {
    
    x.offset <- diff(xlims)*0.1 * offset
    y.offset <- diff(ylims)*0.1 * offset 
    
    if (ylims[1] > ylims[2] ){
      y <- c(max(cells.y), max(cells.y), min(cells.y)) - strheight(names[1]) * cex * c(-1, -1, 1)
    } else {
      y <- c(min(cells.y), min(cells.y), max(cells.y)) + strheight(names[1]) * cex * c(-1, -1, 1)
    }
    if (xlims[1] > xlims[2] ){
      x <- c(max(cells.x), min(cells.x), max(cells.x) + diff(xlims)/2) - strwidth(names[1]) * cex * c(-1, 1, 0)
    } else {
      x <- c(min(cells.x), max(cells.x), min(cells.x) + diff(xlims)/2) + strwidth(names[1]) * cex * c(-1, 1, 0)
    }
    
    text(x + x.offset * c(-1, 1, 0), 
         y + y.offset * c(-1, -1, 1), names, cex = cex)
  }
}

# color map 3-state trait
densityMap.3states <- function(trees, states, colors = c("red", "blue", "green"), 
                               xlab = "Million Years Ago", x.labels = NULL,
                               cex.lab = 3, lwd = 2, cex.axis = 2, legend.pos = "topleft",
                               fill.center = "grey90", legend.cex = 1){
  
  fsize <- 36*par()$pin[2]/par()$pin[1]/Ntip(trees)[1]
  densityTree(trees, method="plotSimmap", fix.depth = TRUE, ftype = "off",
              fsize = fsize, lwd=8, alpha = NULL, nodes="intermediate",
              colors=colors, compute.consensus=FALSE, show.axis = FALSE)
  
  # create a white rectangle to hide anythign going beyond the tips
  rect(-0.01, par()$usr[3], par()$usr[2], par()$usr[4], col = "white", border = "white")
  
  # x axis
  if ( "multiPhylo" %in% class(trees)){  h<-max(nodeHeights(trees[[1]])) } 
  else { h<-max(nodeHeights(trees)) }
  if (is.null(x.labels)) x.labels <- seq(0, h, by=5)
  axis(1, pos=-0.025*h, at=x.labels, labels=x.labels, tick = TRUE,  cex.axis=cex.axis, lwd=lwd, lend=2)
  text(mean(par()$usr[1:2]), -0.175*h, xlab, cex = cex.lab)
  
  # legend
  if (grepl("top", legend.pos, fixed = TRUE)) {
    ylims <- sort(quantile(par()$usr[3:4], c(0.95, 0.85)))
    
  } else if (grepl("bottom", legend.pos, fixed = TRUE)) {
    ylims <- sort(quantile(par()$usr[3:4], c(0.05, 0.15)))
  }
  r <- diff(ylims)/diff(par()$usr[3:4]) * 1.5
  a <- diff(par()$usr[1:2]) * r
  if (grepl("left", legend.pos, fixed = TRUE)) {
    xlims <- quantile(par()$usr[1:2], c(0.95, 0.8))
  } else if (grepl("right", legend.pos, fixed = TRUE)) {
    xlims <- quantile(par()$usr[1:2], c(0.05, 0.2))
  }
  xlims <- c(xlims[1], xlims[1] + a)
  
  make.3state.legend(colors, names = states, fill.center = fill.center, 
                     xlims = xlims, ylims = ylims, size = legend.cex, pc = 1, b = 100,
                     border.col = "transparent", offset = 0, cex = cex.lab, add = T)
  
}

get_state_transitions <- function(tree, states) {
  transitions <- rep(FALSE, length(states))
  names(transitions) <- names(states)
  
  for (edge in 1:nrow(tree$edge)) {
    parent <- tree$edge[edge, 1]
    child  <- tree$edge[edge, 2]
    
    parent_state <- states[as.character(parent)]
    child_state  <- states[as.character(child)]
    
    # Mark TRUE only if the state changes from parent to child
    if (!is.na(child_state) && child_state != parent_state) {
      transitions[as.character(child)] <- TRUE
    }
  }
  
  # Also mark the root as a transition if it has a state
  root <- getRoot(tree)
  if (!is.na(states[as.character(root)])) {
    transitions[as.character(root)] <- TRUE
  }
  
  return(transitions)
}

# function to write a table with regions into BioGeoBEARS geographic data input format
write_BioGeoBEARS_geodata <- function(tab, filename = "damsel"){
  
  ntaxa <- nrow(tab)
  nareas <- ncol(tab)
  areas <- colnames(tab)
  codes <- LETTERS[1:length(areas)]
  
  # create description file
  des <- data.frame(code = codes, area = areas)
  des$ntaxa <- apply(tab, 2, sum, na.rm = T)
  
  # write description file
  write.table(des, file = paste0(filename, "_geodata_descriptionFile.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  data <-  data.frame(species = sub(" ", "_", rownames(tab)),
                      areas = apply(tab, 1, paste, collapse = ""), 
                      sep = "\n")
  # write BioGeoBEARS geographical input data
  cat(ntaxa, nareas, paste0("(", paste(codes, collapse = " "), ")"), file = paste0(filename, "_geodata_input.txt"))
  cat("\n", file = paste0(filename, "_geodata_input.txt"), append = TRUE)
  apply(data, 1, cat, file = paste0(filename, "_geodata_input.txt"), append = TRUE)
  invisible(list(description = des, data = data))
  
}

# adjusted wrapper function for BioGeoBears results
my_plot_BioGeoBEARS_results <- function(results_object, analysis_titletxt = NULL, addl_params = list(), 
                                        plotwhat = "text", label.offset = NULL, tipcex = 0.8, statecex = 0.7, 
                                        splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE, plotlegend = FALSE, 
                                        legend_ncol = NULL, legend_cex = 1, cornercoords_loc = "auto", 
                                        tr = NULL, tipranges = NULL, if_ties = "takefirst", pie_tip_statecex = 0.7, 
                                        juststats = FALSE, xlab = "Millions of years ago", root.edge = TRUE, 
                                        colors_list_for_states = NULL, skiptree = FALSE, show.tip.label = TRUE, 
                                        tipcol = "black", dej_params_row = NULL, plot_max_age = NULL, 
                                        skiplabels = FALSE, plot_stratum_lines = TRUE, include_null_range = NULL, 
                                        plot_null_range = FALSE, simplify_piecharts = FALSE, tipboxes_TF = TRUE, 
                                        tiplabel_adj = c(0.5), no.margin = FALSE, xlims = NULL, ylims = NULL, 
                                        cex.lab = 0.8, cex.axis = 0.8, axis.lwd = 1, title_offset = 1) 
{
  junk = "\n\t# manual_ranges_txt=NULL, \n\t# @manual_ranges_txt If you dont want to use the default text for each range, produced\n\t# by areas_list_to_states_list_new(), specify the list here.\n\n\t\n\tscriptdir = \"/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/a_scripts/\"\n\tplot_BioGeoBEARS_results(results_object, analysis_titletxt=NULL, addl_params=list(), plotwhat=\"text\", label.offset=NULL, tipcex=0.8, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=NULL, tipranges=NULL)\n\t\n\t# Defaults\n\taddl_params=list(\"j\"); plotwhat=\"text\"; label.offset=0.45; tipcex=0.7; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges; juststats = FALSE; plotlegend=FALSE; \txlab=\"Millions of years ago\"; if_ties=\"takefirst\"\n\t\n\t\n\t# Setup\nresults_object = resDEC\nanalysis_titletxt =\"BioGeoBEARS DEC on Mariana M1v4_unconstrained\"\naddl_params=list(\"j\"); plotwhat=\"text\"; label.offset=0.45; tipcex=0.7; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges\njuststats=FALSE; plotlegend=FALSE; \txlab=\"Millions of years ago\"; if_ties=\"takefirst\"\nshow.tip.label=TRUE\ntipcol=\"black\"; dej_params_row=NULL; plot_max_age=NULL; skiplabels=FALSE; \ncolors_list_for_states=NULL\nskiptree=FALSE\ninclude_null_range=NULL\nplot_stratum_lines=TRUE\n\tplot_null_range = FALSE\n\t"
  if (is.null(include_null_range) == TRUE) {
    include_null_range = results_object$inputs$include_null_range
  }
  results_object$inputs$include_null_range = include_null_range
  tmp_fg = par("fg")
  par(fg = "black")
  BioGeoBEARS_run_object = results_object$inputs
  if (is.null(tr)) {
    tr = check_trfn(trfn = BioGeoBEARS_run_object$trfn)
  }
  tr_pruningwise = reorder(tr, "pruningwise")
  tips = 1:length(tr_pruningwise$tip.label)
  nodes = (length(tr_pruningwise$tip.label) + 1):(length(tr_pruningwise$tip.label) + 
                                                    tr_pruningwise$Nnode)
  if (is.null(tipranges)) {
    if (BioGeoBEARS_run_object$use_detection_model == FALSE) {
      tipranges = getranges_from_LagrangePHYLIP(lgdata_fn = np(BioGeoBEARS_run_object$geogfn))
    }
    if (BioGeoBEARS_run_object$use_detection_model == TRUE) {
      if (BioGeoBEARS_run_object$use_detection_model == 
          TRUE) {
        tipranges = tipranges_from_detects_fn(detects_fn = BioGeoBEARS_run_object$detects_fn)
      }
    }
  }
  areas = getareas_from_tipranges_object(tipranges)
  numareas = length(areas)
  if (!is.na(results_object$inputs$max_range_size)) {
    max_range_size = results_object$inputs$max_range_size
  }
  else {
    max_range_size = length(areas)
  }
  max_range_size
  if (is.null(results_object$inputs$states_list)) {
    numstates = numstates_from_numareas(numareas = length(areas), 
                                        maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
    states_list_areaLetters = areas_list_to_states_list_new(areas, 
                                                            maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
    states_list_0based_index = rcpp_areas_list_to_states_list(areas, 
                                                              maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
  }
  else {
    states_list_0based_index = results_object$inputs$states_list
  }
  param_ests = extract_params_from_BioGeoBEARS_results_object(results_object, 
                                                              returnwhat = "table", addl_params = addl_params, paramsstr_digits = 4)
  if (juststats == TRUE) {
    return(param_ests)
  }
  else {
    paramstr = extract_params_from_BioGeoBEARS_results_object(results_object, 
                                                              returnwhat = "string", addl_params = addl_params, 
                                                              paramsstr_digits = 4)
  }
  param_names = extract_params_from_BioGeoBEARS_results_object(results_object, 
                                                               returnwhat = "param_names", addl_params = addl_params, 
                                                               paramsstr_digits = 4)
  if (is.null(analysis_titletxt)) {
    tmptxt = results_object$inputs$description
    if (any(is.null(tmptxt), tmptxt == "", tmptxt == "defaults", 
            tmptxt == "default")) {
      analysis_titletxt = ""
    }
    else {
      analysis_titletxt = results_object$inputs$description
    }
  }
  if (is.null(dej_params_row)) {
    analysis_titletxt = paste(analysis_titletxt, "\n", "ancstates: global optim, ", 
                              max_range_size, " areas max. ", paramstr, sep = "")
  }
  else {
    brate_col_TF = names(dej_params_row) == "brate"
    brate_col = (1:length(dej_params_row))[brate_col_TF]
    biogeog_params = dej_params_row[1:(brate_col - 1)]
    biogeog_param_names = names(dej_params_row)[1:(brate_col - 
                                                     1)]
    equals_col = "="
    tmpcols = cbind(biogeog_param_names, equals_col, unlist(biogeog_params))
    txtrows = apply(X = tmpcols, MARGIN = 1, FUN = paste, 
                    sep = "", collapse = "")
    txtrows
    biogeog_params_txt = paste(txtrows, sep = "", collapse = "; ")
    titletxt2 = bquote(paste(.(max_range_size), " areas max., ", 
                             .(biogeog_params_txt), "; ", lambda, "=", .(dej_params_row$brate), 
                             "; ", mu, "=", .(dej_params_row$drate), "; ", alpha, 
                             "=", .(dej_params_row$rangesize_b_exponent), "; ", 
                             omega, "=", .(dej_params_row$rangesize_d_exponent), 
                             "", sep = ""))
  }
  leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr_pruningwise)
  marprobs = results_object$ML_marginal_prob_each_state_at_branch_bottom_below_node
  left_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 
                                                              2], ]
  right_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 
                                                               1], ]
  right_ML_marginals_by_node
  if (is.null(dim(left_ML_marginals_by_node))) {
    left_ML_marginals_by_node = matrix(data = left_ML_marginals_by_node, 
                                       nrow = 1)
  }
  if (is.null(dim(right_ML_marginals_by_node))) {
    right_ML_marginals_by_node = matrix(data = right_ML_marginals_by_node, 
                                        nrow = 1)
  }
  relprobs_matrix = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
  if (length(nodes) > 1) {
    relprobs_matrix_for_internal_states = relprobs_matrix[nodes, 
    ]
  }
  else {
    relprobs_matrix_for_internal_states = relprobs_matrix[nodes, 
    ]
    relprobs_matrix_for_internal_states = matrix(data = relprobs_matrix_for_internal_states, 
                                                 nrow = 1, ncol = ncol(relprobs_matrix))
  }
  relprobs_matrix
  if (is.null(states_list_0based_index)) {
    statenames = areas_list_to_states_list_new(areas, maxareas = max_range_size, 
                                               include_null_range = results_object$inputs$include_null_range, 
                                               split_ABC = FALSE)
    ranges_list = as.list(statenames)
  }
  else {
    ranges_list = states_list_0based_to_ranges_txt_list(state_indices_0based = states_list_0based_index, 
                                                        areanames = areas)
    statenames = unlist(ranges_list)
  }
  MLprobs = get_ML_probs(relprobs_matrix)
  MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, 
                                         returnwhat = "states", if_ties = if_ties)
  if (is.null(colors_list_for_states)) {
    colors_matrix = get_colors_for_numareas(length(areas))
    colors_list_for_states = mix_colors_for_states(colors_matrix, 
                                                   states_list_0based_index, plot_null_range = results_object$inputs$include_null_range)
  }
  if (is.null(ranges_list)) {
    possible_ranges_list_txt = areas_list_to_states_list_new(areas, 
                                                             maxareas = max_range_size, split_ABC = FALSE, include_null_range = results_object$inputs$include_null_range)
  }
  else {
    possible_ranges_list_txt = ranges_list
  }
  cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, 
                                    colors_list_for_states, MLstates)
  if (plotlegend == TRUE) {
    colors_legend(possible_ranges_list_txt, colors_list_for_states, 
                  legend_ncol = legend_ncol, legend_cex = legend_cex)
  }
  if (root.edge == FALSE) {
    tr$root.edge = 0
  }
  if (root.edge == TRUE) {
    if (is.null(tr$root.edge) == TRUE) {
      tr$root.edge = 0
    }
  }
  if (is.null(label.offset)) {
    label.offset = 0.007 * (get_max_height_tree(tr) + tr$root.edge)
  }
  if (show.tip.label == TRUE) {
    if (is.null(plot_max_age)) {
      max_x = 1.25 * (get_max_height_tree(tr) + tr$root.edge)
      min_x = 0
    }
    else {
      nontree_part_of_x = plot_max_age - (get_max_height_tree(tr) + 
                                            tr$root.edge)
      max_x = 1.25 * (get_max_height_tree(tr) + tr$root.edge)
      min_x = -1 * nontree_part_of_x
    }
  }
  else {
    if (is.null(plot_max_age)) {
      max_x = 1.05 * (get_max_height_tree(tr) + tr$root.edge)
      min_x = 0
    }
    else {
      nontree_part_of_x = plot_max_age - (get_max_height_tree(tr) + 
                                            tr$root.edge)
      max_x = 1.05 * (get_max_height_tree(tr) + tr$root.edge)
      min_x = -1 * nontree_part_of_x
    }
  }
  max_tree_x = 1 * (get_max_height_tree(tr) + tr$root.edge)
  if (is.null(xlims)) {
    xlims = c(min_x, max_x)
  }
  else {
    xlims = xlims
  }
  nodecoords = node_coords(tr, tmplocation = cornercoords_loc, 
                           root.edge = root.edge)
  max_tree_x = max(nodecoords$x)
  if (is.null(plot_max_age)) {
    xticks_desired_lims = c(0, max_tree_x)
  }
  else {
    xticks_desired_lims = c(0, plot_max_age)
  }
  xticks_desired = pretty(xticks_desired_lims)
  xaxis_ticks_locs = max_tree_x - xticks_desired
  if (skiptree != TRUE) {
    plot(tr_pruningwise, x.lim = xlims, y.lim = ylims, show.tip.label = FALSE, 
         label.offset = label.offset, cex = tipcex, no.margin = no.margin, 
         edge.width = axis.lwd, root.edge = root.edge)
    if (show.tip.label == TRUE) {
      tiplabels_to_plot = sapply(X = tr_pruningwise$tip.label, 
                                 FUN = substr, start = 1, stop = 30)
      if (skiplabels == FALSE) {
        tiplabels(text = tiplabels_to_plot, tip = tips, 
                  cex = tipcex, adj = 0, bg = "white", frame = "n", 
                  pos = 4, offset = label.offset, col = tipcol)
      }
    }
    axis(side = 1, at = xaxis_ticks_locs, labels = xticks_desired, 
         cex.axis = cex.axis, lwd = axis.lwd, lwd.ticks = axis.lwd)
    mtext(text = xlab, side = 1, line = 20, cex = cex.lab)
  }
  if (plotwhat == "text") {
    par(fg = tmp_fg)
    if (skiplabels == FALSE) {
      nodelabels(text = MLstates[nodes], node = nodes, 
                 bg = cols_byNode[nodes], cex = statecex)
      tiplabels(text = MLstates[tips], tip = tips, bg = cols_byNode[tips], 
                cex = statecex, adj = tiplabel_adj)
    }
    par(fg = "black")
  }
  if (plotwhat == "pie") {
    par(fg = tmp_fg)
    if (skiplabels == FALSE) {
      if (simplify_piecharts == TRUE) {
        colnums_to_keep_in_probs = NULL
        probs = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
        probs2 = probs
        maxprob = rep(0, nrow(probs))
        other = rep(0, nrow(probs))
        num_to_keep = 1
        cat("\nSince simplify_piecharts==TRUE, reducing prob pie charts to (most probable, other)...\n")
        for (i in 1:nrow(probs)) {
          cat(i, " ", sep = "")
          tmprow = probs[i, ]
          positions_highest_prob_to_lowest = rev(order(tmprow))
          positions_to_keep = positions_highest_prob_to_lowest[1:num_to_keep]
          colnums_to_keep_in_probs = c(colnums_to_keep_in_probs, 
                                       positions_to_keep)
          keepTF = rep(FALSE, length(tmprow))
          keepTF[positions_to_keep] = TRUE
          otherTF = keepTF == FALSE
          other[i] = sum(tmprow[otherTF])
          tmprow[otherTF] = 0
          probs2[i, ] = tmprow
        }
        cat("\n")
        colnums_to_keep_in_probs_in_order = sort(unique(colnums_to_keep_in_probs))
        probs3 = cbind(probs2[, colnums_to_keep_in_probs_in_order], 
                       other)
        probs3 = probs3[nodes, ]
        newcols = c(colors_list_for_states[colnums_to_keep_in_probs_in_order], 
                    "white")
        nodelabels(pie = probs3, node = nodes, piecol = newcols, 
                   cex = statecex)
      }
      else {
        nodelabels(pie = relprobs_matrix_for_internal_states, 
                   node = nodes, piecol = colors_list_for_states, 
                   cex = statecex)
      }
      if (tipboxes_TF == TRUE) {
        tiplabels(text = MLstates[tips], tip = tips, 
                  bg = cols_byNode[tips], cex = pie_tip_statecex, 
                  adj = tiplabel_adj)
      }
    }
    par(fg = "black")
  }
  if (skiptree != TRUE) {
    if (titlecex > 0) {
      par(cex.main = titlecex)
      title(analysis_titletxt, line = title_offset)
      if (!is.null(dej_params_row)) {
        title(titletxt2, line = title_offset)
      }
    }
  }
  if (plotsplits == TRUE) {
    if (cornercoords_loc == "manual") {
      stoptxt = cat("\nNOTE: To plot splits, this function needs to access the function 'plot_phylo3_nodecoords'.\n", 
                    "The function is modified from an APE function, and cannot be directly included in the package,\n", 
                    "due to some C code that does not meet CRAN standards. To solve this, give plot_BioGeoBEARS_results\n", 
                    "a 'cornercoords_loc' string that gives the directory of plot_phylo3_nodecoords.R.  Typically this\n", 
                    "can be found via: ", "tmp=np(system.file(\"extdata/a_scripts\", package=\"BioGeoBEARS\"))\n", 
                    "then: list.files(tmp); print(tmp)\n", sep = "")
      plotsplits = FALSE
    }
  }
  if (plotsplits == TRUE) {
    coords_df = corner_coords(tr, tmplocation = cornercoords_loc, 
                              root.edge = root.edge)
    relprobs_matrix = left_ML_marginals_by_node
    if (plotwhat == "text") {
      MLprobs = get_ML_probs(relprobs_matrix)
      MLstates = get_ML_states_from_relprobs(relprobs_matrix, 
                                             statenames, returnwhat = "states", if_ties = if_ties)
      cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, 
                                        colors_list_for_states, MLstates)
      par(fg = tmp_fg)
      if (skiplabels == FALSE) {
        cornerlabels(text = MLstates, coords = coords_df$leftcorns, 
                     bg = cols_byNode, cex = splitcex)
      }
      par(fg = "black")
    }
    if (plotwhat == "pie") {
      par(fg = tmp_fg)
      cornerpies(pievals = relprobs_matrix, coords = coords_df$leftcorns, 
                 piecol = colors_list_for_states, cex = splitcex)
      par(fg = "black")
    }
    relprobs_matrix = right_ML_marginals_by_node
    if (plotwhat == "text") {
      MLprobs = get_ML_probs(relprobs_matrix)
      MLstates = get_ML_states_from_relprobs(relprobs_matrix, 
                                             statenames, returnwhat = "states", if_ties = if_ties)
      cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, 
                                        colors_list_for_states, MLstates)
      par(fg = tmp_fg)
      if (skiplabels == FALSE) {
        cornerlabels(text = MLstates, coords = coords_df$rightcorns, 
                     bg = cols_byNode, cex = splitcex)
      }
      par(fg = "black")
    }
    if (plotwhat == "pie") {
      par(fg = tmp_fg)
      cornerpies(pievals = relprobs_matrix, coords = coords_df$rightcorns, 
                 piecol = colors_list_for_states, cex = splitcex)
      par(fg = "black")
    }
  }
  if (((is.null(BioGeoBEARS_run_object$timeperiods) == FALSE)) && 
      (plot_stratum_lines == TRUE)) {
    timeperiods = BioGeoBEARS_run_object$timeperiods
    line_positions_on_plot = add_statum_boundaries_to_phylo_plot(tr, 
                                                                 timeperiods = timeperiods, lwd = axis.lwd, lty = "dashed", col = "gray50", 
                                                                 plotlines = TRUE)
  }
  param_ests = matrix(data = param_ests, nrow = 1)
  param_ests = adf2(param_ests)
  param_ests = dfnums_to_numeric(param_ests)
  names(param_ests) = c("LnL", "nparams", param_names)
  return(param_ests)
}

# adjusted function from BioGeoBears
add_statum_boundaries_to_phylo_plot <- function (tr, timeperiods = 1, lwd = 1,  lty = "dashed", col = "gray50", 
                                                 plotlines = TRUE) 
{
  SETUP = "\n\t# Loading the default tree\n\ttrfn = np(paste(addslash(extdata_dir), \"Psychotria_5.2.newick\", sep=\"\"))\n\ttr = read.tree(trfn)\n\t\n\t# Get the tree coordinates (APE 5.0 or higher), i.e. the x and y of each node.\n\ttrcoords = plot_phylo3_nodecoords_APE5(tr, plot=FALSE, root.edge=TRUE)\n\t\n\t# Set reasonable x-limits (unlike the defaults on APE5.0)\n\txlims = c(min(trcoords$xx), 1.42*max(trcoords$xx))\n\t\n\t# Plot the tree\n\ttrplot = plot(tr, cex=1, x.lim=xlims); axisPhylo()\n\t\n\t# Add the stratum boundaries\n\tadd_statum_boundaries_to_phylo_plot(tr, timeperiods=1, lty=\"dashed\", col=\"gray50\", plotlines=TRUE)\n\t"
  ntips = length(tr$tip.label)
  tr_table = prt(tr, printflag = FALSE)
  tr_height = tr_table$time_bp[ntips + 1]
  line_positions_on_plot = tr_height - timeperiods
  if (plotlines == TRUE) {
    abline(v = line_positions_on_plot, lwd = lwd, lty = lty, col = col)
  }
  return(line_positions_on_plot)
}

# make legend for BioGeoBears BSR
make_legend <- function(res, states = "all", names = NULL, max_range_size = 3, plot_NULL = FALSE, 
                        ncol = NULL, cex = 0.8, blank = FALSE,
                        save = TRUE, filename = "colors_legend.pdf", width = 6, height = 6) {
  
  
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=res$inputs$geogfn)
  areanames = names(tipranges@df)
  
  states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=plot_NULL)
  
  statenames = areas_list_to_states_list_new(areas=areanames, maxareas=max_range_size, include_null_range=plot_NULL, split_ABC=FALSE)
  relprobs_matrix = res$ML_marginal_prob_each_state_at_branch_top_AT_node
  MLprobs = get_ML_probs(relprobs_matrix)
  MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties="takefirst")
  
  colors_matrix = get_colors_for_numareas(length(areanames))
  colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index)
  
  possible_ranges_list_txt = areas_list_to_states_list_new(areas=areanames,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=plot_NULL)
  cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)
  
  legend_ncol=ncol
  legend_cex=cex
  
  if (is.numeric(states)) {
    possible_ranges_list_txt <- possible_ranges_list_txt[states]
    if (!is.null(names)){
      possible_ranges_list_txt <- lapply(1:length(states), function(i) paste(possible_ranges_list_txt[[i]], ":", names[i]))
    }
    colors_list_for_states <- colors_list_for_states[states]
  }
  colors_legend(possible_ranges_list_txt, colors_list_for_states, legend_ncol=legend_ncol, legend_cex=legend_cex)
  if (save){
    pdf(filename, width=width, height=height)
    colors_legend(possible_ranges_list_txt, colors_list_for_states, 
                  legend_ncol=legend_ncol, legend_cex=legend_cex,
                  make_blank_plot_first = blank)
    dev.off()
  }  
  
}

# Modified EIC from mvMORPH to avoid getting error when having a warning as a listed element
EIC <- function(object, nboot=100L, nbcores=1L, ...){
  
  # retrieve arguments
  args <- list(...)
  if(is.null(args[["eigSqm"]])) eigSqm <- TRUE else eigSqm <- args$eigSqm
  if(is.null(args[["restricted"]])) restricted <- FALSE else restricted <- args$restricted
  if(is.null(args[["REML"]])) args$forceREML <- FALSE else args$forceREML <- args$REML
  
  # retrieve data to simulate bootstrap samples
  beta <- object$coefficients
  if(eigSqm){ # to follow the scheme in RPANDA
    sqM1 <- mvMORPH:::.sqM1(object$corrSt$phy)
    if(!is.null(object$corrSt$diagWeight)){
      w <- 1/object$corrSt$diagWeight
      Y <- crossprod(sqM1, matrix(w*object$variables$Y, nrow=object$dims$n))
      X <- crossprod(sqM1, matrix(w*object$variables$X, nrow=object$dims$n))
    }else{
      X <- crossprod(sqM1, object$variables$X)
      Y <- crossprod(sqM1, object$variables$Y)
    }
    residuals <- Y - X%*%beta
  }else{
    residuals <- residuals(object, type="normalized")
    X <- object$corrSt$X
    Y <- object$corrSt$Y
  }
  
  N = nrow(Y)
  p = object$dims$p
  if(object$REML & args$forceREML==TRUE) ndimCov = object$dims$n - object$dims$m else ndimCov = object$dims$n
  tuning <- object$tuning
  target <- object$target
  penalty <- object$penalty
  if(is.null(object$corrSt$diagWeight)){
    diagWeight <- 1; is_weight = FALSE
  }else{
    diagWeight <- object$corrSt$diagWeight; is_weight = TRUE
    diagWeightInv <- 1/diagWeight
  }
  Dsqrt <- mvMORPH:::.pruning_general(object$corrSt$phy, trans=FALSE, inv=FALSE)$sqrtM # return warning message if n-ultrametric tree is used with OU?
  # TODO (change to allow n-ultrametric and OU) > just need to standardize the data by the weights
  # if(object$model=="OU" & !is.ultrametric(object$variables$tree)) stop("The EIC method does not handle yet non-ultrametric trees with OU processes")
  
  DsqrtInv <- mvMORPH:::.pruning_general(object$corrSt$phy, trans=FALSE, inv=TRUE)$sqrtM
  modelPerm <- object$call
  modelPerm$grid.search <- quote(FALSE)
  modelPerm$start <- quote(object$opt$par)
  
  # Mean and residuals for the model
  MeanNull <- object$variables$X%*%beta
  
  # Estimate the bias term
  D1 <- function(objectBoot, objectFit, ndimCov, p, sqM, Ccov2){ # LL(Y*|param*) - LL(Y*| param)
    
    # Y*|param*
    residualsBoot <- residuals(objectBoot, type="normalized")
    
    # For boot "i" LL1(Y*|param*)
    if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov1 <- as.numeric(objectBoot$corrSt$det - determinant(crossprod(objectBoot$corrSt$X))$modulus + objectBoot$corrSt$const) else Ccov1 <- as.numeric(objectBoot$corrSt$det)
    Gi1 <- try(chol(objectBoot$sigma$Pinv), silent=TRUE)
    if(inherits(Gi1, 'try-error')) return("error")
    quadprod <- sum(backsolve(Gi1, t(residualsBoot), transpose = TRUE)^2)
    detValue <- sum(2*log(diag(Gi1)))
    llik1 <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov1 + ndimCov*detValue + quadprod)
    
    # Y*|param
    #if(!restricted) residualsBoot <- objectBoot$corrSt$Y - objectBoot$corrSt$X%*%objectFit$coefficients # does not account for the phylo model of the original fit
    if(!restricted){
      if(is_weight){
        residualsBoot <- crossprod(sqM, (objectBoot$variables$Y - objectBoot$variables$X%*%objectFit$coefficients)*diagWeightInv)
      }else{
        residualsBoot <- crossprod(sqM, objectBoot$variables$Y - objectBoot$variables$X%*%objectFit$coefficients)
        
      }
    }
    
    # For boot "i" LL2(Y*|param)
    # if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov2 <- as.numeric(objectFit$corrSt$det - determinant(crossprod(objectFit$corrSt$X))$modulus + objectFit$corrSt$const) else Ccov2 <- as.numeric(objectFit$corrSt$det)
    Gi2 <- try(chol(objectFit$sigma$Pinv), silent=TRUE)
    if(inherits(Gi2, 'try-error')) return("error")
    quadprod <- sum(backsolve(Gi2, t(residualsBoot), transpose = TRUE)^2)
    detValue <- sum(2*log(diag(Gi2)))
    llik2 <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov2 + ndimCov*detValue + quadprod)
    
    # Return the difference in LL for D1
    return(llik1 - llik2)
  }
  D3 <- function(objectBoot, objectFit, loglik, ndimCov, p){ # LL(Y|param) - LL(Y| param*)
    
    # Y|param*
    if(!restricted) {
      sqM_temp <- mvMORPH:::.pruning_general(objectBoot$corrSt$phy, trans=FALSE, inv=TRUE)$sqrtM
      if(is_weight){
        residualsBoot <- try(crossprod(sqM_temp, (objectFit$variables$Y - objectFit$variables$X%*%objectBoot$coefficients)/objectBoot$corrSt$diagWeight), silent=TRUE)
      } else {
        residualsBoot <- try(crossprod(sqM_temp, objectFit$variables$Y - objectFit$variables$X%*%objectBoot$coefficients), silent=TRUE)
        
      }
    }else{ residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectFit$coefficients}
    
    #if(!restricted) residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectBoot$coefficients
    #else residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectFit$coefficients
    
    # For boot "i" LL2(Y|param*)
    if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov1 <- as.numeric(objectBoot$corrSt$det - determinant(crossprod(objectBoot$corrSt$X))$modulus + objectBoot$corrSt$const) else Ccov1 <- as.numeric(objectBoot$corrSt$det)
    Gi1 <- try(chol(objectBoot$sigma$Pinv), silent=TRUE)
    if(inherits(Gi1, 'try-error')) return("error")
    quadprod <- sum(backsolve(Gi1, t(residualsBoot), transpose = TRUE)^2)
    detValue <- sum(2*log(diag(Gi1)))
    llik2 <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov1 + ndimCov*detValue + quadprod)
    
    # Return the difference in LL for D1
    return(loglik - llik2)
  }
  
  # Estimate EIC: LL+bias
  
  # Maximum Likelihood
  if(object$REML==TRUE & args$forceREML==FALSE) Ccov <- as.numeric(object$corrSt$det - determinant(crossprod(object$corrSt$X))$modulus + object$corrSt$const) else Ccov <- as.numeric(object$corrSt$det)
  Gi <- try(chol(object$sigma$Pinv), silent=TRUE)
  if(inherits(Gi, 'try-error')) return("error")
  quadprod <- sum(backsolve(Gi, t(residuals), transpose = TRUE)^2)
  detValue <- sum(2*log(diag(Gi)))
  llik <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*detValue + quadprod)
  
  # Estimate parameters on bootstrap samples
  bias <- pbmcapply:::pbmcmapply(function(i){
    
    # generate bootstrap sample
    Yp <- MeanNull + Dsqrt%*%(residuals[sample(N, replace=TRUE),])*diagWeight # sampling with replacement for bootstrap
    rownames(Yp) <- rownames(object$variables$Y)
    
    modelPerm$response <- quote(Yp);
    estimModelNull <- eval(modelPerm);
    d1res <- D1(objectBoot=estimModelNull, objectFit=object, ndimCov=ndimCov, p=p, sqM=DsqrtInv, Ccov2=Ccov)
    d3res <- D3(objectBoot=estimModelNull, objectFit=object, loglik=llik, ndimCov=ndimCov, p=p)
    d1res+d3res
    
  }, 1:nboot, mc.cores = getOption("mc.cores", nbcores))
  
  # check for errors first?
  bias <- mvMORPH:::.check_samples(bias)
  if (length(bias) == 2){
    bias <- bias[[1]]
  }
  nboot_eff <- length(bias)
  # compute the EIC
  pboot <- mean(bias)
  EIC <- -2*llik + 2*pboot
  
  # standard-error
  se <- sd(bias)/sqrt(nboot_eff)
  
  # concatenate the results
  results <- list(EIC=EIC, bias=bias, LogLikelihood=llik, se=se, p=p, n=N)
  class(results) <- c("eic.mvgls","eic")
  
  return(results)
}

#### from Theo: Create Geo object from BioGeoBEARS
CreateCombGeo <- function(tree, bsm, asm, nmap){
  geo_obj <- CreateGeoObject_BioGeoBEARS(tree, ana.events = bsm$RES_ana_events_tables[[nmap]], clado.events = bsm$RES_clado_events_tables[[nmap]], stratified = FALSE)
  a_map <- asm$simmap[[nmap]]
  n <- length(a_map$tip.label)
  root <- as.numeric(branching.times(a_map)[1])
  ## classify event from most ancient to most recent
  a_events <- do.call(rbind, sapply(1:length(a_map$edge.length), function(t_n){
    desc <- a_map$edge[t_n,2]
    anc <- a_map$edge[t_n,1]
    age_desc <- ifelse(desc > n, branching.times(a_map)[desc-n], 0)
    abs_time <- round(root - age_desc - cumsum(rev(a_map$maps[[t_n]])),4)
    if(anc == n + 1) abs_time[1] <- 0
    if(length(abs_time) > 1) rbind(t(sapply((length(abs_time)-1):1, function(i) c(as.numeric(abs_time[i]), as.numeric(names(abs_time)[c(i,i+1)]), t_n))),
                                   c(as.numeric(abs_time[length(abs_time)]), as.numeric(names(abs_time)[length(abs_time)]),0, t_n))
    else rbind(c(as.numeric(abs_time), as.numeric(names(abs_time)),0, t_n))
  }))
  
  ## deal with weird naming from CreateGeoObject
  nodes <- 1:(2*n-1)
  states <- numeric(2*n-1)
  i_l <- 0
  for(no in order(branching.times(a_map), decreasing = T)+n){
    desc <- a_map$edge[a_map$edge[,1] == no,2]
    for(i in 1:2){
      if(desc[i] > n){
        i_l <- i_l + 1
        names(nodes)[desc[i]] <- names(states)[desc[i]] <- sprintf(".AA%s", LETTERS[i_l])
      } else {
        names(nodes)[desc[i]] <- names(states)[desc[i]] <- a_map$tip.label[desc[i]]
      } 
    }
  }
  ## change matrices
  a_events <- a_events[order(a_events[,1]),]
  a_events <- a_events[!duplicated(a_events[,1]),]
  i_geo <- 1
  i_tot <- 1
  last_stop <- 0
  new_geo <- list()
  new_times <- c()
  
  for(i in 1:nrow(a_events)){
    if(a_events[i,3] == 0){## Cladogenesis : modify existing matrix
      anc <- a_map$edge[a_events[i,4],1]
      states[a_map$edge[a_map$edge[,1] == anc,2]] <- a_events[i,2]
      for(k in 1:(sum(round(geo_obj$times,4) < a_events[i,1] & round(geo_obj$times, 4) > last_stop)+1)){
        new_geo[[i_tot]] <- geo_obj$geography.object[[i_geo]]
        new_times[i_tot] <- geo_obj$times[i_geo]
        for(nd1 in rownames(new_geo[[i_tot]])){
          for(nd2 in rownames(new_geo[[i_tot]])){
            if(new_geo[[i_tot]][nd1, nd2] == 1){
              if(any(states[c(nd1, nd2)] == 4)) new_geo[[i_tot]][nd1, nd2] <- 1
              else if(states[nd1] == states[nd2]) new_geo[[i_tot]][nd1, nd2] <- 1
              else new_geo[[i_tot]][nd1, nd2] <- 0
            }
          }
        }
        i_tot <- i_tot + 1
        i_geo <- i_geo + 1
      }
    } else { ## Add matrix, time and span and update following
      ## first update matrices that are before event
      if(sum(round(geo_obj$times,4) < a_events[i,1] & round(geo_obj$times, 4) > last_stop) > 0){
        for(k in 1:(sum(round(geo_obj$times,4) < a_events[i,1] & round(geo_obj$times, 4) > last_stop))){
          new_geo[[i_tot]] <- geo_obj$geography.object[[i_geo]]
          new_times[i_tot] <- geo_obj$times[i_geo]
          for(nd1 in rownames(new_geo[[i_tot]])){
            for(nd2 in rownames(new_geo[[i_tot]])){
              if(new_geo[[i_tot]][nd1, nd2] == 1){
                if(any(states[c(nd1, nd2)] == 4)) new_geo[[i_tot]][nd1, nd2] <- 1
                else if(states[nd1] == states[nd2]) new_geo[[i_tot]][nd1, nd2] <- 1
                else new_geo[[i_tot]][nd1, nd2] <- 0
              }
            }
          }
          i_tot <- i_tot + 1
          i_geo <- i_geo + 1
        }
      }
      
      desc <- a_map$edge[a_events[i,4],2]
      states[desc] <- a_events[i,3]
      new_geo[[i_tot]] <- new_geo[[i_tot - 1]]
      new_times[i_tot] <- a_events[i,1]
      for(nd1 in rownames(new_geo[[i_tot]])){
        for(nd2 in rownames(new_geo[[i_tot]])){
          if(new_geo[[i_tot]][nd1, nd2] == 1){
            if(any(states[c(nd1, nd2)] == 4)) new_geo[[i_tot]][nd1, nd2] <- 1
            else if(states[nd1] == states[nd2]) new_geo[[i_tot]][nd1, nd2] <- 1
            else new_geo[[i_tot]][nd1, nd2] <- 0
          }
        }
      }
      i_tot <- i_tot + 1
    }
    last_stop <- a_events[i,1]
  }
  list(geography.object = new_geo, times = new_times, spans = diff(new_times))
}

# FROM: cmt2/bgb_to_revgadgets.R
bgb_to_revgadgets <- function(results_path, geo_data_path, tree_path, area_names = NULL) {
  # load biogeobears results object
  res <- readRDS(results_path)
  # change data directories in results object
  res[["inputs"]]$geogfn <- geo_data_path
  res[["inputs"]]$trfn <- tree_path
  
  ##### Process data for plotting ##### 
  
  # read in tree separately
  
  tree <- RevGadgets::readTrees(paths = res[["inputs"]]$trfn)
  states <- res$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
  
  # create a dataframe with results in revgadgets compliant format
  rev_data <- data.frame(matrix(nrow = nrow(res$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS),
                                ncol = 15))
  colnames(rev_data) <- c("end_state_1", "end_state_2", "end_state_3", 
                          "end_state_1_pp", "end_state_2_pp", "end_state_3_pp", 
                          "end_state_other_pp",
                          "start_state_1", "start_state_2", "start_state_3", 
                          "start_state_1_pp", "start_state_2_pp", "start_state_3_pp", 
                          "start_state_other_pp",
                          "node")
  # get end states
  for (i in 1:nrow(res$ML_marginal_prob_each_state_at_branch_top_AT_node)) {
    row <- res$ML_marginal_prob_each_state_at_branch_top_AT_node[i,]
    rev_data[i, 1] <- order(row,decreasing=T)[1]
    rev_data[i, 2] <- order(row,decreasing=T)[2]
    rev_data[i, 3] <- order(row,decreasing=T)[3]
    rev_data[i, 4] <- row[order(row,decreasing=T)[1]]
    rev_data[i, 5] <- row[order(row,decreasing=T)[2]]
    rev_data[i, 6] <- row[order(row,decreasing=T)[3]]
    rev_data[i, 7] <- sum(row[order(row,decreasing=T)[4:length(row)]]) 
  }
  # get start states
  for (i in 1:nrow(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)) {
    row <- res$ML_marginal_prob_each_state_at_branch_bottom_below_node[i,]
    rev_data[i, 8] <- order(row,decreasing=T)[1]
    rev_data[i, 9] <- order(row,decreasing=T)[2]
    rev_data[i, 10] <- order(row,decreasing=T)[3]
    rev_data[i, 11] <- row[order(row,decreasing=T)[1]]
    rev_data[i, 12] <- row[order(row,decreasing=T)[2]]
    rev_data[i, 13] <- row[order(row,decreasing=T)[3]]
    rev_data[i, 14] <- sum(row[order(row,decreasing=T)[4:length(row)]]) 
  }
  rev_data$node <- 1:nrow(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)
  
  # make better labels 
  tipranges <- getranges_from_LagrangePHYLIP(res[["inputs"]]$geogfn)
  
  geo <- res$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
  geo_num <- unlist(lapply(lapply(geo, as.character), paste0, collapse ="_"))
  if (is.null(area_names)) {
    area_names <- colnames(tipranges@df)
  }
  if (length(area_names) != length(colnames(tipranges@df))) stop("Number of specified area names is incorrect. Check your geo data file.")
  number_codes <- 0:(length(area_names)-1)
  geo_letters <- geo_num
  
  recodes <- paste0(paste0(number_codes, " = '", area_names, "'"), collapse = "; ")
  for (i in 1:length(geo_letters)) {
    code <- geo_letters[i]
    code_split <- unlist(strsplit(code, "_"))
    geo_letters[i] <- paste0(car::recode(code_split, recodes), collapse = "")
  }
  
  #area_names <- rev(area_names)
  label_dict <- data.frame(lab_num_short = 1:length(geo),
                           lab_num_long = geo_num,
                           lab_letters = geo_letters)
  
  # replace short number codes with letters 
  not_state_cols <- c(grep("_pp", colnames(rev_data)),
                      grep("node", colnames(rev_data)))
  state_cols <- c(1:ncol(rev_data))[!c(1:ncol(rev_data)) %in% not_state_cols]
  for (i in state_cols) { # loop through by column indices for the state columns
    col <- as.character(rev_data[,i])
    for (j in 1:length(col)) { # loop through each item in the column and replace with letters
      col[j] <- label_dict$lab_letters[which(label_dict$lab_num_short == col[j])]
    }
    rev_data[,i] <- col
  }
  
  #change "NA" to NA 
  rev_data %>%
    naniar::replace_with_na_all(condition = ~.x == "NA") -> rev_data
  # change any NAs in PP columns to 0 
  pp_cols <- grep("_pp", colnames(rev_data))
  for (p in pp_cols) { rev_data[ ,p][is.na(rev_data[ ,p])] <- 0 }
  
  # make treedata object (combine data and tree)
  tibble::as_tibble(tree[[1]][[1]]) %>%
    full_join(rev_data, by = 'node') %>%
    as.treedata() -> rev_treedata
  #add list of states
  attributes(rev_treedata)$state_labels <- as.character(na.omit(unique(unlist(rev_data[,c(1:3, 8:10)]))))
  
  return(rev_treedata)
}

# image transformation based on Procrustes analysis
imageTransformation <- function (sampleList, landList, adjustCoords = F, transformRef = "meanshape", 
                                 crop = FALSE, cropOffset = c(0, 0, 0, 0), res = 300, keep.ASP = T, aspect.ratio = NULL,
                                 drop = NULL, removebg.by = c(FALSE, "color", "landmarks"), 
                                 smooth = FALSE, rescale = F, resampleFactor = NULL, transformType = "tps", 
                                 focal = FALSE, sigma = 3, interpolate = NULL, bgcol = NULL, bg.offset = NULL, 
                                 plot = FALSE, save = FALSE, dir = "./", overwrite = FALSE) 
{
  removebg.by = removebg.by[1]
  rasterList <- list()
  if (length(sampleList) != length(landList)) {
    stop("sampleList is not of the same length as lanArray")
  }
  if (any(!names(landList) %in% names(sampleList))) {
    stop("sampleList names do not match names in lanArray")
  }
  landList <- landList[names(sampleList)]
  lanArray <- patternize::lanArray(landList, adjustCoords = adjustCoords, 
                                   sampleList)
  if (is.matrix(transformRef)) {
    refShape <- transformRef
  }
  if (!is.matrix(transformRef)) {
    if (transformRef == "meanshape") {
      invisible(utils::capture.output(transformed <- Morpho::procSym(lanArray)))
      refShape <- transformed$mshape
    }
    if (transformRef %in% names(landList)) {
      e <- which(names(landList) == transformRef)
      refShape <- lanArray[, , e]
    }
  }
  if (removebg.by == "landmarks") {
    if (is.null(drop)) {
      drop = 1:nrow(refShape)
    }
    sp.ref = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(refShape[-drop, 
    ])), 1)))
    if (is.numeric(smooth)) {
      sp.ref <- smoothr::smooth(sp.ref, method = "ksmooth", 
                                smooth = smooth)
    }
  }
  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    if (!overwrite) {
      files <- gsub("\\.tif", "", list.files(dir))
      idx <- which(names(sampleList) %in% files)
      if (length(idx) > 0) {
        message(length(idx), " images already present in the directory as transformed and thus, removed from the task. Modify parameters if is not the case.")
        sampleList <- sampleList[-idx]
        landList <- landList[-idx]
      }
    }
  }
  
  for (n in 1:length(sampleList)) {
    image <- sampleList[[n]]
    extRasterOr <- raster::extent(image)
    if (!is.null(resampleFactor)) {
      image <- patternize::redRes(image, resampleFactor)
    }
    if (crop) {
      landm <- lanArray[, , n]
      extRaster <- raster::extent(min(landm[, 1]) - min(landm[,1]) * cropOffset[1]/100, 
                                  max(landm[, 1]) + max(landm[,1]) * cropOffset[2]/100,
                                  min(landm[, 2]) - min(landm[, 2]) * cropOffset[3]/100, 
                                  max(landm[, 2]) + max(landm[, 2]) * cropOffset[4]/100)
      imageC <- raster::crop(image, extRaster)
      y <- raster::raster(ncol = dim(image)[2], nrow = dim(image)[1])
      raster::extent(y) <- extRasterOr
      image <- raster::resample(imageC, y)
    }
    if (focal) {
      gf <- raster::focalWeight(image, sigma, "Gauss")
      rrr1 <- raster::focal(image[[1]], gf)
      rrr2 <- raster::focal(image[[2]], gf)
      rrr3 <- raster::focal(image[[3]], gf)
      image <- raster::stack(rrr1, rrr2, rrr3)
    }
    mapDF <- raster::as.data.frame(image, xy = TRUE)
    invisible(utils::capture.output(transMatrix <- Morpho::computeTransform(refShape, 
                                                                            as.matrix(lanArray[, , n]), type = "tps")))
    invisible(utils::capture.output(mapTransformed <- Morpho::applyTransform(as.matrix(mapDF[, 
                                                                                             1:2]), transMatrix)))
    ASP = 1
    if (keep.ASP) {
      ext = as.vector(raster::extent(image))
      if (!is.null(aspect.ratio)) { ASP = aspect.ratio }
      else {ASP = (ext[2] - ext[1])/(ext[4] - ext[3])}
    }
    rRe <- raster::raster(ncol = res * ASP, nrow = res)
    e = raster::extent(min(refShape[, 1]) - 3 * max(refShape[, 
                                                             1]) * cropOffset[3]/100, max(refShape[, 1]) + 3 * 
                         max(refShape[, 1]) * cropOffset[4]/100, min(refShape[, 
                                                                              2]) - 3 * max(refShape[, 2]) * cropOffset[1]/100, 
                       max(refShape[, 2]) + 3 * max(refShape[, 2]) * cropOffset[2]/100)
    margin = c((e[2] - e[1]) * 0.05, (e[4] - e[3]) * 0.05)
    raster::extent(rRe) = raster::extent(e[1] - margin[1], 
                                         e[2] + margin[1], e[3] - margin[2], e[4] + margin[2])
    imgTransformed <- raster::stack(sapply(names(image), 
                                           function(x) raster::rasterize(mapTransformed, field = mapDF[, 
                                                                                                       x], rRe)))
    if (removebg.by == "landmarks") {
      message("Removing background based on landmarks reference...")
      imgTransformed <- raster::mask(imgTransformed, sp.ref)
    }
    else if (removebg.by == "color") {
      message("Removing background based on color...")
      imgTransformed <- removebg(imgTransformed, bgcol = bgcol, 
                                 bg.offset = bg.offset, plot = F)
    }
    if (is.matrix(transformRef)) {
      imgTransformed = raster::flip(raster::flip(imgTransformed, 
                                                 "x"), "y")
    }
    else if (!is.matrix(transformRef)) {
      if (transformRef == "meanshape") {
        imgTransformed = raster::flip(raster::flip(imgTransformed, 
                                                   "x"), "y")
      }
    }
    if (rescale) {
      rRe <- raster::raster(nrow = dim(image)[1], ncol = dim(image)[2])
      raster::crs(rRe) = raster::crs(image)
      raster::extent(rRe) <- raster::extent(image)
      raster::extent(imgTransformed) = raster::extent(image)
      imgTransformed = raster::stack(raster::resample(imgTransformed, 
                                                      rRe, method = "ngb"))
    }
    if (is.numeric(interpolate)) {
      imgTransformed = sp::disaggregate(imgTransformed, 
                                        fact = interpolate, method = "bilinear")
      imgTransformed = raster::resample(imgTransformed, 
                                        rRe, method = "ngb")
    }
    imgTransformed[imgTransformed[] < 0] = 0
    if (plot == "result") {
      if (nlayers(imgTransformed) == 3) {
        raster::plotRGB(imgTransformed)
      }
      if (nlayers(imgTransformed) == 1) {
        plot(imgTransformed)
      }
    }
    if (plot == "compare") {
      graphics::par(mfrow = c(1, 2))
      if (nlayers(imgTransformed) == 3) {
        raster::plotRGB(image)
        text(raster::extent(image)[2], raster::extent(image)[3], 
             "original", adj = c(1, 0))
        raster::plotRGB(imgTransformed, asp = 1)
        text(raster::extent(imgTransformed)[2], raster::extent(imgTransformed)[3], 
             "transformed", adj = c(1, 0))
      }
      if (nlayers(imgTransformed) == 1) {
        plot(image)
        text(raster::extent(image)[2], raster::extent(image)[3], 
             "original", adj = c(1, 0))
        plot(imgTransformed)
        text(raster::extent(imgTransformed)[2], raster::extent(imgTransformed)[3], 
             "transformed", adj = c(1, 0))
      }
    }
    if (save) {
      raster::writeRaster(imgTransformed, filename = file.path(dir, 
                                                               paste0(names(landList)[n], ".tif")), overwrite = overwrite)
    }
    rasterList[[names(landList)[n]]] <- imgTransformed
    print(paste("sample", names(landList)[n], "transformation done and added to rasterList", 
                sep = " "))
  }
  return(rasterList)
}

# plot images along coordinates
plotImages <- function(x, y, images, width = 0.1, height = NULL, interpolate = FALSE,
                       names = NULL, cex = 1, pos = 1, adj = 1,
                       cols  = c("red", "grey90", "blue"), angle = NULL, flip = TRUE, ...){
  
  cols = grDevices::colorRampPalette(cols)(n=100)
  stopifnot(length(x) == length(y))
  if (is.null(height)) {
    asp <- sapply(images, function(i) abs(diff(extent(i)[1:2])/diff(extent(i)[3:4])))
    height <- width / asp
  }
  if (!is.null(angle)) {
    if(length(angle) < length(images)){
      angle <- rep(angle[1], length(images))
    }
  } else {
    angle <- rep(0, length(images))
  }
  if(length(images) < length(x)){
    images <- replicate(length(x), images, simplify=FALSE)
  }
  if(!is.null(names) && length(names) < length(images)){
    names <- rep(names, length(images))
  }
  if(length(pos) < length(images)){
    pos <- rep(pos, length(images))
  }
  if(length(width) < length(x)){
    width <- rep(width, length(images))
  }
  if(length(height) < length(x)){
    height <- rep(height, length(images))
  }
  width = width+width*diff(range(x))
  height = height+height*diff(range(y))
  for (ii in seq_along(x)){
    if (flip){
      if (angle[ii] > 90 | angle[ii] < -90 ){
        images[[ii]] <- raster::flip(images[[ii]], "y")
      }
    }
    if (class(images[[ii]]) %in% c("RasterStack", "RasterBrick") && dim(images[[ii]])[3] == 3) {
      images[[ii]] = sapply(1:3, function(i) images[[ii]][[i]]/255)
      images[[ii]] = as.array(raster::stack(append(images[[ii]], images[[ii]][[1]])))
      images[[ii]][,,4][!is.na(images[[ii]][,,4])] = 1
      for (i in 1:4) {
        images[[ii]][,,i][is.na(images[[ii]][,,i])] = 0
        images[[ii]][,,i][images[[ii]][,,i] >1] = 1
        images[[ii]][,,i][images[[ii]][,,i] <0] = 0
      }
    }
    if (is.array(images[[ii]]) && dim(images[[ii]])[3] == 4){
      graphics::rasterImage(images[[ii]], xleft=x[ii] - 0.5*width[ii],
                            ybottom= y[ii] - 0.5*height[ii],
                            xright=x[ii] + 0.5*width[ii],
                            ytop= y[ii] + 0.5*height[ii], interpolate=interpolate, angle = angle[[ii]], ...)
    }
    if (class(images[[ii]]) == "RasterLayer") {
      e = as.vector( raster::extent(images[[ii]]))
      ratio = (e[2]-e[1])/(e[4]-e[3])
      raster::extent(images[[ii]]) =  c(x[ii] - 0.5*width[ii], x[ii] + 0.5*width[ii], y[ii] - 0.5*height[ii]*ratio, y[ii] + 0.5*height[ii]*ratio)
      raster::image(images[[ii]], interpolate=interpolate, add = T, legend = F, col = cols, ...)
    }
    if(!is.null(names)){
      graphics::text(x[ii], y[ii] - 0.5*height[ii], names[ii], pos = pos, adj = adj)
    }
  }
}

# get the angle relative to a center from a point
get_angle <- function(center_x, center_y, point_x, point_y) {
  dx <- point_x - center_x
  dy <- point_y - center_y
  angle <- atan2(dy, dx) * (180 / pi)
  return(angle)
}

# slight modification from phyloch axisGeo function
axisGeo <- function (GTS, tip.time = 0, unit = c("epoch", "period"), ages = TRUE, offset = -0.2,
                     cex = 1, col = "white", texcol = "black", gridty = 0, gridcol = "black") 
{
  adjustCex <- function(space, string, cex) {
    while (strwidth(string, cex = cex) >= space & cex > 0.001) cex <- cex - 
        0.001
    cex
  }
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  ntips <- lastPP$Ntip
  root <- ntips + 1
  if (lastPP$direction == "rightwards") 
    maxage <- max(lastPP$xx) + tip.time
  if (lastPP$direction == "upwards") 
    maxage <- max(lastPP$yy) + tip.time
  gts <- GTS
  maid <- grep("MA", names(gts))
  gts <- cbind(gts[, 1:maid], c(0, head(gts[, maid], -1)), 
               gts[, ((maid + 1):dim(gts)[2])])
  names(gts)[maid:(maid + 1)] <- c("fromMA", "toMA")
  if (sum(gts[1, maid:(maid + 1)]) == 0) 
    gts <- gts[-1, ]
  ind <- which(gts$fromMA <= maxage & gts$toMA >= tip.time)
  gts <- gts[c(ind, max(ind) + 1), ]
  gts$toMA[1] <- tip.time
  gts$fromMA[dim(gts)[1]] <- maxage
  par(xpd = NA)
  plotGeo <- function(gts, unit, yy) {
    id <- which(names(gts) %in% c(unit, "fromMA", "toMA"))
    gts <- gts[id]
    names(gts) <- c("unit", "from", "to")
    stages <- unique(gts$unit)
    if (length(col) == 1) 
      col <- rep(col, 2)
    col1 <- rep(col, length(stages))
    col1 <- head(col1, length(col1)/2)
    if (length(texcol) == 1) 
      texcol <- rep(texcol, 2)
    col2 <- rep(texcol, length(stages))
    col2 <- head(col2, length(col2)/2)
    xgrid <- NULL
    for (i in seq(along = stages)) {
      cat("\nStage ", i, ": ", stages[i], sep = "")
      from <- maxage - max(gts[gts$unit == stages[i], 2])
      to <- maxage - min(gts[gts$unit == stages[i], 3])
      rect(from, yy[1], to, yy[2], col = col1[i], border = "black", 
           lwd = 0.5)
      xgrid <- c(xgrid, from, to)
      en <- as.character(stages[i])
      if ((to - from) > strwidth(en, cex = cex)) 
        text(mean(c(from, to)), mean(yy), en, cex = cex, 
             col = col2[4])
      else {
        thiscex <- adjustCex(to - from, en, cex) * 0.95
        if (2 * thiscex >= cex) 
          text(mean(c(from, to)), mean(yy), en, cex = thiscex, 
               col = col2[i])
        else {
          while (nchar(en) > 0 & strwidth(en, cex = cex) >= 
                 (to - from)) en <- paste(head(unlist(strsplit(en, 
                                                               "")), -1), collapse = "")
          if (nchar(en) > 1) 
            en <- paste(paste(head(unlist(strsplit(en, 
                                                   "")), -1), collapse = ""), ".", sep = "")
          text(mean(c(from, to)), mean(yy), en, cex = cex, 
               col = col2[i])
        }
      }
    }
    xgrid
  }
  plotGeoToLeft <- function(gts, unit, yy) {
    id <- which(names(gts) %in% c(unit, "fromMA", "toMA"))
    gts <- gts[id]
    names(gts) <- c("unit", "from", "to")
    stages <- unique(gts$unit)
    if (length(col) == 1) 
      col <- rep(col, 2)
    col1 <- rep(col, length(stages))
    col1 <- head(col1, length(col1)/2)
    if (length(texcol) == 1) 
      texcol <- rep(texcol, 2)
    col2 <- rep(texcol, length(stages))
    col2 <- head(col2, length(col2)/2)
    xgrid <- NULL
    for (i in seq(along = stages)) {
      from <- maxage - max(gts[gts$unit == stages[i], 2])
      to <- maxage - min(gts[gts$unit == stages[i], 3])
      rect(yy[2], from, yy[1], to, col = col1[i], border = "black", 
           lwd = 0.5)
      xgrid <- c(xgrid, from, to)
      en <- as.character(stages[i])
      yxr <- (max(lastPP$y.lim) - min(lastPP$y.lim))/(max(lastPP$x.lim) - 
                                                        min(lastPP$x.lim)) * 1.5
      if ((to - from) > strwidth(en, cex = cex * yxr)) 
        text(mean(yy), mean(c(from, to)), en, cex = cex, 
             col = col2[4], srt = 90)
      else {
        asp <- (to - from)/yxr
        thiscex <- adjustCex(asp, en, cex) * 0.95
        if (1.5 * thiscex >= cex) 
          text(mean(yy), mean(c(from, to)), en, cex = thiscex, 
               col = col2[i], srt = 90)
        else {
          while (nchar(en) > 0 & strwidth(en, cex = cex * 
                                          yxr) >= (to - from)) en <- paste(head(unlist(strsplit(en, 
                                                                                                "")), -1), collapse = "")
          if (nchar(en) > 1) 
            en <- paste(paste(head(unlist(strsplit(en, 
                                                   "")), -1), collapse = ""), ".", sep = "")
          text(mean(yy), mean(c(from, to)), en, cex = cex, 
               col = col2[i], srt = 90)
        }
      }
    }
    xgrid
  }
  bh <- -strheight("Ap", cex = cex) * 1.5
  if (lastPP$direction == "rightwards") {
    if (ages) 
      yy <- c(bh, 2 * bh)
    else yy <- c(0, bh)
    for (j in seq_along(unit)) {
      cat("\nPlot unit:", unit[j])
      if (j == 1) 
        xgrid <- plotGeo(gts, unit[j], yy)
      else plotGeo(gts, unit[j], yy)
      yy <- yy + bh
    }
  }
  if (lastPP$direction == "upwards") {
    if (ages) 
      yy <- c(bh, 2 * bh)
    else yy <- c(0, bh)
    for (j in seq(along = unit)) {
      if (j == 1) 
        xgrid <- plotGeoToLeft(gts, unit[j], yy)
      else plotGeoToLeft(gts, unit[j], yy)
      yy <- yy + bh
    }
  }
  xgrid <- unique(sort(xgrid, decreasing = TRUE))
  label <- maxage - xgrid
  id <- TRUE
  for (k in seq(along = xgrid)) {
    if (lastPP$direction == "rightwards") 
      lines(rep(xgrid[k], 2), c(0, ntips + 1), lty = gridty, 
            col = gridcol)
    if (lastPP$direction == "upwards") 
      lines(c(0, ntips + 1), rep(xgrid[k], 2), lty = gridty, 
            col = gridcol)
    if (k < length(xgrid)) {
      spneeded <- strwidth(label[k], cex = cex * 0.8)/2
      spavailable <- xgrid[k] - xgrid[k + 1]
      if (spavailable < spneeded * 1.5) 
        id <- c(id, FALSE)
      else id <- c(id, TRUE)
    }
  }
  id <- c(id, TRUE)
  if (ages) {
    xgrid <- xgrid[id]
    label <- label[id]
    if (lastPP$direction == "rightwards") 
      text(xgrid, offset, round(label, digits = 1), cex = cex * 
             0.8)
    if (lastPP$direction == "upwards") 
      text(offset, xgrid, round(label, digits = 1), cex = cex * 
             0.8)
  }
  par(xpd = FALSE)
}

# preprocess images to standardize extent, dimensions and resolution
preprocessImages <- function(imgList, res = NULL, interpolate = FALSE) {
  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("Please install the 'raster' package.")
  }
  
  # Ensure all are brick
  if (!all(sapply(imgList, inherits, c("RasterBrick", "RasterStack")))) {
    stop("All elements in imgList must be of class 'RasterBrick'. Please check your input data.")
  }
  
  # get a reference image
  ref <- imgList[[1]]
  
  if (is.null(res)) {
    # Extract resolution, extent, and dimensions from the first image
    ref_res <- raster::res(ref)
    ref_ext <- raster::extent(ref)
    ref_dim <- dim(ref)
    
    # Check all images against the reference
    for (i in seq_along(imgList)) {
      img <- imgList[[i]]
      if (!all.equal(raster::res(img), ref_res) ||
          !all.equal(raster::extent(img), ref_ext) ||
          !all.equal(dim(img), ref_dim)) {
        stop(sprintf("Images differ in resolution, extent, or dimensions and no 'res' target was provided. Mismatch found at image index %d.", i))
      }
    }
    
  } else {
    # Resample all images to the specified resolution
    e = as.vector(raster::extent(ref))
    ratio = (e[2] - e[1])/(e[4] - e[3])
    target <- raster::raster(nrow = res, ncol = floor(res * 
                                                        ratio))
    raster::crs(target) = crs(ref)
    raster::extent(target) <- e
    
    imgList <- lapply(imgList, function(img) {
      raster::resample(img, target, method = ifelse(interpolate, "bilinear", "ngb"))
    })
  }
  
  return(imgList)
}

# check imag format (RGB or Raster)
checkFormat <- function(img) {
  
  dims <- dim(img)
  vRange <- range(raster::values(img), na.rm = TRUE)
  
  if (dims[3] >= 3 && all(vRange >= 0 & vRange <= 255)) {
    format <- "RGB"
  } else if (dims[3] == 1) {
    format <- "raster"
  } else {
    stop(
      sprintf(
        "Unsupported image format: expected 1 or 3 (RGB) rasterLayer(s), got dims[3] = %d with value range [%.2f, %.2f].",
        dims[3], min(vRange), max(vRange)
      )
    )
  }
  format
}

# Extract and linearize image data
extractImageData <- function(imgList, impute.NA = TRUE) {
  format <- checkFormat(imgList[[1]])
  n <- length(imgList)
  
  # Extract pixel matrices for each channel
  extractChannel <- function(img, ch) raster::getValues(img[[ch]])
  
  if (format == "RGB") {
    R <- sapply(imgList, extractChannel, ch = 1)
    G <- sapply(imgList, extractChannel, ch = 2)
    B <- sapply(imgList, extractChannel, ch = 3)
    
    if (impute.NA) {
      # Identify background pixels (all NA across all images and channels)
      allNA <- rowSums(is.na(R), na.rm = TRUE) == n
      
      # Remove background pixels
      R <- R[!allNA, , drop = FALSE]
      G <- G[!allNA, , drop = FALSE]
      B <- B[!allNA, , drop = FALSE]
      
      # Impute missing values per pixel across images (per channel)
      impute <- function(channel) {
        apply(channel, 1, function(row) {
          if (all(is.na(row))) return(rep(NA, length(row)))
          m <- mean(row, na.rm = TRUE)
          row[is.na(row)] <- m
          return(row)
        })
      }
      R <- t(impute(R))
      G <- t(impute(G))
      B <- t(impute(B))
    } else {
      # Identify background pixels (any NA across all images and channels)
      allNA <- rowSums(is.na(R)) > 0
      
      # Remove background and NA pixels
      R <- R[!allNA, , drop = FALSE]
      G <- G[!allNA, , drop = FALSE]
      B <- B[!allNA, , drop = FALSE]
    }
    
    # Return flattened image data and background mask
    dataMatrix <- do.call(rbind, lapply(1:n, function(i) {
      c(R[, i], G[, i], B[, i])
    }))
    rownames(dataMatrix) <- names(imgList)
    attr(dataMatrix, "background") <- allNA
    return(dataMatrix)
    
  } else {
    # rasterLayer case
    rL <- sapply(imgList, extractChannel, ch = 1)
    
    allNA <- apply(is.na(rL), 1, all)
    rL <- rL[!allNA, , drop = FALSE]
    
    if (impute.NA) {
      rL <- t(apply(rL, 1, function(row) {
        if (all(is.na(row))) return(rep(NA, length(row)))
        m <- mean(row, na.rm = TRUE)
        row[is.na(row)] <- m
        return(row)
      }))
    }
    
    rownames(rL) <- names(imgList)
    attr(rL, "background") <- allNA
    return(t(rL))
  }
}

# compute PCA (or phylogenetic PCA)
computePCA <- function(mat, tree = NULL, scale = FALSE) {
  if (!is.null(tree)) {
    if (!requireNamespace("phytools")) stop("Please install the 'phytools' package.")
    match_names <- intersect(tree$tip.label, rownames(mat))
    tree <- keep.tip(tree, match_names)
    mat <- mat[tree$tip.label, , drop = FALSE]
    pca <- phytools::as.prcomp(phytools::phyl.pca(tree, mat, method = "lambda",
                                                  mode = ifelse(scale, "corr", "cov")))
  } else {
    pca <- prcomp(mat, center = TRUE, scale. = scale)
  }
  return(pca)
}

is.scaled <- function(x) {
  !is.null(attr(x, "scaled:center")) && !is.null(attr(x, "scaled:scale"))
}

# Reconstruct images based on PC scores
reconstruct_img <- function(pc.vec, imgList, dataMatrix, pca) {
  
  if (!requireNamespace("raster")) stop("Please install the 'raster' package.")
  
  # Reconstruct flattened pixel values from PC scores
  # Get PCA rotation and center
  center <- if (is.null(attr(pca, "center"))) pca$center else attr(pca, "center")
  rotation <- pca$rotation
  
  if (length(pc.vec) < ncol(rotation)) {
    pc.vec <- c(pc.vec, rep(0, ncol(rotation) - length(pc.vec)))
  }
  
  # Reconstrunct phenotype
  reconstructed <- pc.vec %*% t(rotation)
  if (!is.null(center)) {
    reconstructed <- reconstructed + matrix(center, nrow = 1)
  }
  reconstructed <- as.vector(reconstructed)  # flatten
  
  # Determine if image is RGB or signle raster layer
  format <- checkFormat(imgList[[1]])
  
  # Get background indexes
  background <- attr(dataMatrix, "background")
  
  # if data has been scaled
  scaled <- is.scaled(dataMatrix)
  
  if (scaled) {
    center_ <- attr(dataMatrix, "scaled:center")
    scale_ <- attr(dataMatrix, "scaled:scale")
    
    # Unscale
    reconstructed <- reconstructed *  scale_ +  center_
  }
  
  # Get images dimension
  img_dim <- dim(imgList[[1]])
  n_pixels <- length(background)
  
  # Clip function
  clipRGB <- function(x) pmin(pmax(x, 0), 255)
  
  if (format == "RGB") {
    n_valid <- sum(!background)
    R_full <- G_full <- B_full <- rep(NA, n_pixels)
    
    # Assume RGB channels are stacked (R1...Rn, G1...Gn, B1...Bn)
    rgb_mat <- matrix(reconstructed, nrow = 3, byrow = TRUE)
 
    R_full[!background] <- clipRGB(rgb_mat[1, ])
    G_full[!background] <- clipRGB(rgb_mat[2, ])
    B_full[!background] <- clipRGB(rgb_mat[3, ])
    
    r_template <- imgList[[1]]
    names(r_template) <- c("R", "G", "B")
    
    r <- raster::brick(
      raster::setValues(r_template[[1]], R_full),
      raster::setValues(r_template[[2]], G_full),
      raster::setValues(r_template[[3]], B_full)
    )
    return(r)
    
  } else {
    # Single-layer raster
    full <- rep(NA, n_pixels)
    full[!background] <- as.vector(reconstructed)
    
    r_template <- imgList[[1]]
    names(r_template) <- c("layer")
    
    r <- raster::setValues(r_template, full)
    return(r)
  }
}

# Wrapper function to perform image PCA
imagePCA <- function(imgList, tree = NULL, pcs = 1:2, impute.NA = TRUE, 
                     res = NULL, interpolate = FALSE, scale = FALSE, plot = FALSE, ...) {
  
  cat("Preprocessing images...\n")
  imgList <- preprocessImages(imgList, res, interpolate)
  
  cat("Extracting image data...\n")
  dataMatrix <- extractImageData(imgList, impute.NA)
  
  cat("Scaling matrix...\n")
  dataMatrix <- scale(dataMatrix)
  
  cat("Running PCA...\n")
  pca <- computePCA(dataMatrix, tree, scale)
  
  cat("Reconstructing PC origin...\n")
  pc.imgCt <- reconstruct_img(rep(0, ncol(pca$rotation)), imgList, dataMatrix, pca)
  
  cat("Reconstructing PC boundaries...\n")
  pc.minmax <-  unlist(lapply(pcs, function(i) {
    pc.vecMax <- pc.vecMin <- rep(0, ncol(pca$rotation))
    
    pc.vecMin[i] <- min(pca$x[,i])
    pc.vecMax[i] <- max(pca$x[,i])
    
    pc.imgMin <- reconstruct_img(pc.vecMin, imgList, dataMatrix, pca)
    pc.imgMax <- reconstruct_img(pc.vecMax, imgList, dataMatrix, pca)
    
    list(
      setNames(list(pc.imgMin), paste0("PC", i, "min")),
      setNames(list(pc.imgMax), paste0("PC", i, "max"))
    )
  }))
  
  out <- list(images = imgList,
              tree = if (!is.null(tree)) tree else NULL,
              data = dataMatrix,
              pca = pca,
              mask = attr(dataMatrix, "background"),
              ras = c(center = pc.imgCt, pc.minmax),
              type = checkFormat(imgList[[1]]),
              components = paste0("PC", pcs))
  
  return(out)
}

# Function to reconstruct images from PCA contributions
reconstruct_pixel_contributions <- function(pc_object, pc = 1, scale = FALSE, average = FALSE, plot = c("positive", "negative")) {
  
  ctr_img <- pc_object$ras[[1]]
  
  # get cell IDs from object
  background <- pc_object$mask
  NR <- raster::nrow(ctr_img)
  NC <- raster::ncol(ctr_img)
  
  # create raster template
  r <- raster::raster(matrix(NA, nrow = NR, ncol = NC), crs = crs(ctr_img))
  raster::extent(r) <- extent(ctr_img)
  r[!background] = 0
  
  # get pixel contributions
  rotation <- pc_object$pca$rotation[,pc]
  center <- if (is.null(attr(pc_object$pca, "center"))) pc_object$pca$center else attr(pc_object$pca, "center")
  
  # get contribution limits
  min_val <- min(rotation)
  max_val <- max(rotation)
  
  # if PCA was done on RGB images, separate each channel
  if (pc_object$type == "RGB"){
    
    n <- length(center)/3
    rotation <- lapply(1:3, function(i) {
      x <- rotation[(n * i - n + 1):(n * i)]
      x
    })
    names(rotation) <- c("R", "G", "B")
    
    # create raster stack with negative contributions to the selected PC
    rn <- lapply(1:3, function(i) {
      r[!background] = t(rotation[[i]])
      r[r > 0] = 0
      r
    })
    
    # create raster stack with positive contributions to the selected PC
    rp <- lapply(1:3, function(i) {
      r[!background] = t(rotation[[i]])
      r[r < 0] = 0
      r
    })
    
    if (scale) {
      rn <- sapply(rn, function(i) {
        i / raster::cellStats(i, "min", na.rm = TRUE)
      }, simplify = FALSE)
      rp <- sapply(rp, function(i) {
        i / raster::cellStats(i, "max", na.rm = TRUE)
      }, simplify = FALSE)
    }
    
    names(rn) <- c("R", "G", "B")
    names(rp) <- c("R", "G", "B")
    
    rn <- stack(rn)
    rp <- stack(rp)
    
    # average RGB channel contributions
    if (average) {
      rp = mean(rp, na.rm = TRUE)
      rn = mean(rn, na.rm = TRUE)
      names(rp) <- names(rn) <- c("RGB")
    }
    
  } else {
    n <- length(center)
    
    # create raster with negative contributions to the selected PC
    rn <- r
    rn[!background] = t(rotation)
    rn[rn > 0] = 0
    
    # create raster with positive contributions to the selected PC
    rp <- r
    rp[!background] = t(rotation)
    rp[rp < 0] = 0
    
  }
  
  return(list(positive = rp, negative = rn))
  
}

plot_pixel_contributions <- function(r, ...) {
  
  sign <- ifelse(cellStats(r, "min") < 0, "negative", "positive")
  
  if (sign == "positive") {
    if (nlayers(r) == 3){
      par(mfrow = c(3, 1))
      # R
      color_palette_positive <- colorRampPalette(c("grey98", "red"))(100)
      plot(r[["R"]], zlim = c(0, round(max_val, 3)), col = color_palette_positive, main = "R", ...)
      # G
      color_palette_positive <- colorRampPalette(c("grey98", "green"))(100)
      plot(r[["G"]], zlim = c(0, round(max_val, 3)), col = color_palette_positive, main = "G", ...)
      # B
      color_palette_positive <- colorRampPalette(c("grey98", "blue"))(100)
      plot(r[["B"]], zlim = c(0, round(max_val, 3)), col = color_palette_positive, main = "B", ...)
      par(mfrow = c(1,1))
    } else if (nlayers(r) == 1){
      color_palette_positive <- colorRampPalette(c("grey98", "red"))(100)
      plot(r, zlim = c(0, round(max_val, 3)), col = color_palette_positive, ...)
    }
    
  }
  if (sign == "negative") {
    if (nlayers(r) == 3){
      par(mfrow = c(3, 1))
      # R
      color_palette_negative <- colorRampPalette(c("red", "grey98"))(100)
      plot(r[["R"]], zlim = c(round(min_val, 3), 0), col = color_palette_negative, main = "R", ...)
      # G
      color_palette_negative <- colorRampPalette(c("green", "grey98"))(100)
      plot(r[["G"]], zlim = c(round(min_val, 3), 0), col = color_palette_negative, main = "G", ...)
      # B
      color_palette_negative <- colorRampPalette(c("blue", "grey98"))(100)
      plot(r[["B"]], zlim = c(round(min_val, 3), 0), col = color_palette_negative, main = "B", ...)
      par(mfrow = c(1,1))
    } else if (nlayers(r) == 1){
      color_palette_negative <- colorRampPalette(c("blue", "grey98"))(100)
      plot(r, zlim = c(round(min_val, 3), 0), col = color_palette_negative, ...)
    }
  }
} 

reverse.list <- function (ll) 
{
  nms <- unique(unlist(lapply(ll, function(X) names(X))))
  ll <- lapply(ll, function(X) setNames(X[nms], nms))
  ll <- apply(do.call(rbind, ll), 2, as.list)
  lapply(ll, function(X) X[!sapply(X, is.null)])
}

is.color <- function(x) {
  grepl("^#(?:[0-9a-fA-F]{3}){1,2}$", x) | x %in% colors()
}

plot.imgPCA <- function(imgPCA_object, pcs = NULL, group_vector = NULL, group_imgList = NULL, center = TRUE, 
                        alpha = 1, pch = 19, cex = 1, show.sd = TRUE, pixel.contributions = TRUE, 
                        colPCA = NULL, palette = NULL, display = c("points", "images", "labels")[1], 
                        img.size = 0.5, color.contributions = NULL, interpolate = FALSE, ...) {
  
  if (is.null(pcs)) {
    pcs <- imgPCA_object$components
  } else {
    if (is.numeric(pcs)) {
      if (length(pcs) == 1) {
        pcs <- c(pcs, pcs[1] + 1)
      }
      pcs <- paste0("PC", pcs)
    }
  }
  if (length(pcs) > 2) {
    pcs <- pcs[1:2]
  }
  if (any(!pcs %in% colnames(imgPCA_object$pca$x))) {
    stop(
      paste0("The specified Principal Components indices (", 
             paste(pcs, collapse = ", "), 
             ") are not found in the PCA object. "))
  }
  
  group = FALSE
  
  images <- imgPCA_object$images
  comp <- imgPCA_object$pca
  mapList <- imgPCA_object$ras
  
  summ <- summary(comp)
  pcdata <- data.frame(comp$x[,pcs])
  PCx <- pcs[1]
  PCy <- pcs[2]

  if (center){
    xlim <- max(abs(pcdata[, PCx])) * c(-1, 1)
    ylim <- max(abs(pcdata[, PCy])) * c(-1, 1)
    
  } else {
    xlim <- range(pcdata[, PCx])
    ylim <- range(pcdata[, PCy])
  }
  
  xmin <- xlim[1]
  xmax <- xlim[2]
  ymin <- ylim[1]
  ymax <- ylim[2]
  
  if (!is.null(group_vector)){
    if (length(group_vector) != nrow(pcdata)) {
      stop("Length of provided grouping vector is not equal to the number of entries in the imgPCA object.")
    } else {
      group = TRUE
      pcdata$group <- group_vector
      pcdata <- pcdata %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(across(where(is.numeric), list(mean = mean, sd = sd), .names = "{.col}_{.fn}")) %>%
        as.data.frame()
      rownames(pcdata) <- pcdata$group
      
      PCx_sd <- paste0(PCx, "_sd")
      PCy_sd <- paste0(PCy, "_sd")
      
      PCx <- paste0(PCx, "_mean")
      PCy <- paste0(PCy, "_mean")
      
    }
  }
  
  if (is.null(palette)) {
    colRamp <- (grDevices::colorRampPalette(viridis::viridis(9)))(n = 100)
  }
  else {
    colRamp <- (grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, 
                                                                     name = palette)))(n = 100)
  }
  
  if (is.null(colPCA)) {
    colPCA <- lapply(c(PCx, PCy), function(j) sapply(rownames(pcdata), 
                                                     function(i) colRamp[ceiling((pcdata[i, j] + abs(min(pcdata[, 
                                                                                                                j])) + 1) * 100/(diff(range(pcdata[, j])) + 1))]))
    colPCA <- reverse.list(colPCA)
    colPCA <- sapply(colPCA, function(i) mix.colors(unlist(i)))
  }
  colPCA <- scales::alpha(colPCA, alpha)
  
  layout_matrix <- matrix(c(
    2, 6, 3, 8, 
    1, 1, 1, 5, 
    1, 1, 1, 7, 
    1, 1, 1, 4
  ), ncol = 4, nrow = 4, byrow = TRUE)
  
  # Apply the layout
  layout(layout_matrix, widths = c(2.5, 3, 2.5, 2), heights = c(1, 2, 3, 2))
  
  par(mar = c(10,10,0,0), mgp = c(5, 1.5, 0))  
  
  plot(NULL, 
       xlim = c(xmin, xmax),
       ylim = c(ymin, ymax), 
       xlab = paste(pcs[1], "(", round(summ$importance[2, pcs[1]] * 100, 1), " %)"),
       ylab = paste(pcs[2], "(", round(summ$importance[2, pcs[2]] * 100, 1), " %)"),
       cex.axis = 2.5, cex.lab = 3
  )
  
  if (group && show.sd) {
    for (i in 1:nrow(pcdata)) {
      # Horizontal lines for PCx SD
      segments(pcdata[i,PCx] - pcdata[i,PCx_sd], pcdata[i,PCy],
               pcdata[i,PCx] + pcdata[i,PCx_sd], pcdata[i,PCy],
               col = colPCA[rownames(pcdata)])
      
      # Vertical lines for PCy SD
      segments(pcdata[i,PCx], pcdata[i,PCy] - pcdata[i,PCy_sd],
               pcdata[i,PCx], pcdata[i,PCy] + pcdata[i,PCy_sd],
               col = colPCA[rownames(pcdata)])
    }
  }
  if (display == "points") {
    points(pcdata[, PCx], pcdata[, PCy], col = colPCA[rownames(pcdata)], 
           pch = pch, cex = cex)
  }
  if (display == "images"){
    if (group) {
      if (is.null(group_imgList)) {
        group_imgList <- list()
        for (i in unique(group_vector)){
          
          idx <- which(group_vector == as.character(i))
          imgs <- images[idx]
          RR <- mean(stack(sapply(imgs, function(x) x[[1]]), na.rm =  T))
          GG <- mean(stack(sapply(imgs, function(x) x[[2]]), na.rm =  T))
          BB <- mean(stack(sapply(imgs, function(x) x[[3]]), na.rm =  T))
          group_imgList[[i]] <- stack(RR, GG, BB)
          names(group_imgList[[i]]) <- c("R", "G", "B")
        }
      }
      if (any(names(group_imgList) %in% rownames(pcdata))) {
        plotImages(pcdata[names(group_imgList), PCx], pcdata[names(group_imgList), PCy], group_imgList, interpolate = TRUE, 
                   width = img.size)
      } else {
        stop("Image names are not present in the data.")
      }
    } else {
      plotImages(pcdata[, PCx], pcdata[, PCy], images, interpolate = TRUE, 
                 width = img.size)
    }
  }
  if (display == "labels") {
    graphics::text(pcdata[, PCx], pcdata[, PCy], cex = names.cex, 
                   pos = 1, offset = 1, as.character(rownames(pcdata)))
  }
  
  # get original pca data for plotting patterns
  pcdata <- data.frame(comp$x)
  PCx <- pcs[1]
  PCy <- pcs[2]
  ext <- extent(images[[1]])
  
  # plot legend function
  plot_legend <- function(x) {
    plot(NULL, type="n", xlim=ext[1:2], ylim=ext[3:4], xlab="", ylab="", bty="n", axes = FALSE)
    rasterImage(x, ext[1], ext[3], ext[2], ext[4], interpolate = interpolate)
  }

  # Raster images with/without pixel contributions
  if (pixel.contributions){
    # Plot pixel contributions for PCx negative
    pix_contr_pcx <- reconstruct_pixel_contributions(imgPCA_object, pc = PCx, scale = TRUE, average = TRUE)
    
    rxn_cont <- pix_contr_pcx$negative 
    rxn_cont <- rxn_cont / cellStats(rxn_cont, "max", na.rm = TRUE)
    rxn_cont[is.na(rxn_cont)] = 0
    
    rxp_cont <- pix_contr_pcx$positive 
    rxp_cont <- rxp_cont / cellStats(rxp_cont, "max", na.rm = TRUE)
    rxp_cont[is.na(rxp_cont)] = 0
    
    # Pixel contributions for PCy
    pix_contr_pcy <- reconstruct_pixel_contributions(imgPCA_object, pc = PCy, scale = TRUE, average = TRUE)
    
    ryn_cont <- pix_contr_pcy$negative 
    ryn_cont <- ryn_cont / cellStats(ryn_cont, "max", na.rm = TRUE)
    ryn_cont[is.na(ryn_cont)] = 0
    
    ryp_cont <- pix_contr_pcy$positive 
    ryp_cont <- ryp_cont / cellStats(ryp_cont, "max", na.rm = TRUE)
    ryp_cont[is.na(ryp_cont)] = 0
    
  } else {
    r_temp <- mapList[[1]][[1]]
    r_temp[] = 1
    rxn_cont <- rxp_cont <- ryn_cont <- ryp_cont <- r_temp
  }
  
  if (!is.null(color.contributions)) {
    
    if (all(is.color(color.contributions))) {
      color.contributions <- rep(color.contributions, length.out = 4)
    } else {
      stop("'color.contributions' is not a R color coded vector.")
    }
    
    ## PCx
    par(mar = c(0,8,0,0))
    vals <- raster::as.matrix(rxn_cont)
    image(x = xFromCol(rxn_cont, 1:ncol(vals)),
          y = yFromRow(rxn_cont, nrow(vals):1),  # reverse rows for correct orientation
          z = t(vals[nrow(vals):1, ]),           # transpose + flip for proper display
          col = colorRampPalette(c("white", color.contributions[1]))(100),
          xlab = "", ylab = "", axes = FALSE)
    
    par(mar = c(0,8,0,0))
    vals <- raster::as.matrix(rxp_cont)
    image(x = xFromCol(rxp_cont, 1:ncol(vals)),
          y = yFromRow(rxp_cont, nrow(vals):1),  # reverse rows for correct orientation
          z = t(vals[nrow(vals):1, ]),           # transpose + flip for proper display
          col = colorRampPalette(c("white", color.contributions[2]))(100),
          xlab = "", ylab = "", axes = FALSE)
    
    ## PCy
    par(mar = c(10,0,0,5))
    vals <- raster::as.matrix(ryn_cont)
    image(x = xFromCol(ryn_cont, 1:ncol(vals)),
          y = yFromRow(ryn_cont, nrow(vals):1),  # reverse rows for correct orientation
          z = t(vals[nrow(vals):1, ]),           # transpose + flip for proper display
          col = colorRampPalette(c("white", color.contributions[3]))(100),
          xlab = "", ylab = "", axes = FALSE)
    
    par(mar = c(10,0,0,5))
    vals <- raster::as.matrix(ryp_cont)
    image(x = xFromCol(ryp_cont, 1:ncol(vals)),
          y = yFromRow(ryp_cont, nrow(vals):1),  # reverse rows for correct orientation
          z = t(vals[nrow(vals):1, ]),           # transpose + flip for proper display
          col = colorRampPalette(c("white", color.contributions[4]))(100),
          xlab = "", ylab = "", axes = FALSE)
  } else {
    ## PCx
    if (any(grepl(PCx, names(mapList)))) {
      rxn_img <- mapList[[paste0(PCx, "min")]]
      rxp_img <- mapList[[paste0(PCx, "max")]]
    } else{
      rxn_img <- images[[rownames(pcdata)[which.min(pcdata[,PCx] + abs(pcdata[,PCy]))]]]
      rxp_img <- images[[rownames(pcdata)[which.max(pcdata[,PCx] - abs(pcdata[,PCy]))]]]
    }
    
    rxn <- as.array(stack(rxn_img / 255 , rxn_cont))
    rxn[is.na(rxn)] = 1
    par(mar = c(0,8,0,0))
    plot_legend(rxn)
    
    rxp <- as.array(stack(rxp_img / 255 , rxp_cont))
    rxp[is.na(rxp)] = 1
    par(mar = c(0,8,0,0))
    plot_legend(rxp)
    
    ## PCy
    if (any(grepl(PCy, names(mapList)))) {
      ryn_img <- mapList[[paste0(PCy, "min")]]
      ryp_img <- mapList[[paste0(PCy, "max")]]
    } else{
      ryn_img <- images[[rownames(pcdata)[which.min(pcdata[,PCy] + abs(pcdata[,PCx]))]]]
      ryp_img <- images[[rownames(pcdata)[which.max(pcdata[,PCy] - abs(pcdata[,PCx]))]]]
    }
    
    ryn <- as.array(stack(ryn_img / 255 , ryn_cont))
    ryn[is.na(ryn)] = 1
    par(mar = c(10,0,0,5))
    plot_legend(ryn)
    
    
    ryp <- as.array(stack(ryp_img / 255 , ryp_cont))
    ryp[is.na(ryp)] = 1
    par(mar = c(10,0,0,5))
    plot_legend(ryp)
  }
  
  # Add horizontal bidirectional arrow (plot 6)
  par(mar = c(0,10,0,0))  
  plot.new()  # Empty plot for the arrow
  arrows(0.2, 0.25, 1, 0.25, length = 0.1, angle = 30, code = 3, col = "black", lwd = 2)
  text(0.6, 0.6, PCx, cex = 3)
  
  # Add vertical bidirectional arrow (plot 7)
  par(mar = c(10,0,0,0))  
  plot.new()  # Empty plot for the arrow
  arrows(0.25, 0.2, 0.25, 1, length = 0.1, angle = 30, code = 3, col = "black", lwd = 2)
  text(0.6, 0.6, PCy, cex = 3, srt = 270)
  
}

show_extreme_shape <- function(pc_axis = 1, scores, coords, sp_names, flip = NULL, 
                               plot = c("min", "max"), title = TRUE, image_names = NULL, ...) {
  
  plot <- match.arg(plot)
  
  # species representing min and max PC
  min_idx <- which.min(scores[, pc_axis])
  max_idx <- which.max(scores[, pc_axis])
  
  # individuals representing species 
  min_idx <- which(sp_names == names(min_idx))
  max_idx <- which(sp_names == names(max_idx))
  
  # Apply flipping if specified
  if (!is.null(flip)) {
    if ("x" %in% flip) {
      coords[, 1, ] <- -coords[, 1, ]
    }
    if ("y" %in% flip) {
      coords[, 2, ] <- -coords[, 2, ]
    }
  }
  
  # mean shape
  shape_mean <- geomorph::mshape(coords)
  
  # individual min max shapes
  shape_min <- coords[,,min_idx]
  shape_max <- coords[,,max_idx]
  
  # averaged min max shapes
  avg_min_shape <- apply(shape_min, c(1, 2), mean)
  avg_max_shape <- apply(shape_max, c(1, 2), mean)
  
  # Check which plot to display
  if ("min" %in% plot) {
    geomorph::plotRefToTarget(shape_mean, avg_min_shape, method = "TPS", ...)
    if (title) title(paste("Min PC", pc_axis, if (!is.null(image_names)) paste0("(", image_names[min_idx], ")")))
    return(Morpho::tps3d(shape_mean, avg_min_shape, shape_mean))
  }
  
  if ("max" %in% plot) {
    geomorph::plotRefToTarget(shape_mean, avg_max_shape, method = "TPS", ...)
    if (title) title(paste("Max PC", pc_axis, if (!is.null(image_names)) paste0("(", image_names[max_idx], ")")))
    return(Morpho::tps3d(shape_mean, avg_max_shape, shape_mean))
  }
}

draw_labels_with_boxes <- function(names, xlim = c(-1, 1), ylim = c(-1, 1), label_font = 2, label_size = 0.8, padding = 0.05) {
  
  # Get text width and height to position the boxes
  text_width <- strwidth(names, units = "user")  # Get the width of the text labels
  text_height <- strheight(names, units = "user")  # Get the height of the text labels
  
  # Loop through each label to draw the rectangle and then add the text
  for (i in 1:nrow(eigen)) {
    # Draw a white rectangle under the label
    rect(
      eigen[i, 1] - text_width[i] / 2 - padding,  # x position (adjust for label width)
      eigen[i, 2] - text_height[i] / 2 - padding,  # y position (adjust for label height)
      eigen[i, 1] + text_width[i] / 2 + padding,  # x position (adjust for label width)
      eigen[i, 2] + text_height[i] / 2 + padding,  # y position (adjust for label height)
      col = "white", border = NA  # White box, no border
    )
    
    # Add the text label above the rectangle
    text(
      eigen[i, 1], eigen[i, 2], labels = names[i], font = label_font, cex = label_size  # Bold text with size control
    )
  }
}

shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  # draw background text with small shift in x and y in background colour
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  # draw actual text in exact xy position in foreground colour
  text(xy$x, xy$y, labels, col=col, ... )
}

plot_morphoPCA <- function(morphoPCA, morpho_gpa, sp_names,
                           genus_cols, imglist_nobg,
                           pcx = 1, pcy = 2,
                           display = "images", cex.labels = 1.5, ...) {
  
  require(scales)  # for alpha()
  
  summ <- summary(morphoPCA)
  eigen <- morphoPCA$rotation[, c(pcx, pcy)] / max(abs(morphoPCA$rotation[, c(pcx, pcy)]))
  
  pcdata <- morphoPCA$x
  
  xlim <- range(pcdata[, pcx]) * 1.1
  ylim <- range(pcdata[, pcy]) * 1.1
  
  # scale eigen vectors to PCA lims
  x_range <- range(eigen[, 1])
  y_range <- range(eigen[, 2])
  x_scale <- diff(xlim * 0.75) / diff(x_range)
  y_scale <- diff(ylim * 0.75) / diff(y_range)
  
  eigen[, 1] <- (eigen[, 1] - x_range[1]) * x_scale + xlim[1] * 0.75
  eigen[, 2] <- (eigen[, 2] - y_range[1]) * y_scale + ylim[1] * 0.75
  
  names <- gsub("_", " ", rownames(eigen))
  names[names == "relative body length"] <- "relative\nbody length"
  names[names == "peduncle depth factor"] <- "peduncle depth\nfactor"
  
  pos <- ifelse(
    abs(eigen[, 2]) > diff(ylim) / 2,
    ifelse(eigen[, 2] > 0, 3, 1),
    ifelse(eigen[, 1] > 0, 4, 2)
  )
  
  # Layout matrix
  layout_matrix <- matrix(c(
    2, 6, 3, 8,
    1, 1, 1, 5,
    1, 1, 1, 7,
    1, 1, 1, 4
  ), ncol = 4, nrow = 4, byrow = TRUE)
  
  layout(layout_matrix, widths = c(2.5, 3, 2.5, 2), heights = c(1, 2, 3, 2))
  
  par(mar = c(10,10,0,0), mgp = c(5, 1.5, 0))
  plot(NULL,
       xlim = xlim,
       ylim = ylim,
       xlab = paste("PC", pcx, "(", round(summ$importance[2, pcx] * 100, 1), "%)"),
       ylab = paste("PC", pcy, "(", round(summ$importance[2, pcy] * 100, 1), "%)"),
       cex.axis = 2.5, cex.lab = 3
  )
  
  if (display == "points") {
    points(pcdata[, pcx], pcdata[, pcy],
           bg = alpha(genus_cols[rownames(pcdata)], 0.9),
           pch = 21, cex = 2)
  }
  
  if (display == "images") {
    plotImages(pcdata[, pcx], pcdata[, pcy],
               imglist_nobg[rownames(pcdata)],
               interpolate = TRUE, width = 0.1)
  }
  
  if (display == "labels") {
    graphics::text(pcdata[, pcx], pcdata[, pcy],
                   cex = 1.5, pos = 1, offset = 1,
                   labels = as.character(rownames(pcdata)))
  }
  
  # Arrows
  arrows(0, 0, eigen[,1], eigen[,2], length = 0.2, lwd = 4, col = "black")
  arrows(0, 0, eigen[,1], eigen[,2], length = 0.2, lwd = 2, col = "grey92")
  
  # Text with border
  shadowtext(eigen[,1], eigen[,2], labels = names, pos = pos,
             col = "#FFB499", bg = "#4B0000", cex = cex.labels, ...)
  
  # Extreme shapes
  par(mar = c(0,8,0,0))
  show_extreme_shape(pc_axis = pcx, pcdata, morpho_gpa$coords,
                     sp_names, flip = "y", plot = "min", mag = 1,
                     useRefPts = FALSE, title = FALSE)
  
  par(mar = c(0,8,0,0))
  show_extreme_shape(pc_axis = pcx, pcdata, morpho_gpa$coords,
                     sp_names, flip = "y", plot = "max", mag = 1,
                     useRefPts = FALSE, title = FALSE)
  
  par(mar = c(10,0,0,5))
  show_extreme_shape(pc_axis = pcy, pcdata, morpho_gpa$coords,
                     sp_names, flip = "y", plot = "min", mag = 1,
                     useRefPts = FALSE, title = FALSE)
  
  par(mar = c(10,0,0,5))
  show_extreme_shape(pc_axis = pcy, pcdata, morpho_gpa$coords,
                     sp_names, flip = "y", plot = "max", mag = 1,
                     useRefPts = FALSE, title = FALSE)
  
}
  
# make a color bar legend
color.bar <- function(position, cols, xlims = c(0, 1), outline = TRUE, title = NULL, digits = 1, 
                      direction = "rightwards", fsize = 1, size = c(1, 1)) {
  
  if (position == "topleft") {
    x <- par()$usr[1] + 0.05 * (par()$usr[2] - par()$usr[1])
    y <- par()$usr[4] - 0.05 * (par()$usr[4] - par()$usr[3])
  } else if (position == "bottomleft") {
    x <- par()$usr[1] + 0.05 * (par()$usr[2] - par()$usr[1])
    y <- par()$usr[3] + 0.05 * (par()$usr[4] - par()$usr[3])
  } else if (position == "topright") {
    x <- par()$usr[2] - 0.05 * (par()$usr[2] - par()$usr[1])
    y <- par()$usr[4] - 0.05 * (par()$usr[4] - par()$usr[3])
  } else if (position == "bottomright") {
    x <- par()$usr[2] - 0.05 * (par()$usr[2] - par()$usr[1])
    y <- par()$usr[3] + 0.05 * (par()$usr[4] - par()$usr[3])
  }
  
  # Adjust the color bar's length and width using size argument
  leg <- size[1]  
  lwd <- size[2]  
  
  # Define color palette and gradient
  ncols <- length(cols)
  
  # Set direction for the color bar
  if (direction %in% c("rightwards", "leftwards")) {
    X <- x + cbind(0:(length(cols) - 1) / length(cols), 1:length(cols) / length(cols)) * (leg)
    if (direction == "leftwards") {
      X <- X[nrow(X):1, ]
      if (!is.null(lims)) 
        lims <- lims[2:1]
    }
    Y <- cbind(rep(y, length(cols)), rep(y, length(cols)))
  } else if (direction %in% c("upwards", "downwards")) {
    Y <- y + cbind(0:(length(cols) - 1) / length(cols), 1:length(cols) / length(cols)) * (leg)
    if (direction == "downwards") {
      X <- X[nrow(X):1, ]
      if (!is.null(lims)) 
        lims <- lims[2:1]
    }
    X <- cbind(rep(x, length(cols)), rep(x, length(cols)))
  }
  
  # Plot the color bar outline
  if (outline) 
    lines(c(X[1, 1], X[nrow(X), 2]), c(Y[1, 1], Y[nrow(Y), 2]), lwd = lwd + 2, lend = 2)
  
  # Plot the color gradient
  for (i in 1:length(cols)) 
    lines(X[i, ], Y[i, ], col = cols[i], lwd = lwd, lend = 2)
  
  # Add axis labels for xlims
  if (direction %in% c("rightwards", "leftwards")) {
    text(x = c(x, x + leg), y =  rep(Y[nrow(Y), 1], 2), round(xlims, digits), pos = 1, cex = fsize)
  } else if (direction %in% c("upwards", "downwards")) {
    text(x = rep(Y[1, 1], 2), y = c(y, y + leg), labels = round(xlims, digits), pos = 2, cex = fsize)
  }
  
  # Title for the color bar
  if (!is.null(title)) {
    if (direction %in% c("rightwards", "leftwards")) {
      text(x = x + leg / 2, y = y, title, pos = 3, cex = fsize)
    } else if (direction %in% c("upwards", "downwards")) {
      text(x = x - 0.04 * diff(par()$usr[1:2]), y = y + leg / 2, title, pos = 3, cex = fsize, srt = 90)
    }
  }
}

# plot l1ou detected shift optimas
plot_shift_optima <- function(fit,
                              scale = FALSE,
                              normalize = FALSE,
                              color.by = c("value", "regime"), 
                              convergent.only = TRUE, 
                              edge.width = 2,
                              shift.edge.width = 3,
                              root.value = 0,
                              color.palette = colorRampPalette(c("red", "grey85", "blue")),
                              regime.palette = NULL,
                              bg_tree = "grey85",
                              show.tip.label = FALSE,
                              cex.value = 1.2,
                              title = "Shift values",
                              legend.position = "topleft",
                              cex.legend = 3, ...) {
  
  color.by <- match.arg(color.by)
  
  # Check input
  if (!inherits(fit, "l1ou")) {
    stop("fit must be a 'l1ou' object.")
  }
  
  xlim.tree <- c(0, max(ape::branching.times(fit$tree)))
  
  # Detect automatically if univariate or multivariate
  multivariate <- is.matrix(fit$shift.values) && ncol(fit$shift.values) > 1
  
  # Extract and prepare shift configurations
  sc <- fit$shift.configuration
  regimes <- names(sc)
  
  if (convergent.only) {
    is_convergent <- duplicated(regimes) | duplicated(regimes, fromLast = TRUE)
    sc <- sc[is_convergent] 
    regimes <- regimes[is_convergent]
  } else {
    is_convergent <- rep(TRUE, length(sc))
  }
  
  nEdges <- Nedge(fit$tree)
  ew <- rep(edge.width,nEdges) # to set default edge width of 1
  ew[sc] <- shift.edge.width # to widen edges with a shift
  
  if (multivariate) {
    sv <- fit$shift.values[is_convergent, , drop = FALSE]
    rownames(sv) <- sc
    sc <- matrix(rep(sc, ncol(sv)), nrow = nrow(sv), ncol = ncol(sv))
    sv <- sv[order(sc[,1]), , drop = FALSE]
    sc <- sc[order(sc[,1]), , drop = FALSE]
  } else {
    sv <- fit$shift.values[is_convergent]
    names(sv) <- sc
  }
  
  # Apply scaling or normalization if requested
  if (scale) {
    sv <- log(abs(sv)) * sv / abs(sv) # Log-transform preserving sign
    title <- sprintf("Log-scaled\n%s", title)
  } else if (normalize) {
    sv <- scale(sv, center = FALSE) / sqrt(sum(scale(sv, center = FALSE)^2))
    title <- sprintf("Normalized\n%s", title)
  } 

  # Compute optimum values
  optimum.values <- convert_shifts2regions(fit$tree, sc, sv) + root.value
  
  regime.values <- setNames(optimum.values[sc], regimes)
  
  if (color.by == "value") {
    
    # Continuous coloring by trait value
    xlims <- max(abs(optimum.values)) * c(-1, 1)
    if (all(xlims == 0)) { xlims <- c(-1, 1) }
    phytools::plotBranchbyTrait(fit$tree, optimum.values, mode = "edges",
                                edge.width = ew, xlims = xlims,
                                show.tip.label = FALSE, 
                                title = title, legend = FALSE,
                                palette = color.palette)
    # add shift values
    if (length(regime.values) > 0) {
      ape::edgelabels(text = regime.values, edge = sc,
                      bg = "white", frame = "none", adj = c(0.5, -0.5), cex = cex.value)
    }
    color.bar(legend.position, cols = color.palette(n=100), xlims = xlims, size = c(4, 4), title = title, direction = "rightwards")
    
  } else if (color.by == "regime") {
    
    unique.regimes <- sort(unique(regimes))
    n.regimes <- length(unique.regimes)
    
    # Discrete coloring by regime ID
    regime.IDs <- rep(0, nEdges)
    for (i in 1:length(regime.values)) {
      idx <- optimum.values == regime.values[i]
      regime.IDs[idx] = regimes[i]
    }

    if (is.null(regime.palette)) {
      custom_palette <- c(
        "#E69F00", # soft orange
        "#0072B2", # muted strong blue
        "#009E73", # muted green
        "#FFD700", # gold yellow
        "#87CEEB", # lighter blue
        "#32CD32", # vibrant green
        "#CC79A7", # muted pink
        "#8A2BE2"  # vibrant blue-violet
      )
      if (n.regimes - 1 <= 8) {
        regime.palette <- c(bg_tree, custom_palette[1:(n.regimes)])
      } else {
        regime.palette <- c(bg_tree, rainbow(n.regimes))
      }
    } else {
      if (length(regime.palette) < (n.regimes - 1)) {
        stop(paste0("Provided regime.palette is too short. Expected at least ", n.regimes - 1, " colors (excluding background color)."))
      }
      regime.palette <- c(bg_tree, regime.palette[1:(n.regimes - 1)])
    }
    
    # Map regimes to colors
    regime.colors <- setNames(regime.palette, c(0, unique.regimes))
    edge.colors <- regime.colors[as.character(regime.IDs)]  # match regime of each edge
    
    plot(fit$tree, show.tip.label = show.tip.label, edge.color = edge.colors, 
         edge.width = ew, main = title, cex = cex.legend)

    # add shift values
    if (length(regime.values) > 0) {
      ape::edgelabels(text = round(regime.values, 2), edge = sc,
                      bg = "white", frame = "none", adj = c(0.5, -0.5), cex = cex.value)
      legend.text <- c("Baseline", paste0("Regime ", unique.regimes))
    } else {
      legend.text <- "Baseline"
    }
    # add legend
    legend(legend.position, 
           legend = legend.text,
           col = regime.palette, 
           lwd = shift.edge.width, 
           cex = cex.legend,
           box.lwd = 0,
           bg = "white",
           title = "Shift Regimes")
  }
}

# convert l1ou detected regimes into a simmap object
l1ou2simmap <- function(fit, convergent.only = TRUE, plot = TRUE, randomize = FALSE) {
  
  nEdges <- Nedge(fit$tree)
  tree <- fit$tree
  
  # Detect if multivariate
  multivariate <- is.matrix(fit$shift.values) && ncol(fit$shift.values) > 1
  
  # Extract shift configuration and regime names
  sc <- fit$shift.configuration
  regimes <- names(sc)
  
  if (convergent.only) {
    is_convergent <- duplicated(regimes) | duplicated(regimes, fromLast = TRUE)
    sc <- sc[is_convergent] 
    regimes <- regimes[is_convergent]
  } else {
    is_convergent <- rep(TRUE, length(sc))
  }
  
  # Shift values
  if (multivariate) {
    sv <- fit$shift.values[is_convergent, , drop = FALSE]
    rownames(sv) <- sc
    sc <- matrix(rep(sc, ncol(sv)), nrow = nrow(sv), ncol = ncol(sv))
    sv <- sv[order(sc[,1]), , drop = FALSE]
    sc <- sc[order(sc[,1]), , drop = FALSE]
  } else {
    sv <- fit$shift.values[is_convergent]
    names(sv) <- sc
  }
  
  # Replace shift configuration with random nodes if randomize = TRUE
  if (randomize) {
    n_shifts <- length(sc)
    unique_regimes <- unique(regimes)
    
    # Internal node IDs: greater than Ntip(tree)
    internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode - 1)
    
    if (length(internal_nodes) < n_shifts) {
      stop("Not enough internal nodes to assign shifts.")
    }
    
    # Sample new internal nodes
    sampled_nodes <- setNames(sample(internal_nodes, n_shifts, replace = FALSE), regimes)
    
    sc <- sampled_nodes  
    names(sv) <- sc
  }

  # Convert shifts to regime regions
  optimum.values <- convert_shifts2regions(tree, sc, sv)
  regime.values <- setNames(optimum.values[sc], names(sc))  # regime to value mapping
  
  # Assign regimes to edges
  edge_regimes <- rep(0, nEdges)
  for (i in 1:length(regime.values)) {
    idx <- optimum.values == regime.values[i]
    edge_regimes[idx] <- names(regime.values)[i]
  }
  
  # Create maps list for each edge
  maps <- vector("list", nrow(tree$edge))
  names(maps) <- apply(tree$edge, 1, function(x) paste(x[1], x[2], sep = ","))
  
  for (i in 1:nrow(tree$edge)) {
    edge_length <- tree$edge.length[i]
    regime <- as.character(edge_regimes[i])
    maps[[i]] <- setNames(edge_length, regime)
  }
  
  unique.regimes <- sort(unique(edge_regimes))
  n.regimes <- length(unique.regimes)
  
  # Create mapped.edge matrix
  mapped_edge <- matrix(0, nrow = nrow(tree$edge), ncol = length(unique.regimes),
                        dimnames = list(NULL, unique.regimes))
  for (i in 1:nrow(tree$edge)) {
    regime <- as.character(edge_regimes[i])
    mapped_edge[i, regime] <- tree$edge.length[i]
  }
  
  # Create final simmap object
  mapped_tree <- tree
  mapped_tree$maps <- maps
  mapped_tree$mapped.edge <- mapped_edge
  class(mapped_tree) <- c("simmap", "phylo")
  attr(mapped_tree, "map.order") <- "right-to-left"
  mapped_tree <- phytools::reorderSimmap(mapped_tree)
  
  # Optional plotting
  if (plot) {
    plotSimmap(mapped_tree, 
               colors = setNames(c("grey90", rainbow(length(unique.regimes) - 1)), unique.regimes),   
               lwd = 2, fsize = 0.8, pts = FALSE)
  }
  
  return(mapped_tree)
}

# plot brownie.lite fit model versus model simulations
plot.model_rates <- function(model, show.single.rate = TRUE, legend.position = "top", cex = c(1, 1, 1.2, 1.8), title = "Null Distribution of Multi-rate Evolution") {
  
  simulated_mrlogLik <- sapply(model$sim, function(result) result$logL.multiple)
  observed_mrlogLik <- model$obs$logL.multiple
  observed_srlogLik <- model$obs$logL1
  p_value <- model$p.val
  
  if (length(cex) != 3) cex = rep(cex, length.out = 4)
  
  p_value_text <- paste("p-value", 
                        ifelse(p_value < 0.001, "< 0.001", 
                               paste("=", round(p_value, 3))))
  xlim <- range(c(observed_srlogLik, observed_mrlogLik, simulated_mrlogLik)) 
  expansion <- diff(xlim) * 0.1
  xlim <- c(xlim[1] - expansion, xlim[2] + expansion)
  xlab <- "Log-Likelihood"
  
  par(bg = "white", cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.5, mar = c(6,5,5,5))
  # plot simulated MI
  hist(simulated_mrlogLik, breaks = 10, main = title, xlim = xlim, 
       cex.main = cex[3], cex.axis = cex[2], cex.lab = cex[2],
       xlab = xlab,  col = "gray90",  border = "gray70",  freq = TRUE)
  
  if (show.single.rate) {
    legend.text <- c(paste("Single-rate =", round(observed_srlogLik, 3)),
                     paste("Multi-rate =", round(observed_mrlogLik, 3)))
    legend.col <- c("#FF3399", "#1F78B4")
    # add observed MI
    segments(x0 = observed_srlogLik, y0 = 0, x1 = observed_srlogLik, y1 = par("usr")[4]*0.25,  
             col = "#FF3399", lwd = 2.5)
    points(x = observed_srlogLik, y = par("usr")[4]*0.25, 
           col = "#FF3399", pch = 16, cex = cex[4])
  } else {
    legend.text <- paste("Multi-rate =", round(observed_mrlogLik, 3))
    legend.col <-  "#1F78B4"
  }
  
  # add observed MI
  segments(x0 = observed_mrlogLik, y0 = 0, x1 = observed_mrlogLik, y1 = par("usr")[4]*0.25,  
           col = "#1F78B4", lwd = 2.5)
  points(x = observed_mrlogLik, y = par("usr")[4]*0.25, 
         col = "#1F78B4", pch = 16, cex = cex[4])
  
  # add p_value
  text(observed_mrlogLik, par("usr")[4]*0.32, 
       labels = p_value_text, 
       col = "#1F78B4", cex = cex[1], font = 3)
  
  # add legend
  legend(legend.position, legend = legend.text, 
         col = legend.col, pt.cex = cex[4], pch = 16, bty = "n", cex = cex[2])
  
}

# plot dtt (and ltt)
plot.dtt <- function(dttRes, CI = 0.95, phy = NULL, cex.lab = 1.5, 
                     xlab = "Relative time from origin", ylab = "Relative disparity", ...) {
  
  ltt <- dttRes$times
  dtt.data <- dttRes$dtt
  dtt.sims <- dttRes$sim
  MDI <- paste("MDI = ", round(dttRes$MDI, 2))
  MDIpVal <- dttRes$MDIpVal
  
  median.sims <- apply(dtt.sims, 1, median)
  ylim <- round(range(c(dtt.data, dtt.sims)), 1)
  
  plot(ltt, dtt.data, xlab = xlab, ylab = ylab, cex.lab = cex.lab,
       ylim = ylim, bty = "n", type = "n", ...)
  
  if (!is.null(dtt.sims)) {
    poly <- geiger:::.dtt.polygon(dtt.sims, ltt, alpha = 1 - CI)
    polygon(poly[, "x"], poly[, "y"],
            col = geiger:::.transparency("#ADD8E680", 0.5), border = NA)
    lines(ltt, median.sims, lty = 2, col = "#195A8280", lwd = 1.5)
  }
  
  # Add LTT plot if phylogeny is provided
  if (!is.null(phy)) {
    LTT <- phytools::ltt(phy, plot = FALSE)
    max.y <- ylim[2]
    
    # Scale LTT times and ltt values
    ltt.times.rel <- LTT$times / max(LTT$times)  # relative time 0 (present) to 1 (past)
    ltt.y.scaled <- LTT$ltt / max(LTT$ltt) * max.y
    
    # Overlay LTT lines
    lines(ltt.times.rel, ltt.y.scaled, type = "s", lwd = 6, col = "grey40")
    lines(ltt.times.rel, ltt.y.scaled, type = "s", lwd = 4, col = "white")
    
    # Add secondary y-axis
    axis(4, at = seq(0, max.y, length.out = 5),
         labels = seq(0, 1, by = 0.25),
         col.axis = "black", cex.axis = 1.2)
    mtext("Proportion of lineages", side = 4, line = 3, cex = cex.lab * 0.75)
  }
  
  lines(ltt, dtt.data, type = "l", lwd = 2)
  
  if (!is.null(MDIpVal)) {
    MDIpVal <- if (MDIpVal < 0.001) {
      "p-value < 0.001"
    } else if (MDIpVal < 0.01) {
      "p-value < 0.01"
    } else {
      paste("p-value =", format(MDIpVal, digits = 3))
    }
    text(0, ylim[2], label = paste0(MDI, "; ", MDIpVal), adj = 0, font = 3, ...)
  } else {
    text(0, ylim[2], label = MDI, adj = 0, font = 3, ...)
  }
}

# compute Morphological Disparity Index p-value
compute.MDIpval <- function(dttRes, mdi.range = c(0,1), test = c("greater", "lower", "two-tailed")) {
  
  if (is.null(dttRes$sim)) {
    stop("Without simulated values is not possible to compute p-values.")
  }
  test <- match.arg(test)
  
  obs.MDI <- dttRes$MDI
  
  median.sims <- apply(dttRes$sim, 1, median)
  
  sim_MDIs <- apply(dttRes$sim, 2, function(sim) {
    geiger:::.area.between.curves(x = dttRes$times, f1 = sim, f2 = median.sims, sort(mdi.range))
  })
  
  if (test == "greater") {
    p_value <- sum(sim_MDIs >= obs.MDI) / length(sim_MDIs)
  }
  else if (test == "lower") {
    p_value <- sum(sim_MDIs <= obs.MDI) / length(sim_MDIs)
  }  else if (test == "two-tailed") {
    p_value <- sum(abs(sim_MDIs) >= abs(obs.MDI)) / length(sim_MDIs)
  }
  
  return(p_value)
}

make_Q_matrix <- function(model = c("ER", "SYM", "ARD", "DIR1", "DIR2", "DIR3"),
                          n_states, name_states = NULL, sum.zero = FALSE,
                          random = FALSE) {
  model <- match.arg(model)
  Q <- matrix(0, nrow = n_states, ncol = n_states)
  
  rate_id <- 1
  
  if (model %in% c("ER", "SYM", "ARD")) {
    for (i in 1:n_states) {
      for (j in 1:n_states) {
        if (i == j) next
        
        if (model == "ER") {
          Q[i, j] <- 1
        } else if (model == "SYM") {
          if (i < j) {
            rate <- if (random) runif(1, 0.1, 1) else rate_id
            Q[i, j] <- rate
            Q[j, i] <- rate
            rate_id <- rate_id + 1
          }
        } else if (model == "ARD") {
          Q[i, j] <- if (random) runif(1, 0.1, 1) else rate_id
          rate_id <- rate_id + 1
        }
      }
    }
  } else if (model == "DIR1") {
    for (i in 1:(n_states - 1)) {
      Q[i, i + 1] <- if (random) runif(1, 0.1, 1) else rate_id
      rate_id <- rate_id + 1
    }
  } else if (model == "DIR2") {
    for (i in 2:n_states) {
      Q[i, i - 1] <- if (random) runif(1, 0.1, 1) else rate_id
      rate_id <- rate_id + 1
    }
  } else if (model == "DIR3") {
    for (i in 1:(n_states - 1)) {
      rate <- if (random) runif(1, 0.1, 1) else rate_id
      Q[i, i + 1] <- rate
      Q[i + 1, i] <- rate
      rate_id <- rate_id + 1
    }
  }
  
  # Set diagonals so each row sums to zero
  if (sum.zero) {
    for (i in 1:n_states) {
      Q[i, i] <- -sum(Q[i, -i])
    }
  }

  if (!is.null(name_states)) {
    rownames(Q) <- name_states
    colnames(Q) <- name_states
  } else {
    rownames(Q) <- paste0("state_", 1:n_states)
    colnames(Q) <- paste0("state_", 1:n_states)
  }
  
  return(Q)
}

# make consensus simmap object from multiple Stochastic Character Reconstructions
make.consensus.simmap <- function(scm) {
  
  simmap <- scm[[1]]
  simmap$mapped.edge[] = 0
  for (i in 1:length(scm)) {
    simmap$mapped.edge <- simmap$mapped.edge + (scm[[i]]$mapped.edge * exp(scm[[i]]$logL))
  }
  simmap$mapped.edge <- simmap$mapped.edge / sum(sapply(scm, function(x) exp(x$logL)))
  
  return(simmap)
}

# get edge states from a simmap object
get_edge_states <- function(simmap) {
  states <- rep(NA, nrow(simmap$edge))
  for (i in 1:nrow(simmap$edge)) {
    segment_states <- names(simmap$maps[[i]])
    if (length(segment_states) == 1) {
      states[i] <- segment_states
    } else {
      lengths <- simmap$maps[[i]]
      sum_by_state <- tapply(lengths, names(lengths), sum)
      states[i] <- names(sum_by_state)[which.max(sum_by_state)]
    }
  }
  return(states)
}

# Compute Mutual Information contained in two state vectors
mutual_information <- function(x, y, weights = NULL) {
  table_xy <- table(x, y)
  if (is.null(weights)) {
    pxy <- table_xy / sum(table_xy)
  } else {
    # Weighted table
    df <- data.frame(x = x, y = y, w = weights)
    pxy <- with(df, tapply(w, list(x, y), sum))
    pxy[is.na(pxy)] <- 0
    pxy <- pxy / sum(pxy)
  }
  
  px <- rowSums(pxy)
  py <- colSums(pxy)
  
  mi <- 0
  for (i in 1:nrow(pxy)) {
    for (j in 1:ncol(pxy)) {
      if (pxy[i,j] > 0) {
        mi <- mi + pxy[i,j] * log2(pxy[i,j] / (px[i] * py[j]))
      }
    }
  }
  return(mi)
}

# Shannon entropy function
entropy <- function(z, weights = NULL) {
  if (is.null(weights)) {
    pz <- table(z) / length(z)
  } else {
    df <- data.frame(z = z, w = weights)
    pz <- with(df, tapply(w, z, sum))
    pz[is.na(pz)] <- 0
    pz <- pz / sum(pz)
  }
  -sum(pz * log2(pz))
}

# discretize continuous values based on number of trait levels
discretize <- function(values, states, target_proportions = NULL) {
  
  # Set equal proportions as default if no target provided
  if(is.null(target_proportions)) {
    n_states <- length(states)
    target_proportions <- rep(1/n_states, n_states)
  } else {
    # Validate custom proportions
    if(length(target_proportions) != length(states)) {
      stop("Length of target_proportions must match number of states.")
    }
    if(abs(sum(target_proportions) - 1) > .Machine$double.eps^0.5) {
      stop("Target proportions must sum to 1.")
    }
  }
  
  names(target_proportions) <- states
  
  # Calculate empirical quantiles
  breaks <- quantile(values, 
                     probs = cumsum(c(0, target_proportions)),
                     names = FALSE)
  
  # Ensure full coverage
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  
  # Discretize with proportions matching target
  vec.states <- cut(values, 
                       breaks = breaks,
                       labels = names(target_proportions),
                       include.lowest = TRUE)
  
  return(vec.states)
}

# Function to generate the conditioned rate matrix Q_B based on trait A's state
get_Q_B <- function(state_A, Q_A, n_states, corr_strength) {
  base_rate <- 0.1
  Q_B <- matrix(base_rate, n_states, n_states)
  diag(Q_B) <- -base_rate * (n_states - 1)
  
  # Enhance transitions to matching state to induce correlation
  match_state <- which(colnames(Q_A) == state_A)
  Q_B[, match_state] <- Q_B[, match_state] * (1 + 9 * corr_strength)  # up to 10x rate for corr_strength=1
  diag(Q_B) <- -rowSums(Q_B)
  
  return(Q_B)
}

# Update internal node states from simmap$node.states matrix
update_node_states <- function(simmap) {
  tree <- simmap
  n_edges <- nrow(tree$edge)
  
  parent_states <- character(n_edges)
  child_states <- character(n_edges)
  
  for (i in seq_len(n_edges)) {
    edge_map <- simmap$maps[[i]]
    states <- names(edge_map)
    
    # Parent node state: first state on edge
    parent_states[i] <- states[1]
    
    # Child node state: last state on edge
    child_states[i] <- states[length(states)]
  }
  
  node_states_mat <- cbind(parent_states, child_states)
  colnames(node_states_mat) <- c("parent_state", "child_state")
  return(node_states_mat)
}


# Update node states based on simmap$maps
update_states <- function(simmap) {
  tree <- simmap
  n_tips <- ape::Ntip(tree)
  n_nodes_internal <- ape::Nnode(tree)
  n_nodes <- n_tips + n_nodes_internal
  node_states <- character(n_nodes)
  
  # Assign parent node states from first mapped state on each edge
  for (i in seq_along(simmap$maps)) {
    edge <- tree$edge[i, ]
    parent <- edge[1]
    node_states[parent] <- names(simmap$maps[[i]])[1]
  }
  
  # Assign tip states from last mapped state on terminal edges
  for (tip in seq_len(n_tips)) {
    edge_index <- which(tree$edge[, 2] == tip)
    tip_state <- tail(names(simmap$maps[[edge_index]]), 1)
    node_states[tip] <- tip_state
  }
  
  names(node_states) <- as.character(seq_len(n_nodes))
  return(node_states)
}

# Update tip states vector from node states
update_tip_states <- function(simmap, node_states) {
  tree <- simmap
  n_tips <- ape::Ntip(tree)
  tip_states <- node_states[as.character(seq_len(n_tips))]
  names(tip_states) <- tree$tip.label
  return(tip_states)
}

# Update mapped.edge matrix summarizing time spent in each state per edge
update_mapped_edge <- function(simmap) {
  tree <- simmap
  n_edges <- nrow(tree$edge)
  
  all_states <- unique(unlist(lapply(simmap$maps, names)))
  mapped_edge <- matrix(0, nrow = n_edges, ncol = length(all_states),
                        dimnames = list(NULL, all_states))
  
  for (i in seq_len(n_edges)) {
    edge_map <- simmap$maps[[i]]
    for (state in names(edge_map)) {
      mapped_edge[i, state] <- mapped_edge[i, state] + edge_map[state]
    }
  }
  
  return(mapped_edge)
}

sample_state <- function(child_state_A, states_B, corr_strength) {
  n_states <- length(states_B)
  
  # Probability vector initialized to uniform
  prob <- rep(1 / n_states, n_states)
  
  # Find index of matching state in states_B (assuming states_B and states_A share labels or map accordingly)
  match_idx <- which(names(states_B) == child_state_A)
  
  if (length(match_idx) == 1) {
    # Increase probability of matching state proportional to corr_strength
    prob <- prob * (1 - corr_strength)
    prob[match_idx] <- prob[match_idx] + corr_strength
  }
  
  # Normalize to sum to 1 (just in case)
  prob <- prob / sum(prob)
  
  sample(states_B, 1, prob = prob)
}

simulate_correlated_trait <- function(simmap_A, Q_A, corr_strength) {
  
  states_B <- c("X", "Y")
  simmap_B <- simmap_A
  tree <- simmap_A
  edges <- tree$edge
  node_states_A <- tree$node.states
  states_A <- colnames(Q_A)
  states_B <- setNames(states_B, states_A)
  
  root_node <- (Ntip(tree) + 1)
  root_node_idx <- which(edges[,1] == (Ntip(tree) + 1))[1]
  root_state_A <- node_states_A[root_node_idx,1]
  root_state_B <- sample_state(root_state_A, states_B, corr_strength)
  maps_A <- simmap_A$maps
  node_states_B <- matrix(0, nrow = nrow(edges), ncol = ncol(edges))
  for (i in 1:nrow(edges)) {
    durations <- maps_A[[i]]
    states_orig <- names(durations)
    n_segments <- length(durations)
    
    parent_node <- edges[i,1]
    child_node <- edges[i,2]
    
    parent_state_A <- node_states_A[i,1]
    child_state_A <- node_states_A[i,2]
    
    if (i == root_node_idx){
      # if root node --> new root state
      parent_state_B <- root_state_B
    } else {
      if (parent_node == root_node) {
        parent_state_B <- node_states_B[edges[,1] == parent_node,1][1]
        
      } else {
        # keep sampled chid state from parent node
        parent_state_B <- node_states_B[edges[,2] == parent_node,2]
      }
    }
    
    # infer child state based on correlation strength
    #child_state_B <- sample(states_B, 1, prob = ifelse(states_A == child_state_A, corr_strength, 1-corr_strength))
    child_state_B <- sample_state(child_state_A, states_B, corr_strength)
    
    node_states_B[i, 1] <- parent_state_B
    node_states_B[i, 2] <- child_state_B
    
    # Build new state names vector:
    # Intermediate segments keep their original names
    if (n_segments == 1) {
      new_states <- parent_state_B  # single segment: just parent state
    } else if (n_segments == 2) {
      new_states <- c(parent_state_B, child_state_B)
    } else {
      new_states <- c(parent_state_B, states_B[states_orig[2:(n_segments-1)]], child_state_B)
    }
    
    # Assign new state names but keep durations unchanged
    simmap_B$maps[[i]] <- setNames(durations, new_states)
    
  }
  
  simmap_B$node.states <-  node_states_B
  simmap_B$states <- update_tip_states(simmap_B, update_states(simmap_B))
  simmap_B$mapped.edge <- update_mapped_edge(simmap_B)
  
  return(simmap_B)
}

diagnose_trait_correlation <- function(simmap_A, simmap_B) {
  trait_A <- as.numeric(as.factor(simmap_A$states))
  trait_B <- as.numeric(as.factor(simmap_B$states))
  
  picA <- pic(trait_A, simmap_A)
  picB <- pic(trait_B, simmap_B)
  pic_cor <- suppressWarnings(cor(picA, picB))
  
  dist_A <- dist(trait_A)
  dist_B <- dist(trait_B)
  mantel_res <- vegan::mantel(dist_A, dist_B, permutations = 999)
  
  chi_res <- chisq.test(table(trait_A, trait_B))
  
  return(list(
    pic_correlation = pic_cor,
    mantel_r = mantel_res$statistic,
    mantel_p = mantel_res$signif,
    chi_stat = chi_res$statistic,
    chi_p_value = chi_res$p.value
  ))
}

# simulate trait from a given tree and vector of discrete traits
simulate_trait <- function(tree, y, empirical.proportions = TRUE, method="BM", alpha = 1, nsim=1, ...) {
  
  y <- factor(y)
  y_levels <- levels(y)
  if (empirical.proportions) {
    y_props <- table(y) / length(y)
  } else {
    y_props <-NULL
  }
  
  if(method == "WN") {
    return(replicate(nsim, sample(y), simplify= FALSE))
  }

  if(method %in% c("BM", "OU")) {
    # Simulate node states (tips + internal)
    sim_nodes <- if(method == "BM") {
      alpha = NULL
      matrix(fastBM(tree, internal=TRUE, nsim = nsim), ncol = nsim, ...)
    } else {
      # not properly implemented
      # matrix(fastBM(tree, internal=TRUE, nsim = nsim, alpha = alpha, ...), ncol = nsim)
    }

    return(lapply(1:nsim, function(k) {
      # Convert to edge matrix (average parent-child nodes)
      edge_vals <- matrix(nrow=nrow(tree$edge), ncol=1)
      for(i in 1:nrow(tree$edge)){
        parent <- tree$edge[i,1]
        child <- tree$edge[i,2]
        edge_vals[i,] <- (sim_nodes[parent,k] + sim_nodes[child,k])/2  # Midpoint
      }
      # Discretize values following y states proportions
      discretize(edge_vals, y_levels, target_proportions = y_props)
    }))
  }
}

# Perform simulation-based test on Mutual Information between two simmap objects
MI.test <- function(simmaps1, simmaps2, nsim = 100, sim.model = "BM", num.cores = detectCores() - 1,
                    normalize = TRUE, alpha = 1, w.adjustment = c(FALSE, "root-proximal", "proximal-root"),
                    plot = TRUE, title = "Null Distribution of Mutual Information", ...) {
  
  sim.model <- match.arg(sim.model)
  w.adjustment <- match.arg(w.adjustment)
  
  # Helper to detect if object is a single map or list
  is_simmap_list <- function(x) {
    is.list(x) && inherits(x[[1]], "simmap")
  }

  # Helper to detect if simmap1 and simmap2 have the same tree
  validate_simmap_trees <- function(simmap1, simmap2) {
    
    # Basic tree structure
    if (!ape::all.equal.phylo(simmap1, simmap2)) {
      message("Tree topologies differ")
      return(FALSE)
    }

    # Edge matrices (strict comparison)
    if (!identical(simmap1$edge, simmap2$edge)) {
      message("Edge relationships differ")
      return(FALSE)
    }
    
    # Branch lengths (tolerance for rounding)
    if (!all.equal(simmap1$edge.length, simmap2$edge.length)) {
      message("Branch lengths differ")
      return(FALSE)
    }
    
    return(TRUE)
  }
  
  if (inherits(simmaps1, "simmap")) { simmaps1 <- list(simmaps1) }
  if (is_simmap_list(simmaps1)) {
    n_maps <- length(simmaps1) 
  } else {
    stop("First argument must be a simmap or a list of simmaps objects.")
  }
  
  # Prepare trait 2 input
  if (inherits(simmaps2, "simmap")) {
    if (!validate_simmap_trees(simmaps1[[1]], simmaps2)) {
      stop("First and second arguments must have the same tree.")
    }
    multiple = FALSE
  } else if (is_simmap_list(simmaps2)) {
    if (length(simmaps2) != n_maps) {
      stop("If second argument is a list, it must have the same length as the first argument.")
    }
    if (!validate_simmap_trees(simmaps1[[1]], simmaps2[[1]])) {
      stop("First and second arguments must have the same tree.")
    }
    multiple = TRUE
  } else {
    stop("Second argument must be a simmap or a list of simmaps objects.")
  }
  
  tree <- simmaps1[[1]]
  n_edges <- Nedge(tree)
  edge_lengths <- tree$edge.length
  
  weights <- edge_lengths
  if (w.adjustment == "root-proximal"){
    weights <- weights / (1 + node.depth.edgelength(tree)[tree$edge[,1]])
  } else if (w.adjustment == "proximal-root"){
    weights <- weights *  1/(1 + node.depth.edgelength(tree)[tree$edge[,1]])
  }
  
  # Parallel processing setup
  require(parallel)
  results <- mclapply(1:n_maps, function(k) {
    x <- get_edge_states(simmaps1[[k]]) 

    if (!multiple) {
      y <- get_edge_states(simmap2)
    } else {
      y <- get_edge_states(simmaps2[[k]])
    }

    # Observed MI
    obs_MI <- mutual_information(x, y, weights = weights)
    
    if (normalize) {
      entropy_x <- entropy(x, weights = weights)
      entropy_y <- entropy(y, weights = weights)
      obs_MI <- obs_MI / sqrt(entropy_x * entropy_y)
    }
    
    # Simulations
    sim_MI <- numeric(nsim)
    for (i in 1:nsim) {
      sim_y <- simulate_trait(tree, y, method = sim.model, nsim = 1, alpha = alpha)[[1]]
      sim_MI[i] <- mutual_information(x, sim_y, weights = weights)
      if (normalize) {
        entropy_x <- entropy(x, weights = weights)
        entropy_y <- entropy(sim_y, weights = weights)
        sim_MI[i] <- sim_MI[i] / sqrt(entropy_x * entropy_y)
      }
    }
    
    list(obs_MI = obs_MI, sim_MI = sim_MI)
  }, mc.cores = num.cores)
  
  # Collect results
  observed_MI <- numeric(n_maps)
  simulated_MI_all <- matrix(NA, nrow = n_maps, ncol = nsim)
  p_values <- numeric(n_maps)
  
  for (k in 1:n_maps) {
    observed_MI[k] <- results[[k]]$obs_MI
    simulated_MI_all[k, ] <- results[[k]]$sim_MI
    p_values[k] <- sum(results[[k]]$sim_MI >= observed_MI[k]) / nsim
  }
  
  overall_obs_MI <- mean(observed_MI)
  overall_p_val <- mean(p_values)
  
  res <- list(obs = observed_MI, sim = simulated_MI_all, p.val = overall_p_val, normalized = normalize)
  
  if (plot) {
    plot.MI(res, title = title)
  }
  
  return(res)
}

# Plot Mutual Information result
plot.MI <- function(MI, ci_level = 0.95, cex = 1, title = "Null Distribution of Mutual Information") {
  
  simulated_MI <- MI$sim
  observed_MI <- MI$obs
  p_value <- MI$p.val
  normalized = MI$normalized
  
  if (length(cex) == 1) cex = rep(cex, 2)
  
  if (length(observed_MI) > 1) {
    observed_MI_ci <- quantile(observed_MI, probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2))
    observed_MI <- mean(observed_MI)
    p_value <- mean(p_value)
    m = "Mean"
  } else {
    observed_MI_ci <- NULL
    m = ""
  }
  
  p_value_text <- paste(m, "p-value", 
                        ifelse(p_value < 0.001, "< 0.001", 
                               paste("=", round(p_value, 3))))
  xlim <- c(0, max(c(observed_MI, observed_MI_ci, simulated_MI)) * 1.1)
  xlab <- ifelse(normalized, "Normalized Mutual Information (NMI)", "Mutual Information (bits)")
  legend.text <- ifelse(normalized, paste(m, "Observed NMI =", round(observed_MI, 3)), paste(m, "Observed MI =", round(observed_MI, 3)))
  
  par(bg = "white", cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.5, mar = c(6,5,5,5))
  # plot simulated MI
  hist(simulated_MI, breaks = 40,  main = title, xlim = xlim, 
       xlab = xlab,  col = "gray90",  border = "gray70",  freq = TRUE)
  
  if (!is.null(observed_MI_ci)) {
    # add CI
    segments(x0 = observed_MI_ci[1], y0 = 0, x1 = observed_MI_ci[1], y1 = par("usr")[4] * 0.25, col = "gray50", lwd = 2, lty = 2)
    segments(x0 = observed_MI_ci[2], y0 = 0, x1 = observed_MI_ci[2], y1 = par("usr")[4] * 0.25, col = "gray50", lwd = 2, lty = 2)
    rect(xleft = observed_MI_ci[1], xright = observed_MI_ci[2], ybottom = 0, ytop = par("usr")[4] * 0.25,
         col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
  }
  
  # add observed MI
  segments(x0 = observed_MI, y0 = 0, x1 = observed_MI, y1 = par("usr")[4]*0.25,  
           col = "#1F78B4", lwd = 2.5)
  points(x = observed_MI, y = par("usr")[4]*0.25, 
         col = "#1F78B4", pch = 16, cex = 1.8)
  
  # add p_value
  text(observed_MI, par("usr")[4]*0.32, 
       labels = p_value_text, 
       col = "#1F78B4", cex = cex[1], font = 3)
  
  # add legend
  legend("top", legend = legend.text, 
         col = "#1F78B4", pt.cex = 1.8, pch = 16, bty = "n", cex = cex[2])
  
}

