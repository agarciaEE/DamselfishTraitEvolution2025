lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", 
         "lme4", "mgcv", "phytools", "corHMM",
         "lmtest", "ColorAR", "castor", "l1ou",
         "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

# load custom functions
source("./scripts/A1.Input_data.R")

##=============================================##
##                                             ##
##        Plot l1ou bootstraps results         ##
##                                             ##
##=============================================##
dir.create("./figures/l1ou", showWarnings = FALSE)

##### Coloration patterns (PCA data)
#########################################

# Custom plots:
###############
pcs <- 1:ncol(colorPCs) # first plot all, then select those with convergent regimes
pcs <- c(2,4,6,7)
convergent_only = FALSE
criterion = "pBIC" 
bs.min = NULL
bs.thresholds   <- c(0.5, 0.7, 0.9, 0.95)
if (!is.null(bs.min)) bs.thresholds <- bs.thresholds[bs.thresholds > bs.min]
bs.col.palette  <- RColorBrewer::brewer.pal(length(bs.thresholds) + 1, "Blues")
ncol = 4

grDevices::cairo_pdf(paste0("./figures/l1ou/l1ou_color_", ifelse(convergent_only, "conv", "all"), "_pBic_PC1-8_bs_all.pdf"), 
                     height = 11 * length(pcs)/ncol + 4, width = 6 * ncol)

layout(rbind(matrix(1:(length(pcs) * 2), ncol = ncol), rep(length(pcs) * 2 + 1,ncol)),
       widths = rep(1, ncol),
       heights = c(rep(c(0.1,1), length(pcs)/ncol), 0.1))
# par(mfrow = c(length(pcs)+1,1))
# # Plot the tree first
# par(mar = c(2, 2, 3, 2))  
# plot.phylo(damsel_tree_subset, 
#            type = "tidy",
#            direction = "rightwards",  
#            show.tip.label = TRUE,
#            cex = 0.5,                 
#            label.offset = 2)
# Loop through the PCs
col_l1ou_res <- data.frame()
for (i in pcs) {
  
  ind_file <- paste0("./results/l1ou/l1ou_models/colPC_", i, "_damsel_l1ou_estimated_ind_shifts_", criterion, ".rds")
  conv_file <- paste0("./results/l1ou/l1ou_models/colPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds")
  if (file.exists(conv_file)) {
    fit <- readRDS(conv_file)
    conv <- TRUE
  } else {
    cat("No convergent shifts for PC", i, ".\n")
    fit <- readRDS(ind_file)
    conv <- FALSE
  }
  
  xmax <- ceiling(max(ape::branching.times(fit$tree)) / 10) * 10
  xlim.tree <- c(-5, xmax)
  nEdges <- Nedge(fit$tree)
  
  # read bootstrap file
  bs_file <- paste0("./results/l1ou/bootstraps/colPC_", i, "_damsel_l1ou_estimated_", ifelse(conv, "conv", "ind"), "_shifts_", criterion, "bootstraps.rds")
  if (file.exists(bs_file)) {
    fit_bootstrap <- readRDS(bs_file)
    bootstrap_values <- setNames(fit_bootstrap$detection.rate,
                                 1:nEdges)
  } else {
    bootstrap_values <- setNames(rep(0,nEdges),
                                 1:nEdges)
  }

  l1ou_res <- l1ou_extract_shift_tips(fit) %>%
    mutate(bootstrap_support = bootstrap_values[shift_edge])
  if (nrow(l1ou_res) > 0) {
    col_l1ou_res <- rbind(col_l1ou_res, 
                          data.frame(PC = paste0("colorPC", i),
                                     l1ou_res))
  }
  
  img <- png::readPNG(paste0("./figures/pca_legends/colorPC", i, "_legend.png"))  
  cropped <- img[round(dim(img)[1]*0.2):round(dim(img)[1]*0.8), , ]
  par(mar = c(0, 1, 2, 2))  
  plot(0:1, 0:1, type = "n", xlab = "", ylab = "", 
       main = paste("color PC", i), cex.main = 3, axes = FALSE, asp = FALSE) 
  rasterImage(cropped, 0, 0, 1, 1)
  
  par(mar = c(0, 1, 2, 0))  
  
  plot_shift_optima_bs(fit, 
                    color.by = "regime", 
                    convergent.only = convergent_only,
                    shift.edge.width = 5, 
                    shift.label.adj = c(0.5, -0.5),
                    cex.value = 1.5, 
                    show.tip.label = FALSE, 
                    show.bootstrap.support = "circle",
                    bs.values = bootstrap_values,
                    bs.min = bs.min,
                    bs.thresholds = bs.thresholds,
                    bs.col.palette = bs.col.palette,
                    bs.label.cex= 3,
                    bs.label.adj = c(-3, 0),
                    bs.label.pos = 0,
                    bs.legend.position = NULL,
                    bs.pie.border = colorspace::darken(bs.col.palette, 0.25),
                    bs.pie.bg = "white",
                    title = "", 
                    cex.main = 3,
                    legend.position = "topleft",
                    cex.legend = 2,
                    x.lim = xlim.tree)
}

# add legend 
par(mar = c(0, 0, 0, 0))
plot.new()
# safer ASCII labels (no Unicode issues on some PDF devices)
thr_labels <- c(
  ifelse(is.null(bs.min),
         paste0("<",  bs.thresholds[1]),
         paste0("[",  bs.min, ", ",
                bs.thresholds[1], ")")),
  paste0("[",  bs.thresholds[-length(bs.thresholds)], ", ",
         bs.thresholds[-1], ")"),
  paste0(">=", bs.thresholds[length(bs.thresholds)])
)

legend("top",
       legend   = thr_labels,
       horiz    = TRUE,
       pt.bg    = bs.col.palette,
       pch      = 21,
       pt.cex   = 5,
       col      = colorspace::darken(bs.col.palette, 0.25),
       x.intersp = 1.2,
       y.intersp = 1,
       cex      = 2,
       bty      = "n",
       title    = "Bootstrap support")

dev.off()

# Save data set with shifts and bootstraps detected
write.csv(col_l1ou_res, "./results/l1ou/color_shifts_bs.csv")

##### Morphological patterns (PCA data)
#########################################
pcs <- 1:8 # first plot all, then select those with convergent regimes
pcs <- c(1,3,4,6)
convergent_only = FALSE
criterion = "pBIC"
bs.min = NULL
bs.thresholds   <- c(0.5, 0.7, 0.9, 0.95)
if (!is.null(bs.min)) bs.thresholds <- bs.thresholds[bs.thresholds > bs.min]
bs.col.palette  <- RColorBrewer::brewer.pal(length(bs.thresholds) + 1, "Blues")
ncol = 4

# Custom plots:
###############
grDevices::cairo_pdf(paste0("./figures/l1ou/l1ou_morpho_", ifelse(convergent_only, "conv", "all"), "_pBic_PC1-8_bs_all.pdf"), 
                     height = 11 * length(pcs)/ncol + 4, width = 6 * ncol)

layout(rbind(matrix(1:(length(pcs) * 2), ncol = ncol), rep(length(pcs) * 2 + 1,ncol)),
       widths = rep(1, ncol),
       heights = c(rep(c(0.1,1), length(pcs)/ncol), 0.1))
# par(mfrow = c(length(pcs)+1,1))
# # Plot the tree first
# par(mar = c(2, 2, 3, 2))  
# plot.phylo(damsel_tree_subset, 
#            type = "tidy",
#            direction = "rightwards",  
#            show.tip.label = TRUE,
#            cex = 0.5,                 
#            label.offset = 2)
# Loop through the PCs
morpho_l1ou_res <- data.frame()
for (i in pcs) {
  
  ind_file <- paste0("./results/l1ou/l1ou_models/morphoPC_", i, "_damsel_l1ou_estimated_ind_shifts_", criterion, ".rds")
  conv_file <- paste0("./results/l1ou/l1ou_models/morphoPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds")
  if (file.exists(conv_file)) {
    fit <- readRDS(conv_file)
    conv <- TRUE
  } else {
    cat("No convergent shifts for PC", i, ".\n")
    fit <- readRDS(ind_file)
    conv <- FALSE
  }
  
  xmax <- ceiling(max(ape::branching.times(fit$tree)) / 10) * 10
  xlim.tree <- c(-5, xmax)
  nEdges <- Nedge(fit$tree)
  
  # read bootstrap file
  bs_file <- paste0("./results/l1ou/bootstraps/morphoPC_", i, "_damsel_l1ou_estimated_", ifelse(conv, "conv", "ind"), "_shifts_", criterion, "bootstraps.rds")
  if (file.exists(bs_file)) {
    fit_bootstrap <- readRDS(bs_file)
    bootstrap_values <- setNames(fit_bootstrap$detection.rate,
                                 1:nEdges)
  } else {
    bootstrap_values <- setNames(rep(0,nEdges),
                                 1:nEdges)
  }

  l1ou_res <- l1ou_extract_shift_tips(fit) %>%
    mutate(bootstrap_support = bootstrap_values[shift_edge])
  if (nrow(l1ou_res) > 0) {
    morpho_l1ou_res <- rbind(morpho_l1ou_res, 
                          data.frame(PC = paste0("morphoPC", i),
                                     l1ou_res))
  }
  
  img <- png::readPNG(paste0("./figures/pca_legends/morphoPC", i, "_legend.png"))  
  cropped <- img[round(dim(img)[1]*0.1):round(dim(img)[1]*0.9), , ]
  par(mar = c(0, 1, 2, 2))  
  plot(0:1, 0:1, type = "n", xlab = "", ylab = "", 
       main = paste("morpho PC", i), cex.main = 3, axes = FALSE, asp = FALSE)
  rasterImage(cropped, 0, 0, 1, 1)
 
  par(mar = c(0, 1, 2, 0))  
  
  plot_shift_optima_bs(fit, 
                       color.by = "regime", 
                       convergent.only = convergent_only,
                       shift.edge.width = 5, 
                       shift.label.adj = c(0.5, -0.5),
                       cex.value = 1.5, 
                       show.tip.label = FALSE, 
                       show.bootstrap.support = "circle",
                       bs.values = bootstrap_values,
                       bs.min = bs.min,
                       bs.thresholds = bs.thresholds,
                       bs.col.palette = bs.col.palette,
                       bs.label.cex= 3,
                       bs.label.adj = c(-3, 0),
                       bs.label.pos = 0,
                       bs.legend.position = NULL,
                       bs.pie.border = colorspace::darken(bs.col.palette, 0.25),
                       bs.pie.bg = "white",
                       title = "", 
                       cex.main = 3,
                       legend.position = "topleft",
                       cex.legend = 2,
                       x.lim = xlim.tree)
  

}

# add legend 
par(mar = c(2, 0, 0, 0))
plot.new()
# safer ASCII labels (no Unicode issues on some PDF devices)
thr_labels <- c(
  ifelse(is.null(bs.min),
         paste0("<",  bs.thresholds[1]),
         paste0("[",  bs.min, ", ",
                bs.thresholds[1], ")")),
  paste0("[",  bs.thresholds[-length(bs.thresholds)], ", ",
         bs.thresholds[-1], ")"),
  paste0(">=", bs.thresholds[length(bs.thresholds)])
)
legend("top",
       legend   = thr_labels,
       horiz    = TRUE,
       pt.bg    = bs.col.palette,
       pch      = 21,
       pt.cex   = 5,
       col      = colorspace::darken(bs.col.palette, 0.25),
       x.intersp = 1.2,
       y.intersp = 1,
       cex      = 2,
       bty      = "n",
       title    = "Bootstrap support")

dev.off()

# Save data set with shifts and bootstraps detected
write.csv(morpho_l1ou_res, "./results/l1ou/morpho_shifts_bs.csv")
