lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", 
         "lme4", "mgcv", "phytools", "corHMM",
         "lmtest", "ColorAR", "castor",
         "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

dir.create("figures/l1ou", recursive = TRUE, showWarnings = FALSE)
dir.create("results/l1ou", recursive = TRUE, showWarnings = FALSE)

# load input data and custom functions
source("./scripts/A1.Input_data.R")

##=============================================##
##                                             ##
##                Run analysis                 ##
##                                             ##
##=============================================##

##### Coloration patterns (PCA data)
#########################################

pcs <- c("all", as.character(1:8))
criterion = "pBIC" # try BIC too
ncores = 5
max.nShifts = 50
for (i in pcs){
  cat("Carrying out l1ou on PC", i, "...")
  Y <- get(paste0("colorPC", ifelse(i=="all", "s", i)))
  # Adjusts the tree and traits to meet the requirements of estimate_shift_configuration
  #################################################
  trait.data <- adjust_data(tree, Y, normalize = FALSE)
  
  # Detects convergent regimes under an OU model
  #################################################
  ## first fit a model to find individual shifts (no convergence assumed):
  fit_ind <- estimate_shift_configuration(trait.data$tree, trait.data$Y, 
                                          criterion=criterion, nCores = ncores, 
                                          max.nShifts = max.nShifts, quietly = FALSE)

  # Computes the information criterion score for a given configuration
  #################################################
  print(configuration_ic(trait.data$tree, fit_ind$Y, fit_ind$shift.configuration, criterion=criterion))
  
  pdf(paste0("./figures/l1ou/colorPC_", i, "_independent_shifts_", criterion, ".pdf"), width = 20, height = 25)
  ew <- rep(1,nEdges) # to set default edge width of 1
  ew[fit_ind$shift.configuration] <- 3 # to widen edges with a shift
  plot(fit_ind, cex=0.5, label.offset=0.02, edge.width=ew)
  dev.off()
  
  # save model
  saveRDS(fit_ind, paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_ind_shifts_", criterion, ".rds"))
  
  cat(fit_ind$nShifts, "shifts detected\n")
  
  if (fit_ind$nShifts > 1) {
    cat("Estimating convergent regimes... ")
    ## then detect which of these shifts are convergent:
    fit_conv <- estimate_convergent_regimes(fit_ind, criterion="pBIC", nCores = 5)
    print(configuration_ic(trait.data$tree, fit_conv$Y, fit_conv$shift.configuration, criterion=criterion))
  
    pdf(paste0("./figures/l1ou/colorPC_", i, "_convergent_shifts_", criterion, ".pdf"), width = 20, height = 25)
    ew <- rep(1,nEdges) # to set default edge width of 1
    ew[fit_conv$shift.configuration] <- 3 # to widen edges with a shift
    plot(fit_conv, cex=0.5, label.offset=0.02, edge.width=ew)
    dev.off()
    
    # save model
    saveRDS(fit_conv, paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds"))
  }
}

# Custom plots:
###############
pcs <- as.character(1:8) # only univariate analysis
convergent_only = FALSE
criterion = "pBIC" 
pdf(paste0("./figures/l1ou/l1ou_color_", ifelse(convergent_only, "conv", "all"), "_pBic_PC1-8.pdf"), 
    width = 4 * length(pcs), height = 10)

layout(matrix(1:(length(pcs)), nrow  = 1), widths = c(rep(1, length(pcs))))
# # Plot the tree first
# par(mar = c(2, 2, 3, 2))  
# plot.phylo(damsel_tree_subset, type = "tidy", direction = "rightwards",  
#            show.tip.label = TRUE, cex = 0.5, label.offset = 2)
# Loop through the PCs
col_l1ou_res <- data.frame()
for (i in pcs) {
  ind_file <- paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_ind_shifts_", criterion, ".rds")
  conv_file <- paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds")
  if (file.exists(conv_file)) {
    fit <- readRDS(conv_file)
  } else {
    cat("No convergent shifts for PC", i, ".\n")
    fit <- readRDS(ind_file)
  }
  
  if (i == pcs[1]) { 
    par(mar = c(2, 1, 4, 2))  # smaller left margin for first PC plot
  } else {
    par(mar = c(2, 4, 4, 2))  # regular margins for the rest
  }

  l1ou_res <- l1ou_extract_shift_tips(fit)
  if (nrow(l1ou_res) > 0) {
    col_l1ou_res <- rbind(col_l1ou_res, 
                      data.frame(PC = paste0("colorPC", i),
                                 l1ou_res))
  }
  
  plot_shift_optima(fit, 
                    color.by = "regime", 
                    scale=FALSE,
                    normalize=FALSE,
                    convergent.only = convergent_only,
                    shift.edge.width = 5, 
                    cex.value = 1.1, 
                    cex = 0.5, 
                    show.tip.label = FALSE, 
                    title = paste0("color PC", i), 
                    legend.position = "topleft",
                    cex.main = 3,
                    cex.legend = 1.2)
}
dev.off()

as_tibble(col_l1ou_res) %>% print(n = Inf)
col_l1ou_res %>%
  filter(PC == "colorPC7")
##### Morphological traits (PCA data)
#########################################
pcs <- c("all", as.character(1:8))
criterion = "pBIC" # try BIC too
ncores = 5
max.nShifts = 50
for (i in pcs){
  
  Y <- get(paste0("morphoPC", ifelse(i=="all", "s", i)))
  # Adjusts the tree and traits to meet the requirements of estimate_shift_configuration
  #################################################
  trait.data <- adjust_data(tree, Y, normalize = FALSE)
  
  # Detects convergent regimes under an OU model
  #################################################
  ## first fit a model to find individual shifts (no convergence assumed):
  fit_ind <- estimate_shift_configuration(trait.data$tree, trait.data$Y, criterion=criterion, nCores = ncores, max.nShifts = max.nShifts, quietly = FALSE)
  
  # Computes the information criterion score for a given configuration
  #################################################
  print(configuration_ic(trait.data$tree, fit_ind$Y, fit_ind$shift.configuration, criterion=criterion))
  
  pdf(paste0("./figures/l1ou/morphoPC_", i, "_independent_shifts_", criterion, ".pdf"), width = 20, height = 25)
  ew <- rep(1,nEdges) # to set default edge width of 1
  ew[fit_ind$shift.configuration] <- 3 # to widen edges with a shift
  plot(fit_ind, cex=0.5, label.offset=0.02, edge.width=ew)
  dev.off()
  
  # save model
  saveRDS(fit_ind, paste0("./rdata/morphoPC_", i, "_damsel_l1ou_estimated_ind_shifts_", criterion, ".rds"))
  
  cat(fit_ind$nShifts, "shifts detected\n")
  
  if (fit_ind$nShifts > 1) {
    ## then detect which of these shifts are convergent:
    fit_conv <- estimate_convergent_regimes(fit_ind, criterion="pBIC", nCores = 5)
    print(configuration_ic(trait.data$tree, fit_conv$Y, fit_conv$shift.configuration, criterion=criterion))
    
    pdf(paste0("./figures/l1ou/morphoPC_", i, "_convergent_shifts_", criterion, ".pdf"), width = 20, height = 25)
    ew <- rep(1,nEdges) # to set default edge width of 1
    ew[fit_conv$shift.configuration] <- 3 # to widen edges with a shift
    plot(fit_conv, cex=0.5, label.offset=0.02, edge.width=ew)
    dev.off()
    
    # save model
    saveRDS(fit_conv, paste0("./rdata/morphoPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds"))
  }
}

# Custom plots:
###############
convergent_only = FALSE
criterion = "pBIC" 

pdf(paste0("./figures/l1ou/l1ou_morpho_", ifelse(convergent_only, "conv", "all"), "_pBic_PC1-8.pdf"), width = 4 * length(pcs), height = 10)

layout(matrix(1:(length(pcs)), nrow  = 1), widths = c(rep(1, length(pcs))))
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
  ind_file <- paste0("./rdata/morphoPC_", i, "_damsel_l1ou_estimated_ind_shifts_", criterion, ".rds")
  conv_file <- paste0("./rdata/morphoPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds")
  if (file.exists(conv_file)) {
    fit <- readRDS(conv_file)
  } else {
    cat("No convergent shifts for PC", i, ".\n")
    fit <- readRDS(ind_file)
  }
  
  if (i == pcs[1]) { 
    par(mar = c(2, 1, 4, 2))  # smaller left margin for first PC plot
  } else {
    par(mar = c(2, 4, 4, 2))  # regular margins for the rest
  }
  
  l1ou_res <- l1ou_extract_shift_tips(fit)
  if (nrow(l1ou_res) > 0) {
    morpho_l1ou_res <- rbind(morpho_l1ou_res, 
                          data.frame(PC = paste0("morphoPC", i),
                                     l1ou_res))
  }
  
  plot_shift_optima(fit, 
                    color.by = "regime", 
                    scale=FALSE,
                    normalize=FALSE,
                    convergent.only = convergent_only,
                    shift.edge.width = 5, 
                    cex.value = 1.1, 
                    cex = 0.5, 
                    show.tip.label = FALSE, 
                    title = paste0("morpho PC", i), 
                    legend.position = "topleft",
                    cex.main = 3,
                    cex.legend = 1.2)
}

dev.off()

