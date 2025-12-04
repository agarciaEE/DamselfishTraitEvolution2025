library(ape)
library(l1ou)

# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

# load custom functions
# source("./Rscripts/custom_functions.R")

# load tree
tree <- read.tree("./data/damsel_subset.tre")

# post-hoc analysis: Computes bootstrap support for shift positions ############
pcs <- as.character(1:8)
niter = 100
criterion = "pBIC"
for (i in pcs) {
  
  l1ou_file_ind <- paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_ind_shifts_", criterion, ".rds")
  l1ou_file_conv <- paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds")
  
  mode = NULL
  if (file.exists(l1ou_file_conv)) {
    l1ou_file <- l1ou_file_conv
    mode = "conv"
  } else if (file.exists(l1ou_file_ind)) {
    l1ou_file <- l1ou_file_ind
    mode = "ind"
  } else {
    message("Ind and conv file not found for PC", i)
  }
  # read model
  fit <- readRDS(l1ou_file)
  
  # perform bootstraps
  fit_conv_bootstrap <- l1ou_bootstrap_support(fit, nItrs=niter)
  # save new model
  saveRDS(fit_conv_bootstrap, paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_", mode, "_shifts_", criterion, "_bootstraps.rds"))
  
  ## plotting result
  # using only 2 replicates in vastly insufficient in general,
  # but used here to make the illustrative example run faster.
  nEdges <- Nedge(tree)
  e.w <- rep(1,nEdges)
  e.w[fit_conv$shift.configuration] <- 3
  e.l <- round(fit_conv_bootstrap$detection.rate*100, digits=1)
  # to avoid annotating edges with support at or below 10%
  e.l <- ifelse(e.l>10, paste0(e.l,"%"), NA)
  pdf(paste0("./color/PC_", i, "_", mode, "_shifts_", criterion, "_new_bootstrapped.pdf"), width = 20, height = 20)
  plot(fit_conv, edge.label=e.l, edge.ann.cex=0.7, edge.label.adj = c(1.2, 0.5), edge.label.pos = 0, asterisk = FALSE,
       edge.label.ann=TRUE, edge.shift.ann = TRUE, edge.shift.adj = c(0.5, 1.2), cex=0.5, label.offset=0.5, edge.width=e.w)
  dev.off()
}
