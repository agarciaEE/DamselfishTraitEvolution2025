library(ape)
library(l1ou)

# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/")

# load custom functions
source("./Rscripts/custom_functions.R")

# load tree
tree <- read.tree("./data/damsel_subset.tre")

# post-hoc analysis: Computes bootstrap support for shift positions ############
pcs <- c("all", as.character(1:8))
niter = 100
criterion = "pBIC"
for (i in pcs[-3]) {
  
  # read model
  fit_conv <- readRDS(paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, "_new.rds"))
  # perform bootstraps
  fit_conv_bootstrap <- l1ou_bootstrap_support(fit_conv, nItrs=niter)
  # save new model
  saveRDS(fit_conv_bootstrap, paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, "_new_bootstraps.rds"))
  
  ## plotting result
  # using only 2 replicates in vastly insufficient in general,
  # but used here to make the illustrative example run faster.
  nEdges <- Nedge(tree)
  e.w <- rep(1,nEdges)
  e.w[fit_conv$shift.configuration] <- 3
  e.l <- round(fit_conv_bootstrap$detection.rate*100, digits=1)
  # to avoid annotating edges with support at or below 10%
  e.l <- ifelse(e.l>10, paste0(e.l,"%"), NA)
  pdf(paste0("./color/PC_", i, "_convergent_shifts_", criterion, "_new_bootstrapped.pdf"), width = 20, height = 20)
  plot(fit_conv, edge.label=e.l, edge.ann.cex=0.7, edge.label.adj = c(1.2, 0.5), edge.label.pos = 0, asterisk = FALSE,
       edge.label.ann=TRUE, edge.shift.ann = TRUE, edge.shift.adj = c(0.5, 1.2), cex=0.5, label.offset=0.5, edge.width=e.w)
  dev.off()
}