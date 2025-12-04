args <- commandArgs(trailingOnly = TRUE)

library(ape)
library(l1ou)

TREE     <- args[1]
RDS_FILE <- args[2]
NITER    <- ifelse(length(args) < 3, 100, as.integer(args[3]))

# set working directory
setwd(DIR)

OUT_RDS <- sub(".rds", "bootstraps.rds", RDS_FILE)
OUT_PLOT <- sub(".rds", "bootstraps.pdf", RDS_FILE)

# load tree
tree <- read.tree(TREE)

# post-hoc analysis: Computes bootstrap support for shift positions ############

# read model
fit_conv <- readRDS(RDS_FILE)

# perform bootstraps
fit_conv_bootstrap <- l1ou_bootstrap_support(fit_conv, nItrs=NITER)

# save new model
saveRDS(fit_conv_bootstrap, OUT_RDS)

## plotting result
nEdges <- Nedge(tree)
e.w <- rep(1,nEdges)
e.w[fit_conv$shift.configuration] <- 3
e.l <- round(fit_conv_bootstrap$detection.rate*100, digits=1)
# to avoid annotating edges with support at or below 10%
e.l <- ifelse(e.l>10, paste0(e.l,"%"), NA)
pdf(OUT_PLOT, width = 20, height = 20)
plot(fit_conv, edge.label=e.l, edge.ann.cex=0.7, edge.label.adj = c(1.2, 0.5), edge.label.pos = 0, asterisk = FALSE,
     edge.label.ann=TRUE, edge.shift.ann = TRUE, edge.shift.adj = c(0.5, 1.2), cex=0.5, label.offset=0.5, edge.width=e.w)
dev.off()
