#=============================#
# Libraries
#=============================#
lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", 
         "lme4", "mgcv", "phytools", "corHMM",
         "lmtest", "ColorAR", "castor", "mvMORPH",
         "RPANDA","OUwie", "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

#=============================#
# Settings
#=============================#
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/")

res_folder <- "new_results"
dir.create(res_folder)
dir.create("new_figures/multivariate")

col.palette <- sample(unique(as.character(sapply(c("Zissou1", "Darjeeling1", "Darjeeling2", "FantasticFox1"), function(x) wes_palette(x, 5)))))
ncores <- 5

# load input data and custom functions
source("./Rscripts/A0.Input_data.R")

star_tree <- tree
ntips <- length(star_tree$tip.label)
star_tree$edge <- cbind(rep(ntips + 1, ntips), 1:ntips)
star_tree$edge.length <- rep(max(branching.times(tree)), ntips)
star_tree$Nnode <- 1

## Fitting the models
# Multivariate with color patterns
###################################
scale = FALSE
if(scale) { trait_data <- scale(colorPCs[tree$tip.label,]) } else { trait_data <- colorPCs[tree$tip.label,] }
cat("BM...\n")
fitMBM_cov <- mvBM(tree, trait_data, model = "BM1")
fitMBM_ind <- mvBM(tree, trait_data, model = "BM1", param = list(constraint = "diagonal"))
cat("WN...\n")
fitMWN_cov <- mvBM(star_tree, trait_data, model = "BM1")
fitMWN_ind <- mvBM(star_tree, trait_data, model = "BM1", param = list(constraint = "diagonal"))
cat("OU...\n")
fitMOU_cov <- mvOU(tree, trait_data, model = "OU1")
fitMOU_ind <- mvOU(tree, trait_data, model = "OU1", param = list(constraint = "diagonal"))

load(paste0("./rdata/ASR_habitat1.rda"))
cat("OUM habitat...\n")
scm_consensus <- make.consensus.simmap(trait_SCM_trees)
fitMOUM_ind_habitat <- mvOU(scm_consensus, trait_data, model = "OUM", param = list(constraint = "diagonal"))
fitMOUM_cov_habitat <- mvOU(scm_consensus, trait_data, model = "OUM")

load(paste0("./rdata/ASR_DietEcotype.rda"))
cat("OUM diet...\n")
scm_consensus <- make.consensus.simmap(trait_SCM_trees)
fitMOUM_ind_diet <- mvOU(scm_consensus, trait_data, model = "OUM", param = list(constraint = "diagonal"))
fitMOUM_cov_diet <- mvOU(scm_consensus, trait_data, model = "OUM")

load(paste0("./rdata/ASR_Farming.rda"))
cat("OUM symbiosis...\n")
scm_consensus <- make.consensus.simmap(trait_SCM_trees)
fitMOUM_ind_symbiosis <- mvOU(scm_consensus, trait_data, model = "OUM", param = list(constraint = "diagonal"))
fitMOUM_cov_symbiosis <- mvOU(scm_consensus, trait_data, model = "OUM")

# list models
mods_MCPat <- list(BM_cov = fitMBM_cov, WN_cov = fitMWN_cov, OU_cov = fitMOU_cov, OUM_cov_habitat = fitMOUM_cov_habitat, OUM_cov_diet = fitMOUM_cov_diet, OUM_cov_symbiosis = fitMOUM_cov_symbiosis,
                   BM_ind = fitMBM_ind, WN_ind = fitMWN_ind, OU_ind = fitMOU_ind, OUM_ind_habitat = fitMOUM_ind_habitat, OUM_ind_diet = fitMOUM_ind_diet, OUM_ind_symbiosis = fitMOUM_ind_symbiosis)

aic_MCPat <- c(BM_cov = fitMBM_cov$AIC, WN_cov = fitMWN_cov$AIC, OU_cov = fitMOU_cov$AIC, OUM_cov_habitat = fitMOUM_cov_habitat$AIC, OUM_cov_diet = fitMOUM_cov_diet$AIC, OUM_cov_symbiosis = fitMOUM_cov_symbiosis$AIC,
               BM_ind = fitMBM_ind$AIC, WN_ind = fitMWN_ind$AIC, OU_ind = fitMOU_ind$AIC,  OUM_ind_habitat = fitMOUM_ind_habitat$AIC, OUM_ind_diet = fitMOUM_ind_diet$AIC, OUM_ind_symbiosis = fitMOUM_ind_symbiosis$AIC)

saveRDS(list(AIC = aic_MCPat, best = mods_MCPat[[which.min(aic_MCPat)]]), file = "./rdata/mv_colorPCs_evo_res.rds")
rm(list = ls(pattern = "^fitM_"))

# Multivariate with morphological patterns
############################################
scale = FALSE
trait_data <- if (scale) { scale(morphoPCs[tree$tip.label,]) } else { morphoPCs[tree$tip.label,]}
cat("BM...\n")
fitMBM_cov <- mvBM(tree, trait_data, model = "BM1")
fitMBM_ind <- mvBM(tree, trait_data, model = "BM1", param = list(constraint = "diagonal"))
cat("WN...\n")
fitMWN_cov <- mvBM(star_tree, trait_data, model = "BM1")
fitMWN_ind <- mvBM(star_tree, trait_data, model = "BM1", param = list(constraint = "diagonal"))
cat("OU...\n")
fitMOU_cov <- mvOU(tree, trait_data, model = "OU1")
fitMOU_ind <- mvOU(tree, trait_data, model = "OU1", param = list(constraint = "diagonal"))

load(paste0("./rdata/ASR_habitat1.rda"))
cat("OUM habitat...\n")
scm_consensus <- make.consensus.simmap(trait_SCM_trees)
fitMOUM_ind_habitat <- mvOU(scm_consensus, trait_data, model = "OUM", param = list(constraint = "diagonal"))
fitMOUM_cov_habitat <- mvOU(scm_consensus, trait_data, model = "OUM")

load(paste0("./rdata/ASR_DietEcotype.rda"))
cat("OUM diet...\n")
scm_consensus <- make.consensus.simmap(trait_SCM_trees)
fitMOUM_ind_diet <- mvOU(scm_consensus, trait_data, model = "OUM", param = list(constraint = "diagonal"))
fitMOUM_cov_diet <- mvOU(scm_consensus, trait_data, model = "OUM")

load(paste0("./rdata/ASR_Farming.rda"))
cat("OUM symbiosis...\n")
scm_consensus <- make.consensus.simmap(trait_SCM_trees)
fitMOUM_ind_symbiosis <- mvOU(scm_consensus, trait_data, model = "OUM", param = list(constraint = "diagonal"))
fitMOUM_cov_symbiosis <- mvOU(scm_consensus, trait_data, model = "OUM")

# list models
mods_MCPat <- list(BM_cov = fitMBM_cov, WN_cov = fitMWN_cov, OU_cov = fitMOU_cov, OUM_cov_habitat = fitMOUM_cov_habitat, OUM_cov_diet = fitMOUM_cov_diet, OUM_cov_symbiosis = fitMOUM_cov_symbiosis,
                   BM_ind = fitMBM_ind, WN_ind = fitMWN_ind, OU_ind = fitMOU_ind, OUM_ind_habitat = fitMOUM_ind_habitat, OUM_ind_diet = fitMOUM_ind_diet, OUM_ind_symbiosis = fitMOUM_ind_symbiosis)

aic_MCPat <- c(BM_cov = fitMBM_cov$AIC, WN_cov = fitMWN_cov$AIC, OU_cov = fitMOU_cov$AIC, OUM_cov_habitat = fitMOUM_cov_habitat$AIC, OUM_cov_diet = fitMOUM_cov_diet$AIC, OUM_cov_symbiosis = fitMOUM_cov_symbiosis$AIC,
               BM_ind = fitMBM_ind$AIC, WN_ind = fitMWN_ind$AIC, OU_ind = fitMOU_ind$AIC,  OUM_ind_habitat = fitMOUM_ind_habitat$AIC, OUM_ind_diet = fitMOUM_ind_diet$AIC, OUM_ind_symbiosis = fitMOUM_ind_symbiosis$AIC)

saveRDS(list(AIC = aic_MCPat, best = mods_MCPat[[which.min(aic_MCPat)]]), file = "./rdata/mv_morphoPCs_evo_res.rds")
rm(list = ls(pattern = "^fitM_"))