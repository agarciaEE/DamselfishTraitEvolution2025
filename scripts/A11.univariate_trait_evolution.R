lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", "phyloTop",
         "lme4", "mgcv", "phytools", "corHMM", "ggpmisc",
         "lmtest", "ColorAR", "castor", "tibble", "phyloTop",
         "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

#=============================#
#   Settings
#=============================#
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

res_folder <- "results"
dir.create(res_folder)
dir.create("figures/univariate/", recursive = TRUE)
dir.create("results/univariate/", recursive = TRUE)

col.palette <- sample(unique(as.character(sapply(c("Zissou1", "Darjeeling1", "Darjeeling2", "FantasticFox1"), function(x) wes_palette(x, 5)))))
ncores <- 5
geo_model = FALSE

# load input data and custom functions
source("./scripts/A1.Input_data.R")

pomacentridae_regs <- BioGeoBEARS::getranges_from_LagrangePHYLIP("./data/CB2013_pomacentridae_geodata_input.txt")@df
region_names <- setNames(area_names <- c("IO", "IAA", "CPO", "EPO", "AO", "TS"), LETTERS[1:6])

#=============================#
# Trait Evolution Models
#=============================#

# color
trait_evo_res <- list()
trait_evo_fits <- list()

for (v in 1:ncol(colorPCs)) {
  trait_vector <- setNames(scale(colorPCs[tree$tip.label, v])[,1], tree$tip.label)
  
  # Fit BM, WN, OU
  cat("BM...\n")
  fit_BM <- fitContinuous(tree, trait_vector, model = "BM")
  cat("WN...\n")
  fit_WN <- fitContinuous(tree, trait_vector, model = "white")
  cat("OU...\n")
  fit_OU <- fitContinuous(tree, trait_vector, model = "OU")

  fit_OUM_list <- list()
  for (trait in discrete_traits) {
    asr_path <- paste0("./rdata/ASR_", trait, ".rda")
    if (file.exists(asr_path)) {
      load(asr_path)  # Should load `fit.marginal`
      
      # Format data for OUwie (assumes fit.marginal is loaded properly)
      ouwie_data <- data.frame(fit.marginal$data.legend, 
                               X = trait_vector[fit.marginal$data.legend$sp])
      ouwie_data <- ouwie_data[fit.marginal$phy$tip.label,]
      
      cat("OUM:", trait, "...\n")
      fit_OUM_list[[paste0("OUM_", trait)]] <- OUwie(fit.marginal$phy, ouwie_data, model = "OUM", algorithm = "invert", get.root.theta = FALSE)
      
    } else {
      warning("Trait rda file not found!", immediate. = TRUE)
    }
  }
  # Collect AICs and model fits
  aic_vec <- c(BM = fit_BM$opt$aic, WN = fit_WN$opt$aic, OU = fit_OU$opt$aic, sapply(fit_OUM_list, function(m) m$AIC))
  models <- c(list(BM = fit_BM, WN = fit_WN, OU = fit_OU), fit_OUM_list)
  
  trait_evo_res[[paste0("colorPC", v)]] <- list(AIC = aic_vec, best_model = models[[which.min(aic_vec)]])
  trait_evo_fits[[paste0("colorPC", v)]] <- models
}

# Save results
saveRDS(trait_evo_res, "./rdata/uni_colorPC_res.rds")
saveRDS(trait_evo_fits, "./rdata/uni_colorPC_evo_models.rds")

DD_evo_results <- list()
DD_evo_fits <- list()

for (v in 1:ncol(colorPCs)) {
  cat(sprintf("Working on color PC%s...\n", v))
  
  trait_vector <- setNames(scale(colorPCs[tree$tip.label, v])[,1], tree$tip.label)
  
  # Fit models with error handling
  cat("Running Matching Competition model....\n")
  fit_MC <- tryCatch({ fit_t_comp(tree, trait_vector, model = "MC") }, error = function(e) NULL)
  cat("Running Linear Diversity Dependent model....\n")
  fit_DDlin <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDlin") }, error = function(e) NULL)
  cat("Running Exponential Diversity Dependent model....\n")
  fit_DDexp <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDexp") }, error = function(e) NULL)
  
  if (geo_model) {
    # Fit models with error handling
    cat("Running Matching Competition geo model....\n")
    fit_MC_geo <- tryCatch({ fit_t_comp(tree, trait_vector, model = "MC", geography.object=geo_obj) }, error = function(e) NULL)
    cat("Running Linear Diversity Dependent geo model....\n")
    fit_DDlin_geo <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDlin", geography.object=geo_obj) }, error = function(e) NULL)
    cat("Running Exponential Diversity Dependent geo model....\n")
    fit_DDexp_geo <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDexp", geography.object=geo_obj) }, error = function(e) NULL)
  } else {
    fit_MC_geo <- NULL
    fit_DDlin_geo <- NULL
    fit_DDexp_geo <- NULL
  }
  
  fits <- list(MC = fit_MC, DDlin = fit_DDlin, DDexp = fit_DDexp, 
               MC_geo = fit_MC_geo, DDlin_geo = fit_DDlin_geo, DDexp_geo = fit_DDexp_geo)
  
  aics <- sapply(fits, function(fit) if (!is.null(fit)) fit$aic else NA)
  
  DD_evo_results[[paste0("colorPC", v)]] <- list(AIC = aics, best_model = fits[[which.min(aics)]])
  DD_evo_fits[[paste0("colorPC", v)]] <- fits
  
  saveRDS(DD_evo_results, "./rdata/uni_colorPC_DD_res.rds")
  saveRDS(DD_evo_fits, "./rdata/uni_colorPC_DD_models.rds")
  
}

# morpho
trait_evo_res <- list()
trait_evo_fits <- list()

for (v in 1:ncol(morphoPCs)) {
  trait_vector <- setNames(scale(morphoPCs[tree$tip.label, v]), tree$tip.label)
  
  # Fit BM, WN, OU
  cat("BM...\n")
  fit_BM <- fitContinuous(tree, trait_vector, model = "BM")
  cat("WN...\n")
  fit_WN <- fitContinuous(tree, trait_vector, model = "white")
  cat("OU...\n")
  fit_OU <- fitContinuous(tree, trait_vector, model = "OU")
  
  fit_OUM_list <- list()
  for (trait in discrete_traits) {
    asr_path <- paste0("./rdata/ASR_", trait, ".rda")
    if (file.exists(asr_path)) {
      load(asr_path)  # Should load `fit.marginal`
      
      # Format data for OUwie (assumes fit.marginal is loaded properly)
      ouwie_data <- data.frame(fit.marginal$data.legend, 
                               X = trait_vector[fit.marginal$data.legend$sp])
      ouwie_data <- ouwie_data[fit.marginal$phy$tip.label,]
      
      cat("OUM:", trait, "...\n")
      fit_OUM_list[[paste0("OUM_", trait)]] <- OUwie(fit.marginal$phy, ouwie_data, model = "OUM", algorithm = "invert", get.root.theta = FALSE)
      
    } else {
      warning("Trait rda file not found!", immediate. = TRUE)
    }
  }
  # Collect AICs and model fits
  aic_vec <- c(BM = fit_BM$opt$aic, WN = fit_WN$opt$aic, OU = fit_OU$opt$aic, sapply(fit_OUM_list, function(m) m$AIC))
  models <- c(list(BM = fit_BM, WN = fit_WN, OU = fit_OU), fit_OUM_list)
  
  trait_evo_res[[paste0("morphoPC", v)]] <- list(AIC = aic_vec, best_model = models[[which.min(aic_vec)]])
  trait_evo_fits[[paste0("morphoPC", v)]] <- models
}

# Save results per trait
saveRDS(trait_evo_res, "./rdata/uni_morphoPC_res.rds")
saveRDS(trait_evo_fits, "./rdata/uni_morphoPC_evo_model.rds")

DD_evo_results <- list()
DD_evo_fits <- list()

for (v in 1:8) {
  cat(sprintf("Working on morpho PC%s...\n", v))
  
  trait_vector <- setNames(scale(morphoPCs[tree$tip.label, v])[,1], tree$tip.label)
  
  # Fit models with error handling
  cat("MC model...\n")
  fit_MC <- tryCatch({ fit_t_comp(tree, trait_vector, model = "MC") }, error = function(e) NULL)
  cat("DDlin model...\n")
  fit_DDlin <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDlin") }, error = function(e) NULL)
  cat("DDexp model...\n")
  fit_DDexp <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDexp") }, error = function(e) NULL)
  
  if ( geo_model) {
    # Fit models with error handling
    fit_MC_geo <- tryCatch({ fit_t_comp(tree, trait_vector, model = "MC", geography.object=geo_obj) }, error = function(e) NULL)
    fit_DDlin_geo <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDlin", geography.object=geo_obj) }, error = function(e) NULL)
    fit_DDexp_geo <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDexp", geography.object=geo_obj) }, error = function(e) NULL)
  } else {
    fit_MC_geo <- NULL
    fit_DDlin_geo <- NULL
    fit_DDexp_geo <- NULL
  }
  
  fits <- list(MC = fit_MC, DDlin = fit_DDlin, DDexp = fit_DDexp, 
               MC_geo = fit_MC_geo, DDlin_geo = fit_DDlin_geo, DDexp_geo = fit_DDexp_geo)
  
  aics <- sapply(fits, function(fit) if (!is.null(fit)) fit$aic else NA)
  
  DD_evo_results[[paste0("morphoPC", v)]] <- list(AIC = aics, best_model = fits[[which.min(aics)]])
  DD_evo_fits[[paste0("morphoPC", v)]] <- fits
  
  saveRDS(DD_evo_results, "./rdata/uni_morphoPC_DD_res.rds")
  saveRDS(DD_evo_fits, "./rdata/uni_morphoPC_DD_models.rds")
}

color_DD_evo_results <- readRDS("./rdata/uni_colorPC_DD_res.rds")
color_evo_results <- readRDS("./rdata/uni_colorPC_res.rds")

sapply(1:length(color_evo_results), function(i) {
  res <- c(color_evo_results[[i]]$AIC, color_DD_evo_results[[i]]$AIC)
  idx <- which.min(res)
  res[idx]
})

morpho_DD_evo_results <- readRDS("./rdata/uni_morphoPC_DD_res.rds")
morpho_evo_results <- readRDS("./rdata/uni_morphoPC_res.rds")

sapply(1:length(morpho_evo_results), function(i) {
  res <- c(morpho_evo_results[[i]]$AIC)#, morpho_DD_evo_results[[i]]$AIC)
  idx <- which.min(res)
  res[idx]
}) 

rm(DD_evo_results, DD_evo_fits, trait_evo_res, trait_evo_fits)
gc()

### Regional models
# control for tree transformation getting tree depth, PD and other metrics to evaluate tree differences...
DD_evo_reg_results <- list()
DD_evo_reg_fits <- list()
region_stats <- list()
for (v in 1:ncol(colorPCs)) {
  trait_vector <- setNames(scale(colorPCs[tree$tip.label, v])[,1], tree$tip.label)

  region_stats_df <- data.frame()
  fit_reg_List <- list()
  aic_reg_List <- list()
  best_models <- list()
  for (r in colnames(pomacentridae_regs)) {
    
    reg_species <- rownames(pomacentridae_regs)[pomacentridae_regs[,r] == "1"]
    if (length(reg_species) > 1) {
      tree_reg <- ape::keep.tip(tree, reg_species)
      trait_vector_reg <- trait_vector[reg_species]
      
      # Tree-based metrics
      reg_ntip <- Ntip(tree_reg)
      reg_tree_depth <- max(ape::branching.times(tree_reg))
      reg_PD <- picante::pd(samp = as.data.frame(matrix(1, nrow=1, ncol=length(tree_reg$tip.label),
                                              dimnames = list(NULL, tree_reg$tip.label))),
                            tree = tree_reg)
      dist_mat <- cophenetic(tree_reg)
      reg_mpd <- mean(dist_mat[upper.tri(dist_mat)])
      reg_mntd <- mean(apply(dist_mat + diag(Inf, nrow(dist_mat)), 1, min))
      reg_TI <- phyloTop::colless.phylo(tree_reg, normalise = TRUE)
    
      # Trait-based metrics
      trait_mean <- mean(trait_vector_reg) 
      trait_range <- range(trait_vector_reg)
      trait_var <- var(trait_vector_reg, na.rm = TRUE)
      signal_K <- tryCatch({
        phylosig(tree_reg, trait_vector_reg, method = "K", test = TRUE)
      }, error = function(e) NA)

      region_stats_df <- rbind(region_stats_df, data.frame(
        region = region_names[r],
        n_species = reg_ntip, 
        depth = reg_tree_depth,
        PD = reg_PD$PD,
        MPD = reg_mpd,
        MNTD = reg_mntd,
        imbalance = reg_TI,
        trait_mean = trait_mean,
        trait_min = trait_range[1],
        trait_max = trait_range[2],
        trait_var = trait_var,
        K = signal_K$K,
        K_pval = signal_K$P
      ))

      # Fit BM, WN, OU
      cat("BM...\n")
      fit_BM <- fitContinuous(tree_reg, trait_vector_reg, model = "BM")
      cat("WN...\n")
      fit_WN <- fitContinuous(tree_reg, trait_vector_reg, model = "white")
      cat("OU...\n")
      fit_OU <- fitContinuous(tree_reg, trait_vector_reg, model = "OU", bounds = list(alpha=c(1e-8, 100)))

      fit_OUM_list <- list()
      for (trait in discrete_traits) {
        asr_path <- paste0("./rdata/ASR_", trait, ".rda")
        if (file.exists(asr_path)) {
          load(asr_path)  # Should load `fit.marginal`

          # Format data for OUwie (assumes fit.marginal is loaded properly)
          ouwie_data <- data.frame(fit.marginal$data.legend,
                                   X = trait_vector[fit.marginal$data.legend$sp])
          # subset
          ouwie_data_reg <- ouwie_data[reg_species,]
          fit.marginal_reg <- ape::keep.tip(fit.marginal$phy, reg_species)

          cat("OUM:", trait, "...\n")
          fit_OUM_list[[paste0("OUM_", trait)]] <- OUwie(fit.marginal_reg, ouwie_data_reg, model = "OUM",
                                                         algorithm = "invert", get.root.theta = FALSE)

        } else {
          warning("Trait rda file not found!", immediate. = TRUE)
        }
      }

      # Fit denso-dependent models with error handling
      fit_MC <- tryCatch({ fit_t_comp(tree_reg, trait_vector_reg, model = "MC") }, error = function(e) NULL)
      fit_DDlin <- tryCatch({ fit_t_comp(tree_reg, trait_vector_reg, model = "DDlin") }, error = function(e) NULL)
      fit_DDexp <- tryCatch({ fit_t_comp(tree_reg, trait_vector_reg, model = "DDexp") }, error = function(e) NULL)

      # Collect AICs and model fits
      fits_reg <- c(list(BM = fit_BM, WN = fit_WN, OU = fit_OU),
                    fit_OUM_list,
                    list(MC = fit_MC, DDlin = fit_DDlin, DDexp = fit_DDexp))

      aics_reg <- c(BM = fit_BM$opt$aic, WN = fit_WN$opt$aic, OU = fit_OU$opt$aic,
                    sapply(fit_OUM_list, function(m) m$AIC),
                    sapply(fits_reg[c("MC", "DDlin", "DDexp")], function(fit) if (!is.null(fit)) fit$aic else NA))

      fit_reg_List[[r]] <- fits_reg
      aic_reg_List[[r]] <- aics_reg
      best_models[[r]] <- c(model = names(which.min(aics_reg)), fits_reg[[names(which.min(aics_reg))]])
    }
  }
  region_stats[[paste0("colorPC", v)]] <- region_stats_df
  DD_evo_reg_results[[paste0("colorPC", v)]] <- list(AIC = aic_reg_List, best_models = best_models)
  DD_evo_reg_fits[[paste0("colorPC", v)]] <- fit_reg_List
}

region_stats_df <- plyr::ldply(region_stats, .id = "trait")

saveRDS(region_stats_df, "./rdata/uni_colorPC_tree_region_stats.rds")
saveRDS(DD_evo_reg_results, "./rdata/uni_colorPC_DD_reg_res.rds")
saveRDS(DD_evo_reg_fits, "./rdata/uni_colorPC_DD_reg_models.rds")

DD_evo_reg_results <- list()
DD_evo_reg_fits <- list()
region_stats <- list()
for (v in 1:8) {
  trait_vector <- setNames(scale(morphoPCs[tree$tip.label, v])[,1], tree$tip.label)
  
  region_stats_df <- data.frame()
  fit_reg_List <- list()
  aic_reg_List <- list()
  best_models <- list()
  for (r in colnames(pomacentridae_regs)) {
    
    reg_species <- rownames(pomacentridae_regs)[pomacentridae_regs[,r] == "1"]
    if (length(reg_species) > 1) {
      tree_reg <- ape::keep.tip(tree, reg_species)
      trait_vector_reg <- trait_vector[reg_species]
      
      # Tree-based metrics
      reg_ntip <- Ntip(tree_reg)
      reg_tree_depth <- max(ape::branching.times(tree_reg))
      reg_PD <- picante::pd(samp = as.data.frame(matrix(1, nrow=1, ncol=length(tree_reg$tip.label),
                                                        dimnames = list(NULL, tree_reg$tip.label))),
                            tree = tree_reg)
      dist_mat <- cophenetic(tree_reg)
      reg_mpd <- mean(dist_mat[upper.tri(dist_mat)])
      reg_mntd <- mean(apply(dist_mat + diag(Inf, nrow(dist_mat)), 1, min))
      reg_TI <- phyloTop::colless.phylo(tree_reg, normalise = TRUE)
      
      # Trait-based metrics
      trait_mean <- mean(trait_vector_reg) 
      trait_range <- range(trait_vector_reg)
      trait_var <- var(trait_vector_reg, na.rm = TRUE)
      signal_K <- tryCatch({
        phylosig(tree_reg, trait_vector_reg, method = "K", test = TRUE)
      }, error = function(e) NA)
      
      region_stats_df <- rbind(region_stats_df, data.frame(
        region = region_names[r],
        n_species = reg_ntip, 
        depth = reg_tree_depth,
        PD = reg_PD$PD,
        MPD = reg_mpd,
        MNTD = reg_mntd,
        imbalance = reg_TI,
        trait_mean = trait_mean,
        trait_min = trait_range[1],
        trait_max = trait_range[2],
        trait_var = trait_var,
        K = signal_K$K,
        K_pval = signal_K$P
      ))
      
      # Fit BM, WN, OU
      cat("BM...\n")
      fit_BM <- fitContinuous(tree_reg, trait_vector_reg, model = "BM")
      cat("WN...\n")
      fit_WN <- fitContinuous(tree_reg, trait_vector_reg, model = "white")
      cat("OU...\n")
      fit_OU <- fitContinuous(tree_reg, trait_vector_reg, model = "OU", bounds = list(alpha=c(1e-8, 100)))

      fit_OUM_list <- list()
      for (trait in discrete_traits) {
        asr_path <- paste0("./rdata/ASR_", trait, ".rda")
        if (file.exists(asr_path)) {
          load(asr_path)  # Should load `fit.marginal`

          # Format data for OUwie (assumes fit.marginal is loaded properly)
          ouwie_data <- data.frame(fit.marginal$data.legend,
                                   X = trait_vector[fit.marginal$data.legend$sp])
          # subset
          ouwie_data_reg <- ouwie_data[reg_species,]
          fit.marginal_reg <- ape::keep.tip(fit.marginal$phy, reg_species)

          cat("OUM:", trait, "...\n")
          fit_OUM_list[[paste0("OUM_", trait)]] <- OUwie(fit.marginal_reg, ouwie_data_reg, model = "OUM",
                                                         algorithm = "invert", get.root.theta = FALSE)

        } else {
          warning("Trait rda file not found!", immediate. = TRUE)
        }
      }

      # Fit denso-dependent models with error handling
      fit_MC <- tryCatch({ fit_t_comp(tree_reg, trait_vector_reg, model = "MC") }, error = function(e) NULL)
      fit_DDlin <- tryCatch({ fit_t_comp(tree_reg, trait_vector_reg, model = "DDlin") }, error = function(e) NULL)
      fit_DDexp <- tryCatch({ fit_t_comp(tree_reg, trait_vector_reg, model = "DDexp") }, error = function(e) NULL)

      # Collect AICs and model fits
      fits_reg <- c(list(BM = fit_BM, WN = fit_WN, OU = fit_OU),
                    fit_OUM_list,
                    list(MC = fit_MC, DDlin = fit_DDlin, DDexp = fit_DDexp))

      aics_reg <- c(BM = fit_BM$opt$aic, WN = fit_WN$opt$aic, OU = fit_OU$opt$aic,
                    sapply(fit_OUM_list, function(m) m$AIC),
                    sapply(fits_reg[c("MC", "DDlin", "DDexp")], function(fit) if (!is.null(fit)) fit$aic else NA))

      fit_reg_List[[r]] <- fits_reg
      aic_reg_List[[r]] <- aics_reg
      best_models[[r]] <- c(model = names(which.min(aics_reg)), fits_reg[[names(which.min(aics_reg))]])
    }
  }
  
  region_stats[[paste0("morphoPC", v)]] <- region_stats_df
  DD_evo_reg_results[[paste0("morphoPC", v)]] <- list(AIC = aic_reg_List, best_models = best_models)
  DD_evo_reg_fits[[paste0("morphoPC", v)]] <- fit_reg_List
}

region_stats_df <- plyr::ldply(region_stats, .id = "trait")

saveRDS(region_stats_df, "./rdata/uni_morphoPC_tree_region_stats.rds")
saveRDS(DD_evo_reg_results, "./rdata/uni_morphoPC_DD_reg_res.rds")
saveRDS(DD_evo_reg_fits, "./rdata/uni_morphoPC_DD_reg_models.rds")

