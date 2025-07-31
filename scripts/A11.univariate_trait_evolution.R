lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", "phyloTop",
         "lme4", "mgcv", "phytools", "corHMM", "ggpmisc",
         "lmtest", "ColorAR", "castor", "tibble",
         "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

#=============================#
#   Settings
#=============================#
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/")

# load custom functions
source("./Rscripts/custom_functions.R")

res_folder <- "new_results"
dir.create(res_folder)
dir.create("new_figures/univariate/")

col.palette <- sample(unique(as.character(sapply(c("Zissou1", "Darjeeling1", "Darjeeling2", "FantasticFox1"), function(x) wes_palette(x, 5)))))
ncores <- 5
geo_model = FALSE

# load input data and custom functions
source("./Rscripts/A0.Input_data.R")

pomacentridae_regs <- BioGeoBEARS::getranges_from_LagrangePHYLIP("./BioGeo/input/CB2021_pomacentridae_geodata_input.txt")@df

#=============================#
# Trait Evolution Models
#=============================#

# color
trait_evo_res <- list()
trait_evo_fits <- list()

for (v in 1:ncol(colorPCs)) {
  trait_vector <- setNames(scale(colorPCs[tree$tip.label, v]), tree$tip.label)
  
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
  trait_vector <- setNames(scale(colorPCs[tree$tip.label, v]), tree$tip.label)
  
  # Fit models with error handling
  fit_MC <- tryCatch({ fit_t_comp(tree, trait_vector, model = "MC") }, error = function(e) NULL)
  fit_DDlin <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDlin") }, error = function(e) NULL)
  fit_DDexp <- tryCatch({ fit_t_comp(tree, trait_vector, model = "DDexp") }, error = function(e) NULL)
  
  if (geo_model) {
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

### Regional models
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

pomacentridae_regs <- BioGeoBEARS::getranges_from_LagrangePHYLIP("./BioGeo/input/CB2021_pomacentridae_geodata_input.txt")@df
region_names <- setNames(area_names <- c("IO", "IAA", "CPO", "EPO", "AO", "TS"), LETTERS[1:6])

# control for tree transformation getting tree depth, PD and other metrics to evaluate tree differences...
DD_evo_reg_results <- list()
DD_evo_reg_fits <- list()
region_stats <- list()
for (v in 1:ncol(colorPCs)) {
  trait_vector <- setNames(scale(colorPCs[tree$tip.label, v]), tree$tip.label)

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
      reg_TI <- colless.phylo(tree_reg, normalise = TRUE)
    
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
  trait_vector <- setNames(scale(morphoPCs[tree$tip.label, v]), tree$tip.label)
  
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
      reg_TI <- colless.phylo(tree_reg, normalise = TRUE)
      
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

### ================================= ### 
###           Review results          ###
### ================================= ### 
region_names <- setNames(area_names <- c("IO", "IAA", "CPO", "EPO", "AO", "TS"), LETTERS[1:6])

col_region_stats_df <- readRDS("./rdata/uni_colorPC_tree_region_stats.rds")
for (i in 1:nrow(col_region_stats_df)) {
  trait <- col_region_stats_df$trait[i]
  region <- col_region_stats_df$region[i]
  best_model <- col_region_stats_df$best_model[i]
  reg_df <- color_res_global[[sub("color", "", trait)]][[2]]
  col_region_stats_df[i, "model_support"] <- reg_df$Weight[reg_df$region == as.character(region) & reg_df$DeltaAIC == 0]
}

col_region_stats_long <- tidyr::gather(col_region_stats_df, metric, value, 3:13)


# Custom labels for facets
metric_labels <- c(
  #"depth" = "Tree depth",
  "imbalance" = "Tree imbalance",
  "K" = "Phylogenetic signal (Blomberg's K)",
  "MNTD" = "Mean Nearest Taxon Distance (MNTD)",
  "MPD" = "Mean Pairwise Distance (MPD)",
  "PD" = "Phylogenetic Diversity (PD)",
  "trait_var" = "Trait variance"
)

ggplot(col_region_stats_long[col_region_stats_long$metric %in% names(metric_labels), ],
       aes(x = value, y = model_support)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "skyblue", fill = "lightblue", se = TRUE) +
  scale_y_continuous(breaks = seq(0,1, 0.2)) +
  coord_cartesian(ylim = c(0,1.1)) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), ..p.value.label.., sep = "~~~")),
               formula = y ~ x,
               parse = TRUE,
               size = 3.5, color = "black") +
  facet_wrap(. ~ metric, scales = "free", labeller = labeller(metric = metric_labels)) +
  theme_light(base_size = 14) +
  labs(x = "Metric value", y = "Best Model Support (Akaike weight)") +
  theme(strip.background = element_rect(fill = "grey95", color = "grey80"),
        strip.text = element_text(size = 13, face = "bold", color = "grey20"),
        plot.margin = margin(35, 20, 25, 20),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20, vjust = 0, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(0.8, "lines"),
        panel.spacing = unit(1, "lines"),
  )
ggsave("./new_figures/Struct_bias_regional_model_support.pdf", width = 15, height = 10)


## Color

col_evo_results <- readRDS("./rdata/uni_colorPC_res.rds")
col_DD_evo_results <- readRDS("./rdata/uni_colorPC_DD_res.rds")
col_DD_evo_reg_results <- readRDS("./rdata/uni_colorPC_DD_reg_res.rds")

col_evo_models <- readRDS("./rdata/uni_colorPC_evo_models.rds")
col_DD_evo_models <- readRDS("./rdata/uni_colorPC_DD_models.rds")
col_DD_evo_reg_models <- readRDS("./rdata/uni_colorPC_DD_reg_models.rds")

col_region_stats_df <- readRDS("./rdata/uni_colorPC_tree_region_stats.rds")
for (i in 1:nrow(col_region_stats_df)) {
  trait <- col_region_stats_df$trait[i]
  region <- names(region_names)[region_names == col_region_stats_df$region[i]]
  best_model <- col_DD_evo_reg_results[[trait]]$best_models[[region]]$model
  AICc <- col_DD_evo_reg_results[[trait]]$best_models[[region]]$AICc
  if (best_model %in% c("WN", "OU")) {AICc <- col_DD_evo_reg_results[[trait]]$best_models[[region]]$opt$aicc}
  else{ 
    if (is.null(AICc)) {AICc <- col_DD_evo_reg_results[[trait]]$best_models[[region]]$aicc}
  }
  col_region_stats_df[i, "best_model"] = best_model
  col_region_stats_df[i, "AICc"] = AICc
}

col_region_stats_df$best_model[grep("OU", col_region_stats_df$best_model)] = "OU-type"
col_region_stats_df$best_model <- relevel(as.factor(col_region_stats_df$best_model), ref = "WN")
col_region_stats_df$region <- relevel(as.factor(col_region_stats_df$region), ref = "CPO")

mm_colPC <- nnet::multinom(best_model ~ region + scale(PD) + scale(K) + scale(trait_var), data = col_region_stats_df)
summary(mm_colPC)

unique(col_region_stats_df$best_model)
cor(col_region_stats_df[c("PD", "K", "trait_var", "depth", "n_species")])

color_res_global <- list()
for (pc in colnames(colorPCs)) {
  
  ## global models
  m_global <- c(col_evo_results[[paste0("color", pc)]]$AIC, col_DD_evo_results[[paste0("col", pc)]]$AIC)
  m_global_df <- data.frame(model = names(m_global), AIC = m_global)
  
  # rename some models for clarity
  m_global_df$model[m_global_df$model == "OUM_DietEcotype"] = "OUM_Diet"
  m_global_df$model[m_global_df$model == "OUM_Farming"] = "OUM_Symbiosis"
  m_global_df$model[m_global_df$model == "OUM_habitat1"] = "OUM_Habitat"
  
  # Compute delta AIC (ΔAIC)
  m_global_df$DeltaAIC <- m_global_df$AIC - min(m_global_df$AIC)
  
  # Compute Akaike weights
  m_global_df$Weight <- exp(-0.5 * m_global_df$DeltaAIC)
  m_global_df$Weight <- m_global_df$Weight / sum(m_global_df$Weight)
  
  # Sort by support
  m_global_df <- m_global_df[order(m_global_df$Weight, decreasing = TRUE), ]
  
  ggplot(m_global_df, aes(x = reorder(model, Weight), y = Weight)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = sprintf("Color %s Model Support", pc), x = "Model", y = "Akaike Weight") +
    scale_y_continuous(limits = c(0,1)) +
    theme_minimal(base_size = 14) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 20)),  
          axis.title.y = element_text(margin = margin(r = 20))) +
    coord_flip()
  ggsave(sprintf("./global_color%s_univariate_model_support.pdf", pc), width = 7, height = 7)
  
  ## regional models
  m_reg <- col_DD_evo_reg_results[[paste0("color", pc)]]$AIC
  m_reg_df <- data.frame()
  for (r in 1:length(m_reg)) {
    m_reg_r <- m_reg[[r]]
    m_reg_rdf <- data.frame(model = names(m_reg_r), AIC = m_reg_r, region = region_names[names(m_reg)[r]])
    
    # rename some models for clarity
    m_reg_rdf$model[m_reg_rdf$model == "OUM_DietEcotype"] = "OUM_Diet"
    m_reg_rdf$model[m_reg_rdf$model == "OUM_Farming"] = "OUM_Symbiosis"
    m_reg_rdf$model[m_reg_rdf$model == "OUM_habitat1"] = "OUM_Habitat"
    
    # Compute delta AIC (ΔAIC)
    m_reg_rdf$DeltaAIC <- m_reg_rdf$AIC - min(m_reg_rdf$AIC)
    
    # Compute Akaike weights
    m_reg_rdf$Weight <- exp(-0.5 * m_reg_rdf$DeltaAIC)
    m_reg_rdf$Weight <- m_reg_rdf$Weight / sum(m_reg_rdf$Weight)
    
    # Sort by support
    m_reg_rdf <- m_reg_rdf[order(m_reg_rdf$Weight, decreasing = TRUE), ]
    m_reg_df <- rbind(m_reg_df, m_reg_rdf)
  }
  
  ggplot(m_reg_df, aes(x = reorder(model, Weight), y = Weight)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = sprintf("Color %s Model Support", pc), x = "Model", y = "Akaike Weight") +
    scale_y_continuous(limits = c(0,1)) +
    theme_minimal(base_size = 14) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 20)),  
          axis.title.y = element_text(margin = margin(r = 20))) +
    coord_flip() +
    facet_wrap(.~region)
  ggsave(sprintf("./regional_color%s_univariate_model_support.pdf", pc), width = 7, height = 7)
  
  color_res_global[[pc]] <- list(m_global_df, m_reg_df)
  
}


## morpho

morpho_evo_results <- readRDS("./rdata/uni_morphoPC_res.rds")
morpho_DD_evo_results <- readRDS("./rdata/uni_morphoPC_DD_res.rds")
morpho_DD_evo_reg_results <- readRDS("./rdata/uni_morphoPC_DD_reg_res.rds")

morpho_evo_models <- readRDS("./rdata/uni_morphoPC_evo_models.rds")
morpho_DD_evo_models <- readRDS("./rdata/uni_morphoPC_DD_models.rds")
morpho_DD_evo_reg_models <- readRDS("./rdata/uni_morphoPC_DD_reg_models.rds")

morpho_region_stats_df <- readRDS("./rdata/uni_morphoPC_tree_region_stats.rds")
for (i in 1:nrow(morpho_region_stats_df)) {
  trait <- morpho_region_stats_df$trait[i]
  region <- names(region_names)[region_names == morpho_region_stats_df$region[i]]
  best_model <- morpho_DD_evo_reg_results[[trait]]$best_models[[region]]$model
  AICc <- morpho_DD_evo_reg_results[[trait]]$best_models[[region]]$AICc
  if (best_model %in% c("WN", "OU")) {AICc <- morpho_DD_evo_reg_results[[trait]]$best_models[[region]]$opt$aicc}
  else{ 
    if (is.null(AICc)) {AICc <- morpho_DD_evo_reg_results[[trait]]$best_models[[region]]$aicc}
  }
  morpho_region_stats_df[i, "best_model"] = best_model
  morpho_region_stats_df[i, "AICc"] = AICc
}

morpho_region_stats_df$best_model[grep("OU", morpho_region_stats_df$best_model)] = "OU-type"
morpho_region_stats_df$best_model <- relevel(as.factor(morpho_region_stats_df$best_model), ref = "WN")
morpho_region_stats_df$region <- relevel(as.factor(morpho_region_stats_df$region), ref = "CPO")

mm_colPC <- nnet::multinom(best_model ~ region + scale(PD) + scale(K) + scale(trait_var), data = morpho_region_stats_df)
summary(mm_colPC)

unique(morpho_region_stats_df$best_model)
cor(morpho_region_stats_df[c("PD", "K", "trait_var", "depth", "n_species")])

morpho_res_global <- list()
for (pc in colnames(morphoPCs)[1:8]) {
  
  ## global models
  m_global <- c(morpho_evo_results[[paste0("morpho", pc)]]$AIC, morpho_DD_evo_results[[paste0("col", pc)]]$AIC)
  m_global_df <- data.frame(model = names(m_global), AIC = m_global)
  
  # rename some models for clarity
  m_global_df$model[m_global_df$model == "OUM_DietEcotype"] = "OUM_Diet"
  m_global_df$model[m_global_df$model == "OUM_Farming"] = "OUM_Symbiosis"
  m_global_df$model[m_global_df$model == "OUM_habitat1"] = "OUM_Habitat"
  
  # Compute delta AIC (ΔAIC)
  m_global_df$DeltaAIC <- m_global_df$AIC - min(m_global_df$AIC)
  
  # Compute Akaike weights
  m_global_df$Weight <- exp(-0.5 * m_global_df$DeltaAIC)
  m_global_df$Weight <- m_global_df$Weight / sum(m_global_df$Weight)
  
  # Sort by support
  m_global_df <- m_global_df[order(m_global_df$Weight, decreasing = TRUE), ]
  
  ggplot(m_global_df, aes(x = reorder(model, Weight), y = Weight)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = sprintf("morpho %s Model Support", pc), x = "Model", y = "Akaike Weight") +
    scale_y_continuous(limits = c(0,1)) +
    theme_minimal(base_size = 14) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 20)),  
          axis.title.y = element_text(margin = margin(r = 20))) +
    coord_flip()
  ggsave(sprintf("./global_morpho%s_univariate_model_support.pdf", pc), width = 7, height = 7)
  
  ## regional models
  m_reg <- morpho_DD_evo_reg_results[[paste0("morpho", pc)]]$AIC
  m_reg_df <- data.frame()
  for (r in 1:length(m_reg)) {
    m_reg_r <- m_reg[[r]]
    m_reg_rdf <- data.frame(model = names(m_reg_r), AIC = m_reg_r, region = region_names[names(m_reg)[r]])
    
    # rename some models for clarity
    m_reg_rdf$model[m_reg_rdf$model == "OUM_DietEcotype"] = "OUM_Diet"
    m_reg_rdf$model[m_reg_rdf$model == "OUM_Farming"] = "OUM_Symbiosis"
    m_reg_rdf$model[m_reg_rdf$model == "OUM_habitat1"] = "OUM_Habitat"
    
    # Compute delta AIC (ΔAIC)
    m_reg_rdf$DeltaAIC <- m_reg_rdf$AIC - min(m_reg_rdf$AIC)
    
    # Compute Akaike weights
    m_reg_rdf$Weight <- exp(-0.5 * m_reg_rdf$DeltaAIC)
    m_reg_rdf$Weight <- m_reg_rdf$Weight / sum(m_reg_rdf$Weight)
    
    # Sort by support
    m_reg_rdf <- m_reg_rdf[order(m_reg_rdf$Weight, decreasing = TRUE), ]
    m_reg_df <- rbind(m_reg_df, m_reg_rdf)
  }
  
  ggplot(m_reg_df, aes(x = reorder(model, Weight), y = Weight)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = sprintf("morpho %s Model Support", pc), x = "Model", y = "Akaike Weight") +
    scale_y_continuous(limits = c(0,1)) +
    theme_minimal(base_size = 14) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 20)),  
          axis.title.y = element_text(margin = margin(r = 20))) +
    coord_flip() +
    facet_wrap(.~region)
  ggsave(sprintf("./regional_morpho%s_univariate_model_support.pdf", pc), width = 7, height = 7)
  
  morpho_res_global[[pc]] <- list(m_global_df, m_reg_df)
  
}
