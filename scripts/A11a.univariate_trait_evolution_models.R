# =============================
#  Libraries
# =============================
libs <- c("ape","picante","sp","geiger","wesanderson","OUwie",
          "mvMORPH","phytools","phyloTop","corHMM","parallel","tibble")
invisible(sapply(libs, require, character.only = TRUE))

setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")
dir.create("rdata", showWarnings = FALSE)
dir.create("results/univariate", recursive = TRUE)
dir.create("figures/univariate", recursive = TRUE)

ncores <- 5
geo_model <- FALSE
source("./scripts/A1.Input_data.R")

pomacentridae_regs <- BioGeoBEARS::getranges_from_LagrangePHYLIP(
  "./data/CB2013_pomacentridae_geodata_input.txt"
)@df
region_names <- setNames(c("IO","IAA","CPO","EPO","AO","TS"), LETTERS[1:6])


# =============================
#  Helper functions
# =============================

fit_trait_models <- function(tree, trait_vector, discrete_traits) {
  # Fit baseline models
  fits <- list(
    BM  = fitContinuous(tree, trait_vector, model = "BM"),
    WN  = fitContinuous(tree, trait_vector, model = "white"),
    EB  = fitContinuous(tree, trait_vector, model = "EB"),
    OU  = fitContinuous(tree, trait_vector, model = "OU")
  )
  
  # Fit OUM by discrete traits
  for (trait in discrete_traits) {
    asr_path <- paste0("./rdata/ASR_", trait, ".rda")
    if (!file.exists(asr_path)) next
    load(asr_path) # loads fit.marginal
    ouwie_data <- data.frame(fit.marginal$data.legend,
                             X = trait_vector[fit.marginal$data.legend$sp])
    ouwie_data <- ouwie_data[fit.marginal$phy$tip.label,]
    fits[[paste0("OUM_", trait)]] <- tryCatch(
      OUwie(fit.marginal$phy, ouwie_data, model = "OUM",
            algorithm = "invert", get.root.theta = FALSE),
      error = function(e) NULL
    )
    rm(ouwie_data, fit.marginal); gc()
  }
  fits
}


fit_DD_models <- function(tree, trait_vector, geo_model = FALSE) {
  base_models <- c("MC","DDlin","DDexp")
  fits <- lapply(base_models, function(m)
    tryCatch(fit_t_comp(tree, trait_vector, model = m), error = function(e) NULL))
  names(fits) <- base_models
  
  if (geo_model && exists("geo_obj")) {
    geo_fits <- lapply(base_models, function(m)
      tryCatch(fit_t_comp(tree, trait_vector, model = m, geography.object = geo_obj),
               error = function(e) NULL))
    names(geo_fits) <- paste0(base_models, "_geo")
    fits <- c(fits, geo_fits)
    rm(geo_fits); gc()
  }
  fits
}


compute_region_stats <- function(tree_reg, trait_vector_reg, region_code) {
  reg_ntip <- Ntip(tree_reg)
  reg_tree_depth <- max(branching.times(tree_reg))
  reg_PD <- picante::pd(
    samp = as.data.frame(matrix(1, nrow=1, ncol=reg_ntip,
                                dimnames = list(NULL, tree_reg$tip.label))),
    tree = tree_reg)
  dist_mat <- cophenetic(tree_reg)
  reg_mpd  <- mean(dist_mat[upper.tri(dist_mat)])
  reg_mntd <- mean(apply(dist_mat + diag(Inf, nrow(dist_mat)), 1, min))
  reg_TI   <- phyloTop::colless.phylo(tree_reg, normalise = TRUE)
  
  signal_K <- tryCatch(phylosig(tree_reg, trait_vector_reg, method="K", test=TRUE),
                       error = function(e) list(K=NA, P=NA))
  
  data.frame(
    region = region_names[region_code],
    n_species = reg_ntip,
    depth = reg_tree_depth,
    PD = reg_PD$PD,
    MPD = reg_mpd,
    MNTD = reg_mntd,
    imbalance = reg_TI,
    trait_mean = mean(trait_vector_reg),
    trait_min = min(trait_vector_reg),
    trait_max = max(trait_vector_reg),
    trait_var = var(trait_vector_reg),
    K = signal_K$K,
    K_pval = signal_K$P
  )
}


fit_models_per_region <- function(tree, trait_vector, pomacentridae_regs, discrete_traits) {
  reg_fits <- list()
  reg_aics <- list()
  reg_stats <- data.frame()
  
  for (r in colnames(pomacentridae_regs)) {
    reg_species <- rownames(pomacentridae_regs)[pomacentridae_regs[,r] == "1"]
    if (length(reg_species) <= 1) next
    
    tree_reg <- keep.tip(tree, reg_species)
    trait_vector_reg <- trait_vector[reg_species]
    
    # Tree and trait stats
    reg_stats <- rbind(reg_stats, compute_region_stats(tree_reg, trait_vector_reg, r))
    
    # Fit models
    fits <- c(fit_trait_models(tree_reg, trait_vector_reg, discrete_traits),
              fit_DD_models(tree_reg, trait_vector_reg, geo_model))
    reg_fits[[r]] <- fits
    reg_aics[[r]] <- sapply(fits, function(fit)
      if (is.null(fit)) NA else if (!is.null(fit$opt$aic)) fit$opt$aic else fit$AIC)
    
    # clean small objects per region
    rm(tree_reg, trait_vector_reg, fits); gc(verbose = FALSE)
  }
  
  list(fits = reg_fits, AIC = reg_aics, stats = reg_stats)
}


# =============================
#  Global trait model fitting
# =============================
run_global_trait_models <- function(trait_matrix, label) {
  
  tmp_file <- paste0("./rdata/uni_", tolower(label), "PC_global_traitevo_models_partial.rds")
  rds_file <- paste0("./rdata/uni_", tolower(label), "PC_global_traitevo_models.rds")

  fits_list <- list()
  for (v in seq_len(ncol(trait_matrix))) {
    cat("Running global fits for", label, "PC", v, "...\n")
    trait_vector <- setNames(scale(trait_matrix[tree$tip.label, v])[,1], tree$tip.label)
    fits <- c(fit_trait_models(tree, trait_vector, discrete_traits),
              fit_DD_models(tree, trait_vector, geo_model))
    fits_list[[paste0(label, "PC", v)]] <- fits
    
    # Incremental save and cleanup
    saveRDS(fits_list, file = tmp_file)
    rm(fits, trait_vector); gc(verbose = FALSE)
  }
  
  saveRDS(fits_list, file = rds_file)
  if (file.exists(rds_file)) file.remove(tmp_file)
  rm(fits_list); gc()
}


# =============================
#  Regional model fitting
# =============================
run_regional_trait_models <- function(trait_matrix, label) {

  csv_file_stats <- paste0("./rdata/uni_", tolower(label), "PC_tree_region_stats.csv")
  rds_file_stats <- paste0("./rdata/uni_", tolower(label), "PC_tree_region_stats.rds")
  
  tmp_file <- paste0("./rdata/uni_", tolower(label), "PC_reg_traitevo_models_partial.rds")
  rds_file <- paste0("./rdata/uni_", tolower(label), "PC_reg_traitevo_models.rds")

  # Initialize CSV file with header (once)
  if (file.exists(csv_file_stats)) file.remove(csv_file_stats)
  write.table(
    data.frame(), csv_file_stats,
    sep = ",", row.names = FALSE, col.names = TRUE
  )
  
  reg_results <- list()
  reg_stats_list <- list()
  for (v in seq_len(ncol(trait_matrix))) {
    cat("Running regional fits for", label, "PC", v, "...\n")
    trait_vector <- setNames(scale(trait_matrix[tree$tip.label, v])[,1], tree$tip.label)
    
    # Run model fits per region
    res <- fit_models_per_region(tree, trait_vector, pomacentridae_regs, discrete_traits)
    
    # Store results
    reg_results[[paste0(label, "PC", v)]] <- list(AIC = res$AIC, fits = res$fits)
    reg_stats_list[[paste0(label, "PC", v)]] <- res$stats
    
    # Add PC label before writing
    res$stats$trait <- paste0(label, "PC", v)
    
    # Append new rows to CSV incrementally
    suppressWarnings(
      write.table(
        res$stats, file = csv_file_stats,
        sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE
      )
    )
    
    # Incremental save and cleanup
    saveRDS(reg_results, file = tmp_file)
    rm(res, trait_vector)
    gc(verbose = FALSE)
  }
  
  saveRDS(reg_results, file = rds_file)
  saveRDS(do.call(rbind, reg_stats_list), rds_file_stats)
  
  cat("✔ Regional model runs completed for", label, 
      "(", ncol(trait_matrix), "PCs ).\n")
  
  if (file.exists(rds_file)) file.remove(tmp_file)
  rm(reg_results, reg_stats_list); gc()
}

# =============================
#  Run all analyses
# =============================
run_global_trait_models(colorPCs,  "color")
run_global_trait_models(morphoPCs, "morpho")
run_regional_trait_models(colorPCs,  "color")
run_regional_trait_models(morphoPCs, "morpho")

gc()
