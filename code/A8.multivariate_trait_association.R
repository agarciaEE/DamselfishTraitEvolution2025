#=============================#
# Libraries
#=============================#
lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", "knitr",
         "lme4", "mgcv", "phytools", "corHMM", "kableExtra",
         "lmtest", "ColorAR", "castor", "mvMORPH",
         "RPANDA","OUwie", "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

#=============================#
# Settings
#=============================#
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

res_folder <- "results"
dir.create(res_folder)
dir.create("results/multivariate", recursive = TRUE)
dir.create("figures/multivariate", recursive = TRUE)

col.palette <- sample(unique(as.character(sapply(c("Zissou1", "Darjeeling1", "Darjeeling2", "FantasticFox1"), function(x) wes_palette(x, 5)))))
ncores <- 5
color_analysis = TRUE
morpho_analysis = TRUE

# load input data and custom functions
source("./scripts/A1.Input_data.R")

## Multivariate analyses testing colPCA and morphoPCA association to a trait
#################################################################################
## helper: extract the name and value of the process parameter for an mvgls fit
.mv_param <- function(mod, mname) {
  if (!inherits(mod, "mvgls") || is.null(mod$param)) {
    return(list(param_name = NA_character_, param_value = NA_real_))
  }
  if (grepl("lambda", mname, ignore.case = TRUE)) {
    return(list(param_name = "lambda", param_value = as.numeric(mod$param)))
  }
  if (grepl("^OU|_OU", mname, ignore.case = TRUE)) {
    return(list(param_name = "alpha",  param_value = as.numeric(mod$param)))
  }
  if (grepl("^EB|_EB", mname, ignore.case = TRUE)) {
    return(list(param_name = "a",      param_value = as.numeric(mod$param)))
  }
  # BM has no process parameter
  list(param_name = NA_character_, param_value = NA_real_)
}

# =======
# COLOUR
# =======
dat <- list(color = as.matrix(colorPCs))

model_BM1 <- mvgls(color ~ 1, data = dat, tree = tree, model = "BM")
model_lambda1 <- mvgls(color ~ 1, data = dat, tree = tree, model = "lambda")
u <- as.numeric(model_lambda1$start_values[2])  # Starting value for 'upper'
max_attempts <- 10
attempt <- 1
repeat {
  model_OU1 <- tryCatch(
    mvgls(color ~ 1, data = dat, tree = tree, model = "OU", upper = u),
    warning = function(w) "warning"
  )
  
  if (!is.character(model_OU1)) break  # Model fit successful, exit loop
  
  if (attempt >= max_attempts) {
    stop("Model failed after ", max_attempts, " attempts.")
  }
  
  cat("Warning caught: increasing 'upper' to", u, "\n")
  u <- u * 1.5  # Increase 'upper' by 50%
  attempt <- attempt + 1
}
model_EB1 <- mvgls(color~1, data = dat, tree = tree, model = "EB")

col_mvgls_models <- list(
  BM1 = model_BM1, 
  lambda1 = model_lambda1,
  OU1 = model_OU1, 
  EB1 = model_EB1
)

res <- lapply(col_mvgls_models, EIC, nbcores = ncores)

# ============================
## Color vs discrete traits
# ============================
for (trait in discrete_traits){
  if (file.exists(paste0("./rdata/ASR_", trait, ".rda"))){
    ## Distance matrices
    load(paste0("./rdata/ASR_", trait, ".rda"))
    
    trait.vec = setNames(eco_traits_subset[, trait], rownames(eco_traits_subset))[tree$tip.label]
    
    cat("Working on trait:", trait, "...\n")
    tree <- fit.marginal$phy
    tree$node.label <- levels(trait.vec)[tree$node.label]
    n <- Ntip(tree)
    
    scm_consensus <- make.consensus.simmap(trait_SCM_trees)
    
    cat("Running Multivariate Phylogenetic Generalized Least Squares on color PCA, testing ", trait, "...\n")
    
    ## MVPGLS
    dat <- list(color = as.matrix(colorPCs), trait = trait.vec)
    
    cat("BM...\n")
    model_BM_trait <- mvgls(color~trait, data = dat, tree = tree, model = "BM")
    
    cat("Lambda...\n")
    model_lambda_trait <- mvgls(color~trait, data = dat, tree = tree, model = "lambda")
    
    cat("OU...\n")
    u <- as.numeric(model_lambda1$start_values[2]) # Start 'upper' from a reasonable value
    attempt <- 1
    repeat {
      model_OU_trait <- tryCatch(
        mvgls(color ~ trait, data = dat, tree = tree, model = "OU", upper = u),
        warning = function(w) "warning"
      )
      
      if (!is.character(model_OU_trait)) break  # If no warning, model fit is successful
      
      if (attempt >= max_attempts) {
        stop("Model failed after ", max_attempts, " attempts.")
      }
      
      cat("Warning caught: increasing 'upper' to", u, "\n")
      u <- u * 1.5  # Increase upper limit by 50%
      attempt <- attempt + 1
    }
    
    cat("EB...\n")
    model_EB_trait <- mvgls(color~trait, data = dat, tree = tree, model = "EB")
    
    cat("BMM...\n")
    model_BMM1 <- mvgls(color~1, data = dat, tree = scm_consensus, model = "BMM")
    
    cat("OUM...\n")
    u <- as.numeric(model_lambda1$start_values[2]) # Start 'upper' from a reasonable value
    attempt <- 1
    repeat {
      model_OUM1 <- tryCatch(
        mvgls(color ~ 1, data = dat, tree = scm_consensus, model = "OUM", upper = u),
        warning = function(w) "warning"
      )
      
      if (!is.character(model_OUM1)) break  # If no warning, model fit is successful
      
      if (attempt >= max_attempts) {
        stop("Model failed after ", max_attempts, " attempts.")
      }
      
      cat("Warning caught: increasing 'upper' to", u, "\n")
      u <- u * 1.5  # Increase upper limit by 50%
      attempt <- attempt + 1
    }
    
    col_mvgls_models[[paste0("model_BM_", trait)]] <- model_BM_trait
    col_mvgls_models[[paste0("model_BMM_", trait)]] <- model_BMM1
    col_mvgls_models[[paste0("model_lambda_", trait)]] <- model_lambda_trait
    col_mvgls_models[[paste0("model_OU_", trait)]] <- model_OU_trait
    col_mvgls_models[[paste0("model_OUM_", trait)]] <- model_OUM1
    col_mvgls_models[[paste0("model_EB_", trait)]] <- model_EB_trait
    
    res <- c(res, lapply(col_mvgls_models[!names(col_mvgls_models) %in% names(res)], EIC, nbcores = ncores))
    
    # clean model objects
    rm(model_BMM1, model_OUM1, model_BM_trait, model_lambda_trait, model_OU_trait, model_EB_trait)

  } else {
    warning(trait, "rda file not found!! Skipping analysis...\n")
  }
}
saveRDS(col_mvgls_models, file = paste0("./rdata/colPCs_mvgls.models.rds"))

model_names <- names(col_mvgls_models)
param_list <- lapply(model_names, function(nm) .mv_param(col_mvgls_models[[nm]], nm))

mvgls.col.restab <- data.frame(model = model_names,
                               param = sapply(param_list, function(i) i[[1]]),
                               value = sapply(param_list, function(i) i[[2]]),
                               predictor = ifelse(
                                 str_detect(word(model_names, -1, sep = "_"), "DietEcotype"),
                                 paste0("DietEcotype"), 
                                 ifelse(str_detect(word(model_names, -1, sep = "_"), "Symbiosis"),
                                        paste0("Symbiosis"), 
                                        ifelse(str_detect(word(model_names, -1, sep = "_"), "Habitat"),
                                               paste0("Habitat"),"1"))),
                               logLik = sapply(res, function(i) i$LogLikelihood),
                               k = sapply(res, function(i) i$p),
                               se = sapply(res, function(i) i$se),
                               EIC = eic <- sapply(res, function(i) i$EIC),
                               delta.EIC = eic-min(eic),
                               weight = unclass(aic.w(eic)))

mvgls.col.restab$predictor[grepl("BMM|OUM", mvgls.col.restab$model)] <- "1"

best_model_name <- mvgls.col.restab$model[which.min(mvgls.col.restab$delta.EIC)]
best_model.col <- col_mvgls_models[[best_model_name]]

# Run MANOVA test
for (trait in discrete_traits) {
  best_trait_model <- mvgls.col.restab[grepl(trait, mvgls.col.restab$predictor),]
  manova_test <- sapply(best_trait_model$model, function(m) manova.gls(col_mvgls_models[[m]], 
                                                                       nperm = 999, 
                                                                       test = "Pillai", 
                                                                       verbose = TRUE), simplify = FALSE)
  mvgls.col.restab$Pillai[grepl(trait, mvgls.col.restab$predictor)] <- lapply(manova_test, function(i) i$stat)
  mvgls.col.restab$pvalue[grepl(trait, mvgls.col.restab$predictor)] <- sapply(manova_test, function(i) i$pvalue)
}

# write table to hard drive
cat("Writting discrete color models result table on hard disk...\n")
write.table(mvgls.col.restab, file = file.path(res_folder, paste0("multivariate/colorPCs_mvgls.models.fit.restable.txt")), 
            quote = FALSE, sep = "\t", row.names = T, col.names = T)

# remove all models
rm(list = ls(pattern = "^model_"))

# =======
# MORPHO
# =======
dat <- list(morpho = as.matrix(morphoPCs))

model_BM1 <- mvgls(morpho ~ 1, data = dat, tree = tree, model = "BM")
model_lambda1 <- mvgls(morpho ~ 1, data = dat, tree = tree, model = "lambda")
u <- as.numeric(model_lambda1$start_values[2]) # Start 'upper' from a reasonable value
max_attempts <- 10
attempt <- 1
repeat {
  model_OU1 <- tryCatch(
    mvgls(morpho ~ 1, data = dat, tree = tree, model = "OU", upper = u),
    warning = function(w) "warning"
  )
  
  if (!is.character(model_OU1)) break  # Model fit successful, exit loop
  
  if (attempt >= max_attempts) {
    stop("Model failed after ", max_attempts, " attempts.")
  }
  
  cat("Warning caught: increasing 'upper' to", u, "\n")
  u <- u * 1.5  # Increase 'upper' by 50%
  attempt <- attempt + 1
}
model_EB1 <- mvgls(morpho~1, data = dat, tree = tree, model = "EB")

morpho_mvgls_models <- list(
  BM1 = model_BM1, 
  lambda1 = model_lambda1,
  OU1 = model_OU1, 
  EB1 = model_EB1
)

res <- lapply(morpho_mvgls_models, EIC, nbcores = ncores)

# ==============================
## Morphology vs discrete traits
# ==============================
for (trait in discrete_traits){
  if (file.exists(paste0("./rdata/ASR_", trait, ".rda"))){
    ## Distance matrices
    load(paste0("./rdata/ASR_", trait, ".rda"))
    
    trait.vec = setNames(eco_traits_subset[, trait], rownames(eco_traits_subset))[tree$tip.label]
    
    cat("Working on trait:", trait, "...\n")
    tree <- fit.marginal$phy
    tree$node.label <- levels(trait.vec)[tree$node.label]
    n <- Ntip(tree)
    
    scm_consensus <- make.consensus.simmap(trait_SCM_trees)
    
  cat("Running Multivariate Phylogenetic Generalized Least Squares on morpho PCA, testing ", trait, "...\n")
  ## MVPGLS
  dat <- list(morpho = as.matrix(morphoPCs), trait = trait.vec)
  
  cat("BM...\n")
  model_BM_trait <- mvgls(morpho~trait, data = dat, tree = tree, model = "BM")
  
  cat("Lambda...\n")
  model_lambda_trait <- mvgls(morpho~trait, data = dat, tree = tree, model = "lambda")
  
  cat("OU...\n")
  u <- as.numeric(model_lambda1$start_values[2]) # Start 'upper' from a reasonable value
  attempt <- 1
  repeat {
    model_OU_trait <- tryCatch(
      mvgls(morpho ~ trait, data = dat, tree = tree, model = "OU", upper = u),
      warning = function(w) "warning"
    )
    
    if (!is.character(model_OU_trait)) break  # If no warning, model fit is successful
    
    if (attempt >= max_attempts) {
      stop("Model failed after ", max_attempts, " attempts.")
    }
    
    cat("Warning caught: increasing 'upper' to", u, "\n")
    u <- u * 1.5  # Increase upper limit by 50%
    attempt <- attempt + 1
  }
  
  cat("EB...\n")
  model_EB_trait <- mvgls(morpho~trait, data = dat, tree = tree, model = "EB")
  
  cat("BMM...\n")
  model_BMM1 <- mvgls(morpho~1, data = dat, tree = scm_consensus, model = "BMM")
  
  cat("OUM...\n")
  u <- as.numeric(model_lambda1$start_values[2]) # Start 'upper' from a reasonable value
  attempt <- 1
  repeat {
    model_OUM1 <- tryCatch(
      mvgls(morpho ~ 1, data = dat, tree = scm_consensus, model = "OUM", upper = u),
      warning = function(w) "warning"
    )
    
    if (!is.character(model_OUM1)) break  # If no warning, model fit is successful
    
    if (attempt >= max_attempts) {
      stop("Model failed after ", max_attempts, " attempts.")
    }
    
    cat("Warning caught: increasing 'upper' to", u, "\n")
    u <- u * 1.5  # Increase upper limit by 50%
    attempt <- attempt + 1
  }
  
  morpho_mvgls_models[[paste0("model_BM_", trait)]] <- model_BM_trait
  morpho_mvgls_models[[paste0("model_BMM_", trait)]] <- model_BMM1
  morpho_mvgls_models[[paste0("model_lambda_", trait)]] <- model_lambda_trait
  morpho_mvgls_models[[paste0("model_OU_", trait)]] <- model_OU_trait
  morpho_mvgls_models[[paste0("model_OUM_", trait)]] <- model_OUM1
  morpho_mvgls_models[[paste0("model_EB_", trait)]] <- model_EB_trait
  
  res <- c(res, lapply(morpho_mvgls_models[!names(morpho_mvgls_models) %in% names(res)], EIC, nbcores = ncores))
  
  # clean model objects
  rm(model_BMM1, model_OUM1, model_BM_trait, model_lambda_trait, model_OU_trait, model_EB_trait)
  } else {
    warning(trait, "rda file not found!! Skipping analysis...\n")
  }
} 

saveRDS(morpho_mvgls_models, file = paste0("./rdata/morphoPCs_mvgls.models.rds"))

model_names <- names(morpho_mvgls_models)
param_list <- lapply(model_names, function(nm) .mv_param(morpho_mvgls_models[[nm]], nm))

mvgls.morpho.restab <- data.frame(model = model_names,
                               param = sapply(param_list, function(i) i[[1]]),
                               value = sapply(param_list, function(i) i[[2]]),
                               predictor = ifelse(
                                 str_detect(word(model_names, -1, sep = "_"), "DietEcotype"),
                                 paste0("DietEcotype"), 
                                 ifelse(str_detect(word(model_names, -1, sep = "_"), "Symbiosis"),
                                        paste0("Symbiosis"), 
                                        ifelse(str_detect(word(model_names, -1, sep = "_"), "Habitat"),
                                               paste0("Habitat"),"1"))),
                               logLik = sapply(res, function(i) i$LogLikelihood),
                               k = sapply(res, function(i) i$p),
                               se = sapply(res, function(i) i$se),
                               EIC = eic <- sapply(res, function(i) i$EIC),
                               delta.EIC = eic-min(eic),
                               weight = unclass(aic.w(eic)))

mvgls.morpho.restab$predictor[grepl("BMM|OUM", mvgls.morpho.restab$model)] <- "1"

best_model_name <- mvgls.morpho.restab$model[which.min(mvgls.morpho.restab$delta.EIC)]
best_model.morpho <- morpho_mvgls_models[[best_model_name]]

# Run MANOVA test
for (trait in discrete_traits) {
  best_trait_model <- mvgls.morpho.restab[grepl(trait, mvgls.morpho.restab$predictor),]
  manova_test <- sapply(best_trait_model$model, function(m) manova.gls(morpho_mvgls_models[[m]], 
                                                                       nperm = 999, 
                                                                       test = "Pillai", 
                                                                       verbose = TRUE), simplify = FALSE)
  
  mvgls.morpho.restab$Pillai[grepl(trait, mvgls.morpho.restab$predictor)] <- sapply(manova_test, function(i) i$stat)
  mvgls.morpho.restab$pvalue[grepl(trait, mvgls.morpho.restab$predictor)] <- sapply(manova_test, function(i) i$pvalue)
}


# write table to hard drive
cat("Writting discrete morpho models result table on hard disk...\n")
write.table(mvgls.morpho.restab, file = file.path(res_folder, paste0("multivariate/morphoPCs_mvgls.models.fit.restable.txt")), 
            quote = FALSE, sep = "\t", row.names = T, col.names = T)


# =====================
# COLOUR + MORPHO
# =====================
colPCs <- colorPCs[tree$tip.label,1:8]
colnames(colPCs) <- paste0("col", colnames(colPCs))
mrphPCs <- morphoPCs[tree$tip.label,1:8]
colnames(mrphPCs) <- paste0("morpho", colnames(mrphPCs))

dat <- list(colmorpho = as.matrix(cbind(colPCs, mrphPCs)))

model_BM1 <- mvgls(colmorpho ~ 1, data = dat, tree = tree, model = "BM")
model_lambda1 <- mvgls(colmorpho ~ 1, data = dat, tree = tree, model = "lambda")
u <- as.numeric(model_lambda1$start_values[2])  # Starting value for 'upper'
max_attempts <- 10
attempt <- 1
repeat {
  model_OU1 <- tryCatch(
    mvgls(colmorpho ~ 1, data = dat, tree = tree, model = "OU", upper = u),
    warning = function(w) "warning"
  )
  
  if (!is.character(model_OU1)) break  # Model fit successful, exit loop
  
  if (attempt >= max_attempts) {
    stop("Model failed after ", max_attempts, " attempts.")
  }
  
  cat("Warning caught: increasing 'upper' to", u, "\n")
  u <- u * 1.5  # Increase 'upper' by 50%
  attempt <- attempt + 1
}
model_EB1 <- mvgls(colmorpho~1, data = dat, tree = tree, model = "EB")

colmorpho_mvgls_models <- list(
  BM1 = model_BM1, 
  lambda1 = model_lambda1,
  OU1 = model_OU1, 
  EB1 = model_EB1
)

res <- lapply(colmorpho_mvgls_models, EIC, nbcores = ncores)

# ============================
## Color vs discrete traits
# ============================
for (trait in discrete_traits){
  if (file.exists(paste0("./rdata/ASR_", trait, ".rda"))){
    ## Distance matrices
    load(paste0("./rdata/ASR_", trait, ".rda"))
    
    trait.vec = setNames(eco_traits_subset[, trait], rownames(eco_traits_subset))[tree$tip.label]
    
    cat("Working on trait:", trait, "...\n")
    tree <- fit.marginal$phy
    tree$node.label <- levels(trait.vec)[tree$node.label]
    n <- Ntip(tree)
    
    scm_consensus <- make.consensus.simmap(trait_SCM_trees)
    
    cat("Running Multivariate Phylogenetic Generalized Least Squares on color PCA, testing ", trait, "...\n")
    
    ## MVPGLS
    dat <- list(colmorpho = as.matrix(cbind(colPCs, mrphPCs)), trait = trait.vec)
    
    cat("BM...\n")
    model_BM_trait <- mvgls(colmorpho~trait, data = dat, tree = tree, model = "BM")
    
    cat("Lambda...\n")
    model_lambda_trait <- mvgls(colmorpho~trait, data = dat, tree = tree, model = "lambda")
    
    cat("OU...\n")
    u <- as.numeric(model_lambda1$start_values[2]) # Start 'upper' from a reasonable value
    attempt <- 1
    repeat {
      model_OU_trait <- tryCatch(
        mvgls(colmorpho ~ trait, data = dat, tree = tree, model = "OU", upper = u),
        warning = function(w) "warning"
      )
      
      if (!is.character(model_OU_trait)) break  # If no warning, model fit is successful
      
      if (attempt >= max_attempts) {
        stop("Model failed after ", max_attempts, " attempts.")
      }
      
      cat("Warning caught: increasing 'upper' to", u, "\n")
      u <- u * 1.5  # Increase upper limit by 50%
      attempt <- attempt + 1
    }
    
    cat("EB...\n")
    model_EB_trait <- mvgls(colmorpho~trait, data = dat, tree = tree, model = "EB")
    
    cat("BMM...\n")
    model_BMM1 <- mvgls(colmorpho~1, data = dat, tree = scm_consensus, model = "BMM")
    
    cat("OUM...\n")
    u <- as.numeric(model_lambda1$start_values[2]) # Start 'upper' from a reasonable value
    attempt <- 1
    repeat {
      model_OUM1 <- tryCatch(
        mvgls(colmorpho ~ 1, data = dat, tree = scm_consensus, model = "OUM", upper = u),
        warning = function(w) "warning"
      )
      
      if (!is.character(model_OUM1)) break  # If no warning, model fit is successful
      
      if (attempt >= max_attempts) {
        stop("Model failed after ", max_attempts, " attempts.")
      }
      
      cat("Warning caught: increasing 'upper' to", u, "\n")
      u <- u * 1.5  # Increase upper limit by 50%
      attempt <- attempt + 1
    }
    
    colmorpho_mvgls_models[[paste0("model_BM_", trait)]] <- model_BM_trait
    colmorpho_mvgls_models[[paste0("model_BMM_", trait)]] <- model_BMM1
    colmorpho_mvgls_models[[paste0("model_lambda_", trait)]] <- model_lambda_trait
    colmorpho_mvgls_models[[paste0("model_OU_", trait)]] <- model_OU_trait
    colmorpho_mvgls_models[[paste0("model_OUM_", trait)]] <- model_OUM1
    colmorpho_mvgls_models[[paste0("model_EB_", trait)]] <- model_EB_trait
    
    res <- c(res, lapply(colmorpho_mvgls_models[!names(colmorpho_mvgls_models) %in% names(res)], EIC, nbcores = ncores))
    
    # clean model objects
    rm(model_BMM1, model_OUM1, model_BM_trait, model_lambda_trait, model_OU_trait, model_EB_trait)
    
  } else {
    warning(trait, "rda file not found!! Skipping analysis...\n")
  }
}
saveRDS(colmorpho_mvgls_models, file = paste0("./rdata/colmorphoPCs_mvgls.models.rds"))

model_names <- names(colmorpho_mvgls_models)
param_list <- lapply(model_names, function(nm) .mv_param(colmorpho_mvgls_models[[nm]], nm))

mvgls.colmorpho.restab <- data.frame(model = model_names,
                               param = sapply(param_list, function(i) i[[1]]),
                               value = sapply(param_list, function(i) i[[2]]),
                               predictor = ifelse(
                                 str_detect(word(model_names, -1, sep = "_"), "DietEcotype"),
                                 paste0("DietEcotype"), 
                                 ifelse(str_detect(word(model_names, -1, sep = "_"), "Symbiosis"),
                                        paste0("Symbiosis"), 
                                        ifelse(str_detect(word(model_names, -1, sep = "_"), "Habitat"),
                                               paste0("Habitat"),"1"))),
                               logLik = sapply(res, function(i) i$LogLikelihood),
                               k = sapply(res, function(i) i$p),
                               se = sapply(res, function(i) i$se),
                               EIC = eic <- sapply(res, function(i) i$EIC),
                               delta.EIC = eic-min(eic),
                               weight = unclass(aic.w(eic)))

mvgls.colmorpho.restab$predictor[grepl("BMM|OUM", mvgls.colmorpho.restab$model)] <- "1"

best_model_name <- mvgls.colmorpho.restab$model[which.min(mvgls.colmorpho.restab$delta.EIC)]
best_model.col <- colmorpho_mvgls_models[[best_model_name]]

# Run MANOVA test
for (trait in discrete_traits) {
  best_trait_model <- mvgls.colmorpho.restab[grepl(trait, mvgls.colmorpho.restab$predictor),]
  manova_test <- sapply(best_trait_model$model, function(m) manova.gls(colmorpho_mvgls_models[[m]], 
                                                                       nperm = 999, 
                                                                       test = "Pillai", 
                                                                       verbose = TRUE), simplify = FALSE)
  mvgls.colmorpho.restab$Pillai[grepl(trait, mvgls.colmorpho.restab$predictor)] <- sapply(manova_test, function(i) i$stat)
  mvgls.colmorpho.restab$pvalue[grepl(trait, mvgls.colmorpho.restab$predictor)] <- sapply(manova_test, function(i) i$pvalue)
}

# write table to hard drive
cat("Writting discrete color models result table on hard disk...\n")
write.table(mvgls.colmorpho.restab, file = file.path(res_folder, paste0("multivariate/color-morphoPCs_mvgls.models.fit.restable.txt")), 
            quote = FALSE, sep = "\t", row.names = T, col.names = T)

# remove all models
rm(list = ls(pattern = "^model_"))

# ======================= #
# === OUTPUT TABLES
# ======================= #
col_mvgls_models <- readRDS(paste0("./rdata/colPCs_mvgls.models.rds"))
morpho_mvgls_models <- readRDS(paste0("./rdata/morphoPCs_mvgls.models.rds"))
colmorpho_mvgls_models <- readRDS(paste0("./rdata/colmorphoPCs_mvgls.models.rds"))

## read tables
mvgls_col.tab <- read.table(file = file.path(res_folder, paste0("multivariate/colorPCs_mvgls.models.fit.restable.txt")))
mvgls_morpho.tab <- read.table(file = file.path(res_folder, paste0("multivariate/morphoPCs_mvgls.models.fit.restable.txt")))
mvgls_colmorpho.tab <- read.table(file = file.path(res_folder, paste0("multivariate/color-morphoPCs_mvgls.models.fit.restable.txt")))

## Turn a data frame into a LaTeX table
make_latex_from_tab <- function(tab,
                                caption = "Model comparison",
                                label   = "tab:model_comparison",
                                include_param = TRUE,
                                support_delta = 2, 
                                support_weight = 0.1,
                                digits = 3) {
  rename_trait <- function(x) {
    x <- sub("DietEcotype", "Diet", x)
    x <- sub("Symbiosis", "Symbiosis", x)
    x <- sub("Habitat", "Habitat", x)
    x
  }
  stopifnot(all(c("model","param","value","logLik","k","se","EIC","delta.EIC","weight") %in% colnames(tab)))
  
  # Best model (null or ecological)
  best_overall_idx <- which.min(tab$EIC)
  best_overall     <- tab[best_overall_idx, ]
  
  is_trait <- grepl("DietEcotype|Symbiosis|Habitat", best_overall$model)
  if (any(is_trait)) {
    best_trait     <- tab[best_overall_idx, ]
    trait_supported <- (best_trait$delta.EIC <= support_delta) || (best_trait$weight >= support_weight)
    pred <- "trait"
  } else {
    best_trait <- NULL
    trait_supported <- FALSE
    pred <- "NULL"
  }
  
  
  # Predictor (~1 vs ~trait) and model family (BM, BMM, lambda, OU, OUM, EB)
  out <- tab %>%
    mutate(
      Predictor = ifelse(
        str_detect(word(model, -1, sep = "_"), "DietEcotype"),
        paste0("$\\sim$Diet"), ifelse(str_detect(word(model, -1, sep = "_"), "Symbiosis"),
        paste0("$\\sim$Symbiosis"), ifelse(str_detect(word(model, -1, sep = "_"), "Habitat"),
        paste0("$\\sim$Habitat"),
        "$\\sim$1"
      ))),
      Model = case_when(
        grepl("BMM", model)     ~ sub("_", "|", sub("model_", "", rename_trait(model))),
        grepl("BM", model)      ~ "BM",
        grepl("lambda", model)  ~ "Pagel's $\\lambda$",
        grepl("OUM", model) ~ sub("_", "|", sub("model_", "", rename_trait(model))),
        grepl("OU", model)     ~ "OU",
        grepl("EB", model)      ~ "EB",
        TRUE ~ model
      ),
      param_disp = case_when(
        param == "lambda" ~ "$\\lambda$",
        param == "alpha"  ~ "$\\alpha$",
        param == "a"      ~ "$a$",
        TRUE              ~ "---"
      ),
      value_disp = ifelse(is.na(value), "---", sprintf(paste0("%.", digits, "f"), value)),
      logLik = sprintf(paste0("%.", digits, "f"), logLik),
      se     = sprintf(paste0("%.", digits, "f"), se),
      EIC    = sprintf(paste0("%.", digits, "f"), EIC),
      DeltaEIC = sprintf(paste0("%.", digits, "f"), delta.EIC),
      # (re)compute weights for safety from ΔEIC
      weight_num = exp(-0.5 * as.numeric(delta.EIC)),
      weight_num = weight_num / sum(weight_num, na.rm = TRUE),
      weight = sprintf(paste0("%.", digits, "f"), weight_num),
      Pillai = ifelse(is.na(Pillai), "---", sprintf(paste0("%.", digits, "f"), Pillai)),
      pvalue = ifelse(is.na(pvalue), "---", fmt_p(pvalue))
    )
  out$Predictor[grepl("BMM|OUM", out$model)] <- "$\\sim$1"
  
  # Order rows to mirror your examples
  out <- out %>%
    mutate(
      Model     = factor(Model, levels = c("BM","Pagel's $\\lambda$","OU","EB","BMM|Diet","BMM|Habitat","BMM|Symbiosis","OUM|Diet", "OUM|Habitat","OUM|Symbiosis")),
      Predictor = factor(Predictor, levels = c("$\\sim$1","$\\sim$Diet", "$\\sim$Habitat", "$\\sim$Symbiosis"))
    ) %>%
    arrange(Predictor, Model)
  
  # Choose columns
  if (include_param) {
    disp <- out %>%
      dplyr::select(Predictor, Model, value_disp,
             logLik, k, se, EIC, DeltaEIC, weight, Pillai, pvalue)
    colnames(disp) <- c("Predictor","Model", "$\\lambda$/$\\alpha$/$a$",
                        "logLik","k","se","EIC","$\\Delta$EIC","weight", "Pillai", "\\textit{p}-value")
  } else {
    disp <- out %>%
      dplyr::select(Predictor, Model, logLik, k, se, EIC, `ΔEIC`, weight, Pillai, pvalue)
    colnames(disp) <- c("Predictor","Model","logLik","k","se","EIC","$\\Delta$EIC","weight", "Pillai", "\\textit{p}-value")
  }
  
  # Build the LaTeX table
  kb <- kable(disp, format = "latex", booktabs = TRUE, escape = FALSE, row.names = FALSE,
              caption = caption, label = label, align = "llccccccccc")
  
  # Bold the best row (global min ΔEIC)
  best_rows <- which(out$delta.EIC < support_delta | out$weight > support_weight)
  kb <- row_spec(kb, best_rows, bold = TRUE) 
  
  kb %>% kable_styling(full_width = FALSE)
}

# Color
# --------
make_latex_from_tab(
  mvgls_col.tab,
  caption = "",
  label   = "mvgls_col",
  include_param = TRUE
)

# Morpho
# ---------
make_latex_from_tab(
  mvgls_morpho.tab,
  caption = "",
  label   = "mvgls_morpho",
  include_param = TRUE
)

# Color-Morpho
# ---------
make_latex_from_tab(
  mvgls_colmorpho.tab,
  caption = "",
  label   = "mvgls_col_morpho",
  include_param = TRUE
)

best_model_name <- mvgls_colmorpho.tab %>% filter(delta.EIC == 0) %>%
  pull(model)

fit <- colmorpho_mvgls_models[[best_model_name]]

S <- fit$sigma$S

## Long format for plotting ---
df_cov <- as.data.frame(as.table(S))
names(df_cov) <- c("Trait1", "Trait2", "Cov")

df <- df_cov %>%
  mutate(
    Type1  = ifelse(grepl("^colPC", Trait1), "Colour", "Morphology"),
    Type2  = ifelse(grepl("^colPC", Trait2), "Colour", "Morphology")
  ) %>%
  filter(Type1 == "Colour", Type2 == "Morphology")

# Symmetric limits for diverging colour scale
zlim <- max(abs(df$Cov), na.rm = TRUE)

## Covariance heatmap 
## ===================
gg_cov <- ggplot(df, aes(x = Trait1, y = Trait2)) +
  geom_tile(aes(fill = Cov), colour = "grey90", size = 0.2) +
  scale_fill_gradient2(
    limits = c(-zlim, zlim),
    low    = "#D73027", mid = "white", high = "#4575B4",
    name   = "Covariance"
  ) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title  = element_blank(),
    panel.grid  = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA)
  )

gg_cov

ggsave(sprintf("figures/col_morpho/col_morpho_mvgls_%s_covariances_heatmap.pdf", best_model_name),
       gg_cov, width = 7, height = 6)

## Table of covariances 
sig_pairs_cross <- df %>%
  arrange(desc(abs(Cov))) %>%
  dplyr::select(Trait1, Trait2, Cov)

write.csv(sig_pairs_cross, 
          sprintf("results/col_morpho/col_morpho_mvgls_%s_covariances_tab.csv", best_model_name),
          row.names = FALSE)

sig_pairs_cross %>%
  filter(abs(Cov) > quantile(Cov, 0.95))

#=============================#
# Plots
#=============================#
col_traitList_best_models <- list()
morpho_traitList_best_models <- list()

palette <- c("Zissou1", "Darjeeling1", "Darjeeling2", "FantasticFox1")[4]
for (trait in discrete_traits){
  
  trait.vec = setNames(eco_traits_subset[, trait], rownames(eco_traits_subset))[tree$tip.label]
  
  if (file.exists(paste0("./rdata/colorPCs.vs.", trait, "_mvgls.models.rds"))) {
    ## load models
    col_mvgls_models <- readRDS(file.path("rdata", paste0("colorPCs.vs.", trait, "_mvgls.models.rds")))
    trait_morph.models <- readRDS(file.path("rdata", paste0("morphoPCs.vs.", trait, "_mvgls.models.rds")))
    
    ## load res tables
    col_restable <- read.table( file.path(res_folder, paste0("colorPCs.vs.", trait, "_mvgls.models.fit.restable.txt")))
    morph_restable <- read.table( file.path(res_folder, paste0("morphoPCs.vs.", trait, "_mvgls.models.fit.restable.txt")))
    
    ## get best models
    col_best_model_name <- col_restable$model[which.min(col_restable$EIC)]
    morpho_best_model_name <- morph_restable$model[which.min(morph_restable$EIC)]
    
    best_model.col <- col_mvgls_models[[col_best_model_name]]
    best_model.morph <- trait_morph.models[[morpho_best_model_name]]
   
    if (nrow(best_model.col$coefficients) > 1) { rownames(best_model.col$coefficients) <- levels(trait.vec) }
    if (nrow(best_model.morph$coefficients) > 1) { rownames(best_model.morph$coefficients) <- levels(trait.vec) }
    
    col_traitList_best_models[[trait]] <- best_model.col
    morpho_traitList_best_models[[trait]] <- best_model.morph
    
    col_formula <- as.character(best_model.col$call)[2]
    morpho_formula <- as.character(best_model.morph$call)[2]
    
    col_param <- round(as.numeric(best_model.col$param), 2)
    morpho_param <- round(as.numeric(best_model.morph$param), 2)
    
    col_model_name <- sub("_trait|1", "", col_best_model_name)
    morpho_model_name <- sub("_trait|1", "", morpho_best_model_name)
    
    col_param_symbol <- if (grepl("lambda", col_best_model_name, ignore.case = TRUE)) {
      expression(lambda)
    } else if (grepl("BM", col_best_model_name, ignore.case = TRUE)) {
      expression(gamma)
    } else if (grepl("OU", col_best_model_name, ignore.case = TRUE)) {
      expression(alpha)
    } else if (grepl("EB", col_best_model_name, ignore.case = TRUE)) {
      expression(r)
    }
    
    morpho_param_symbol <- if (grepl("lambda", morpho_best_model_name, ignore.case = TRUE)) {
      expression(lambda)
    } else if (grepl("BM", morpho_best_model_name, ignore.case = TRUE)) {
      expression(gamma)
    } else if (grepl("OU", morpho_best_model_name, ignore.case = TRUE)) {
      expression(alpha)
    } else if (grepl("EB", morpho_best_model_name, ignore.case = TRUE)) {
      expression(r)
    }
    
    col_label <- bquote(atop(paste("Best model:" ~ .(col_model_name) ~ "(" * .(col_formula) * ");" ~ lambda * ":" ~ .(col_param))))
    morpho_label <- bquote(atop(paste("Best model:" ~ .(morpho_model_name) ~ "(" * .(morpho_formula) * ");" ~ lambda * ":" ~ .(morpho_param))))
    
    #trait.cols <- setNames(viridis::viridis(length(levels(trait.vec))), levels(trait.vec))
    trait.cols <- setNames(as.character(wes_palette(palette, length(levels(trait.vec)))),  levels(trait.vec))
    
    pdf(paste0("./figures/multivariate/col-morph.vs.", trait, ".pdf"), height = 15, width = 10)
    par(mar = c(5.1, 6.1, 2.1, 2.1))
    layout(matrix(c(1:9,9), nrow = 5, ncol = 2, byrow = TRUE), heights = c(1,1,1,1,0.5))
    for (i in c(1,3,5,7)){
      pc1 <- i
      pc2 <- i + 1
    
      ## color patterns
      pervarPCs <- round(imgPCA_1183damsel$pca$sdev^2 / sum(imgPCA_1183damsel$pca$sdev^2) * 100,2)
      xlim = range(c(best_model.col$coefficients[,pc1], colorPCs[,pc1]))
      expansion <- diff(xlim) * 0.1
      xlim <- c(xlim[1] - expansion, xlim[2] + expansion)
      
      ylim = range(c(best_model.col$coefficients[,pc2], colorPCs[,pc2]))
      expansion <- diff(ylim) * 0.1
      ylim <- c(ylim[1] - expansion, ylim[2] + expansion)
                
      plot(colorPCs[,c(pc1,pc2)], col = trait.cols[trait.vec[rownames(colorPCs)]], 
           xlim = xlim, ylim = ylim, cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.7,
           type = "n", xlab = paste0("PC", pc1, " (",pervarPCs[pc1], "%)"), ylab = paste0("PC2", pc2, " (",pervarPCs[pc2], "%)"), 
           main = ifelse(i == 1, "Color patterns", ""))
      for(sp in rownames(colorPCs)){
        if (nrow(best_model.col$coefficients) == 1) {
          rw <- 1
        } else {
          rw <- trait.vec[sp]
        }
        lines(c(colorPCs[sp,pc1], best_model.col$coefficients[rw,pc1]), 
              c(colorPCs[sp,pc2], best_model.col$coefficients[rw,pc2]),
              col = adjustcolor(trait.cols[trait.vec[sp]], .5), lwd = 2)
      }
      p.col <- if (nrow(best_model.col$coefficients) == 1) "grey95" else colorspace::darken(trait.cols[rownames(best_model.col$coefficients)], 0.1)
      points(best_model.col$coefficients[,pc1], best_model.col$coefficients[,pc2], 
             col = colorspace::darken(p.col, 0.25), bg = p.col, pch = 21, cex = 2)
      if (i == 1) {
        usr <- par("usr") 
        text(x = usr[2] - diff(usr[1:2]) * 0.05,  y = usr[3] + diff(usr[3:4]) * 0.05,
             label =  col_label,           
             adj = c(1,1), cex = 0.9, col = "grey40")
      }
      
      ## Morphology
      pervarPCs <- round(morphoPCA$sdev^2 / sum(morphoPCA$sdev^2) * 100, 2)
      xlim = range(c(best_model.morph$coefficients[,pc1], morphoPCs[,pc1]))
      expansion <- diff(xlim) * 0.1
      xlim <- c(xlim[1] - expansion, xlim[2] + expansion)
      
      ylim = range(c(best_model.morph$coefficients[,pc2], morphoPCs[,pc2]))
      expansion <- diff(ylim) * 0.1
      ylim <- c(ylim[1] - expansion, ylim[2] + expansion)
      
      plot(morphoPCs[,c(pc1,pc2)], col = trait.cols[trait.vec[rownames(morphoPCs)]], 
           xlim = xlim, ylim = ylim, cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.7,
           type = "n", xlab = paste0("PC", pc1, " (",pervarPCs[pc1], "%)"), ylab = paste0("PC", pc2, " (",pervarPCs[pc2], "%)"), 
           main = ifelse(i==1, "Morphology", ""))
      for(sp in rownames(morphoPCs)){
        if (nrow(best_model.morph$coefficients) == 1) {
          rw <- 1
        } else {
          rw <- trait.vec[sp]
        }
        lines(c(morphoPCs[sp,pc1], best_model.morph$coefficients[rw,pc1]), 
              c(morphoPCs[sp,pc2], best_model.morph$coefficients[rw,pc2]),
              col = adjustcolor(trait.cols[trait.vec[sp]], .25), lwd = 2)
      }
      p.col <- if (nrow(best_model.morph$coefficients) == 1) "grey95" else colorspace::darken(trait.cols[rownames(best_model.morph$coefficients)], 0.1)
      points(best_model.morph$coefficients[,pc1], best_model.morph$coefficients[,pc2], 
             col = colorspace::darken(p.col, 0.25), bg = p.col, pch = 21, cex = 2)
      if (i == 1) {
        usr <- par("usr") 
        text(x = usr[2] - diff(usr[1:2]) * 0.05, y = usr[3] + diff(usr[3:4]) * 0.05,
             label = morpho_label,
             adj = c(1,1), cex = 0.9, col = "grey40")
      }
    }
    legend_text <- setNames(levels(trait.vec), levels(trait.vec))
    plot.new()
    par(mar = c(0,1,1,1))
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    legend("top", legend = legend_text, 
           fill = trait.cols, cex = 2, ncol = ifelse(length(trait.vec) > 3, 3, length(levels(trait.vec))),
           box.lwd = 0, bg = "transparent")
    
    dev.off()
    
    ## plot PCA + density plots
    dens <- lapply(c("x", "y"), function(co){
      sapply(sprintf("PC%s", 1:8), function(pc){
        a <- lapply(levels(trait.vec), function(an){
          if (length(colorPCs[rownames(colorPCs) %in% names(trait.vec)[trait.vec == an],1]) > 1){
            density(colorPCs[rownames(colorPCs) %in% names(trait.vec)[trait.vec == an],1])[[co]]
          } else {
            NA
          }
          
        } )
        names(a) <- levels(trait.vec)
        return(a)
      })
    })
    
    pdf(paste0("./figures/multivariate/ColorPCA1-4.vs.", trait, ".pdf"), width = 20, height = 12)
    cex.pt = 1.1
    alpha.dens <- .5
    par(fig= c(0,0.4,0,0.75), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 5, 1, 1), mgp = c(2, 0.5, 0), tck = -0.01)
    downmar <- upmar <- leftmar <- par()$mar
    downmar[3:4] <- 0
    par(mar = downmar)
    pervarPCs <- round(imgPCA_1183damsel$pca$sdev^2 / sum(imgPCA_1183damsel$pca$sdev^2) * 100,2)
    plot(colorPCs[,c(1,2)], cex = cex.pt, col = trait.cols[trait.vec[rownames(colorPCs)]], 
         xlim = max(abs(colorPCs[,1])) * c(-1,1), ylim = max(abs(colorPCs[,2])) * c(-1,1), cex.axis = 1,
         pch = 16, xlab = paste0("PC1 (",pervarPCs[1], "%)"), ylab = paste0("PC2 (",pervarPCs[2], "%)"))
    lims <- par()$usr
    upmar[c(1,4)] <- 0
    par(fig = c(0.0, 0.4, 0.75, 1), new = T, mar = upmar)
    plot(1, type="n", xlab="", ylab="", xlim=lims[1:2], 
         ylim=c(0, max(do.call(c,dens[[2]][,1]), na.rm = T)), xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[1]][an,1][[1]], dens[[2]][an,1][[1]], xlim=lims[1:2],
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    legend("topleft", legend = legend_text, bty = "n",
           pch = 16, col = trait.cols[levels(trait.vec)], cex = 1.2)
    leftmar[c(2,3)] <- 0
    #leftmar[1] <- 9.25
    par(fig = c(0.4, 0.5, 0, 0.75), new = T, mar = leftmar)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, max(do.call(c,dens[[2]][,2]), na.rm = T)), 
         ylim=lims[3:4], xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[2]][an,2][[1]], dens[[1]][an,2][[1]], 
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    
    par(fig= c(0.5,0.9,0,0.75), new = T)
    par(mar = downmar)
    plot(colorPCs[,c(3,4)],cex = cex.pt,  col = trait.cols[trait.vec[rownames(colorPCs)]], pch = 16,
         xlim = max(abs(colorPCs[,3])) * c(-1,1), ylim = max(abs(colorPCs[,4])) * c(-1,1), cex.axis = 1,
         xlab = paste0("PC3 (",pervarPCs[3], "%)"), ylab = paste0("PC4 (",pervarPCs[4], "%)"))
    lims <- par()$usr
    #upmar[2] = 9.25
    par(fig = c(0.5, 0.9, 0.75, 1), new = T, mar = upmar)
    plot(1, type="n", xlab="", ylab="", xlim=lims[1:2], ylim=c(0, max(do.call(c,dens[[2]][,3]), na.rm = T)), 
         xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[1]][an,3][[1]], dens[[2]][an,3][[1]], 
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    leftmar[c(2,3)] <- 0
    #leftmar[1] <- 13.25
    par(fig = c(0.9, 1, 0.0, 0.75), new = T, mar = leftmar)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, max(do.call(c,dens[[2]][,4]), na.rm = T)), ylim=lims[3:4], xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[2]][an,4][[1]], dens[[1]][an,4][[1]], col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    dev.off()
    
    pdf(paste0("./figures/multivariate/ColorPCA5-8.vs.", trait, ".pdf"), width = 20, height = 12)
    cex.pt = 1.1
    alpha.dens <- .5
    par(fig= c(0,0.4,0,0.75), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 5, 1, 1), mgp = c(2, 0.5, 0), tck = -0.01)
    downmar <- upmar <- leftmar <- par()$mar
    downmar[3:4] <- 0
    par(mar = downmar)
    pervarPCs <- round(imgPCA_1183damsel$pca$sdev^2 / sum(imgPCA_1183damsel$pca$sdev^2) * 100,2)
    plot(colorPCs[,c(5,6)], cex = cex.pt, col = trait.cols[trait.vec[rownames(colorPCs)]], 
         xlim = max(abs(colorPCs[,5])) * c(-1,1), ylim = max(abs(colorPCs[,6])) * c(-1,1), cex.axis = 1,
         pch = 16, xlab = paste0("PC5 (",pervarPCs[5], "%)"), ylab = paste0("PC6 (",pervarPCs[6], "%)"))
    lims <- par()$usr
    upmar[c(1,4)] <- 0
    #upmar[2] <- 11.75
    par(fig = c(0.0, 0.4, 0.75, 1), new = T, mar = upmar)
    plot(1, type="n", xlab="", ylab="", xlim=lims[1:2], 
         ylim=c(0, max(do.call(c,dens[[2]][,5]), na.rm = T)), xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[1]][an,5][[1]], dens[[2]][an,5][[1]], xlim=lims[1:2],
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    legend("topleft", legend = legend_text, bty = "n",
           pch = 16, col = trait.cols[levels(trait.vec)], cex = 1.2)
    leftmar[c(2,3)] <- 0
    #leftmar[1] <- 11.25
    par(fig = c(0.4, 0.5, 0.0, 0.75), new = T, mar = leftmar)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, max(do.call(c,dens[[2]][,6]), na.rm = T)), 
         ylim=lims[3:4], xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[2]][an,6][[1]], dens[[1]][an,6][[1]], 
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    
    par(fig= c(0.5,0.9,0,0.75), new = T)
    par(mar = downmar)
    plot(colorPCs[,c(7,8)],cex = cex.pt,  col = trait.cols[trait.vec[rownames(colorPCs)]], pch = 16,
         xlim = max(abs(colorPCs[,7])) * c(-1,1), ylim = max(abs(colorPCs[,8])) * c(-1,1), cex.axis = 1,
         xlab = paste0("PC7 (",pervarPCs[7], "%)"), ylab = paste0("PC8 (",pervarPCs[8], "%)"))
    lims <- par()$usr
    #upmar[2] <- 13.25
    par(fig = c(0.5, 0.9, 0.75, 1), new = T, mar = upmar)
    plot(1, type="n", xlab="", ylab="", xlim=lims[1:2], ylim=c(0, max(do.call(c,dens[[2]][,7]), na.rm = T)), 
         xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[1]][an,7][[1]], dens[[2]][an,7][[1]], 
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    leftmar[c(2,3)] <- 0
    #leftmar[1] <- 15
    par(fig = c(0.9, 1, 0.0, 0.75), new = T, mar = leftmar)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, max(do.call(c,dens[[2]][,8]), na.rm = T)), ylim=lims[3:4], xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[2]][an,8][[1]], dens[[1]][an,8][[1]], col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    dev.off()
    
    dens <- lapply(c("x", "y"), function(co){
      sapply(sprintf("PC%s", 1:ncol(morphoPCs)), function(pc){
        a <- lapply(levels(trait.vec), function(an){
          if (length(morphoPCs[rownames(morphoPCs) %in% names(trait.vec)[trait.vec == an],1]) > 1){
            density(morphoPCs[rownames(morphoPCs) %in% names(trait.vec)[trait.vec == an],1])[[co]]
          } else {
            NA
          }
          
        } )
        names(a) <- levels(trait.vec)
        return(a)
      })
    })
    
    pdf(paste0("./figures/multivariate/MorphPCA1-4.vs.", trait, ".pdf"), width = 20, height = 12)
    cex.pt = 1.1
    alpha.dens <- .5
    par(fig= c(0,0.4,0,0.75), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 5, 1, 1), mgp = c(2, 0.5, 0), tck = -0.01)
    downmar <- upmar <- leftmar <- par()$mar
    downmar[3:4] <- 0
    par(mar = downmar)
    pervarPCs <- round(morphoPCA$sdev^2 / sum(morphoPCA$sdev^2) * 100,2)
    plot(morphoPCs[,c(1,2)], cex = cex.pt, col = trait.cols[trait.vec[rownames(colorPCs)]], 
         xlim = max(abs(morphoPCs[,1])) * c(-1,1), ylim = max(abs(morphoPCs[,2])) * c(-1,1), cex.axis = 1,
         pch = 16, xlab = paste0("PC1 (",pervarPCs[1], "%)"), ylab = paste0("PC2 (",pervarPCs[2], "%)"))
    lims <- par()$usr
    upmar[c(1,4)] <- 0
    par(fig = c(0.0, 0.4, 0.75, 1), new = T, mar = upmar)
    plot(1, type="n", xlab="", ylab="", xlim=lims[1:2], 
         ylim=c(0, max(do.call(c,dens[[2]][,1]), na.rm = T)), xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[1]][an,1][[1]], dens[[2]][an,1][[1]], xlim=lims[1:2],
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    legend("topleft", legend = legend_text, bty = "n",
           pch = 16, col = trait.cols[levels(trait.vec)], cex = 1.2)
    leftmar[c(2,3)] <- 0
    #leftmar[1] <- 9.25
    par(fig = c(0.4, 0.5, 0, 0.75), new = T, mar = leftmar)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, max(do.call(c,dens[[2]][,2]), na.rm = T)), 
         ylim=lims[3:4], xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[2]][an,2][[1]], dens[[1]][an,2][[1]], 
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    
    par(fig= c(0.5,0.9,0,0.75), new = T)
    par(mar = downmar)
    plot(morphoPCs[,c(3,4)],cex = cex.pt,  col = trait.cols[trait.vec[rownames(morphoPCs)]], pch = 16,
         xlim = max(abs(morphoPCs[,3])) * c(-1,1), ylim = max(abs(morphoPCs[,4])) * c(-1,1), cex.axis = 1,
         xlab = paste0("PC3 (",pervarPCs[3], "%)"), ylab = paste0("PC4 (",pervarPCs[4], "%)"))
    lims <- par()$usr
    #upmar[2] = 9.25
    par(fig = c(0.5, 0.9, 0.75, 1), new = T, mar = upmar)
    plot(1, type="n", xlab="", ylab="", xlim=lims[1:2], ylim=c(0, max(do.call(c,dens[[2]][,3]), na.rm = T)), 
         xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[1]][an,3][[1]], dens[[2]][an,3][[1]], 
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    leftmar[c(2,3)] <- 0
    #leftmar[1] <- 13.25
    par(fig = c(0.9, 1, 0.0, 0.75), new = T, mar = leftmar)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, max(do.call(c,dens[[2]][,4]), na.rm = T)), ylim=lims[3:4], xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[2]][an,4][[1]], dens[[1]][an,4][[1]], col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    dev.off()
    
    pdf(paste0("./figures/multivariate/MorphPCA5-8.vs.", trait, ".pdf"), width = 20, height = 12)
    cex.pt = 1.1
    alpha.dens <- .5
    par(fig= c(0,0.4,0,0.75), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 5, 1, 1), mgp = c(2, 0.5, 0), tck = -0.01)
    downmar <- upmar <- leftmar <- par()$mar
    downmar[3:4] <- 0
    par(mar = downmar)
    pervarPCs <- round(imgPCA_1183damsel$pca$sdev^2 / sum(imgPCA_1183damsel$pca$sdev^2) * 100,2)
    plot(morphoPCs[,c(5,6)], cex = cex.pt, col = trait.cols[trait.vec[rownames(morphoPCs)]], 
         xlim = max(abs(morphoPCs[,5])) * c(-1,1), ylim = max(abs(morphoPCs[,6])) * c(-1,1), cex.axis = 1,
         pch = 16, xlab = paste0("PC5 (",pervarPCs[5], "%)"), ylab = paste0("PC6 (",pervarPCs[6], "%)"))
    lims <- par()$usr
    upmar[c(1,4)] <- 0
    #upmar[2] <- 11.75
    par(fig = c(0.0, 0.4, 0.75, 1), new = T, mar = upmar)
    plot(1, type="n", xlab="", ylab="", xlim=lims[1:2], 
         ylim=c(0, max(do.call(c,dens[[2]][,5]), na.rm = T)), xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[1]][an,5][[1]], dens[[2]][an,5][[1]], xlim=lims[1:2],
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    legend("topleft", legend = legend_text, bty = "n",
           pch = 16, col = trait.cols[levels(trait.vec)], cex = 1.2)
    leftmar[c(2,3)] <- 0
    #leftmar[1] <- 11.25
    par(fig = c(0.4, 0.5, 0.0, 0.75), new = T, mar = leftmar)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, max(do.call(c,dens[[2]][,6]), na.rm = T)), 
         ylim=lims[3:4], xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[2]][an,6][[1]], dens[[1]][an,6][[1]], 
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    
    par(fig= c(0.5,0.9,0,0.75), new = T)
    par(mar = downmar)
    plot(morphoPCs[,c(7,8)],cex = cex.pt,  col = trait.cols[trait.vec[rownames(morphoPCs)]], pch = 16,
         xlim = max(abs(morphoPCs[,7])) * c(-1,1), ylim = max(abs(morphoPCs[,8])) * c(-1,1), cex.axis = 1,
         xlab = paste0("PC7 (",pervarPCs[7], "%)"), ylab = paste0("PC8 (",pervarPCs[8], "%)"))
    lims <- par()$usr
    #upmar[2] <- 13.25
    par(fig = c(0.5, 0.9, 0.75, 1), new = T, mar = upmar)
    plot(1, type="n", xlab="", ylab="", xlim=lims[1:2], ylim=c(0, max(do.call(c,dens[[2]][,7]), na.rm = T)), 
         xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[1]][an,7][[1]], dens[[2]][an,7][[1]], 
              col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    leftmar[c(2,3)] <- 0
    #leftmar[1] <- 15
    par(fig = c(0.9, 1, 0.0, 0.75), new = T, mar = leftmar)
    plot(1, type="n", xlab="", ylab="", xlim=c(0, max(do.call(c,dens[[2]][,8]), na.rm = T)), ylim=lims[3:4], xaxt = "n", yaxt = "n", bty = "n")
    for(an in unique(trait.vec)){
      polygon(dens[[2]][an,8][[1]], dens[[1]][an,8][[1]], col = alpha(trait.cols[an], alpha.dens), density = NA)
    }
    dev.off()
  }
}

saveRDS(col_traitList_best_models, "./rdata/mv_colPCA_vstrait_bestmodels.rds")
saveRDS(morpho_traitList_best_models, "./rdata/mv_morphoPCA_vstrait_bestmodels.rds")

sapply(col_traitList_best_models, function(x) data.table::last(strsplit(as.character(x$call), '"')))
sapply(morpho_traitList_best_models, function(x) data.table::last(strsplit(as.character(x$call), '"')))

# lambdas
sapply(col_traitList_best_models[c("DietEcotype", "Symbiosis", "Habitat")], function(x) x$param)
sapply(morpho_traitList_best_models[c("DietEcotype", "Symbiosis", "Habitat")], function(x) x$param)
