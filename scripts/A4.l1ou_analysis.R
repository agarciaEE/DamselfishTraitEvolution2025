lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", 
         "lme4", "mgcv", "phytools", "corHMM",
         "lmtest", "ColorAR", "castor",
         "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/")

dir.create("new_figures/l1ou")

# load custom functions
source("./Rscripts/A0.Input_data.R")

##=============================================##
##                                             ##
##                Run analysis                 ##
##                                             ##
##=============================================##
dir.create("./new_figures/l1ou", showWarnings = FALSE)

##### Coloration patterns (PCA data)
#########################################

pcs <- c("all", as.character(1:8))
criterion = "pBIC" # try BIC too
ncores = 5
max.nShifts = 50
for (i in pcs){
  
  Y <- get(paste0("colorPC", ifelse(i=="all", "s", i)))
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
  
  pdf(paste0("./new_figures/l1ou/colorPC_", i, "_independent_shifts_", criterion, "_new.pdf"), width = 20, height = 25)
  ew <- rep(1,nEdges) # to set default edge width of 1
  ew[fit_ind$shift.configuration] <- 3 # to widen edges with a shift
  plot(fit_ind, cex=0.5, label.offset=0.02, edge.width=ew)
  dev.off()
  
  # save model
  saveRDS(fit_ind, paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_ind_shifts_", criterion, "_new.rds"))
  
  ## then detect which of these shifts are convergent:
  fit_conv <- estimate_convergent_regimes(fit_ind, criterion="pBIC", nCores = 5)
  print(configuration_ic(trait.data$tree, fit_conv$Y, fit_conv$shift.configuration, criterion=criterion))
  
  pdf(paste0("./new_figures/l1ou/colorPC_", i, "_convergent_shifts_", criterion, "_new.pdf"), width = 20, height = 25)
  ew <- rep(1,nEdges) # to set default edge width of 1
  ew[fit_conv$shift.configuration] <- 3 # to widen edges with a shift
  plot(fit_conv, cex=0.5, label.offset=0.02, edge.width=ew)
  dev.off()
  
  # save model
  saveRDS(fit_conv, paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, "_new.rds"))
  
}

# Custom plots:
###############
pcs <- 1:ncol(colorPCs) # first plot all, then select those with convergent regimes
pcs <- c(1,2,6)
convergent_only = FALSE
criterion = "pBIC" 
pdf(paste0("./new_figures/l1ou/l1ou_color_", ifelse(convergent_only, "conv", "all"), "_pBic_PC1-8.pdf"), width = 4 * length(pcs), height = 10)

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
for (i in pcs) {
  fit_conv <- readRDS(paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds"))
  
  if (i == 1) { 
    par(mar = c(2, 1, 4, 2))  # smaller left margin for first PC plot
  } else {
    par(mar = c(2, 4, 4, 2))  # regular margins for the rest
  }
  
  plot_shift_optima(fit_conv, 
                    color.by = "regime", 
                    convergent.only = convergent_only,
                    shift.edge.width = 5, 
                    cex.value = 1.1, 
                    cex = 0.5, 
                    show.tip.label = FALSE, 
                    title = paste0("color PC", i), 
                    legend.position = "topleft",
                    cex.legend = 1)
}

dev.off()

# test differential trait evolutionary rate across l1ou regimes
###############################################################
pcs <- 1:8
pcs <- c(1,2,6)
convergent_only = FALSE
init_iter <- 2000
max_iter <- 10000
num_cores <- 5
n_simulations <- 100
color_rates <- list()
for (i in pcs) {
  
  maxit <- init_iter
  
  cat("Processing color PC", i, "...\n")
  
  fit <- readRDS(paste0("./rdata/colPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds"))
  
  mapped_tree <- l1ou2simmap(fit, convergent.only = TRUE, plot = FALSE, randomize = FALSE)

  trait <- scale(colorPCs[mapped_tree$tip.label,i])

  model_res <- brownie.lite(mapped_tree, trait, maxit = maxit)
  
  model_res$convergence <- model_res$convergence == "Optimization has converged."
  
  repeat {
    
    if (model_res$convergence || maxit >= max_iter) {
      break
    }
  
    maxit <- maxit * 2
    cat("Model did not converge. Increasing max iterations to", maxit, "\n")
    
    if (maxit > max_iter) {
      maxit <- max_iter
    }
    
    model_res <- brownie.lite(mapped_tree, trait, maxit = maxit)
    model_res$convergence <- model_res$convergence == "Optimization has converged."
  }
  
  # Final convergence status
  if (!model_res$convergence) {
    cat("Warning: Model did not converge even after reaching max iterations.\n")
  } else {
    if (model_res$P.chisq < 0.05) {
      cat("The multi-rate model is significant! Performing simulations...\n")
      
      sim_res <- parallel::mclapply(1:n_simulations, function(j) {
        
        # Simulate tree with l1ou2simmap
        sim_tree <- l1ou2simmap(fit, convergent.only = TRUE, plot = FALSE, randomize = TRUE)
        
        # Run brownie.lite on the simulated tree
        res <- brownie.lite(sim_tree, trait)
        
        return(res)
      }, mc.cores = num_cores)  
      
      observed_logLik <- model_res$logL.multiple
      simulated_logLik <- sapply(sim_res, function(result) result$logL.multiple)
      
      p_value <- sum(simulated_logLik >= observed_logLik) / length(simulated_logLik)
      
      if (p_value < 0.05) {
        cat("The observed model is significantly better than the simulated models.\n")
      } else {
        cat("The observed model is not significantly better than the simulated models.\n")
      }
      
    } else {
      cat("The multi-rate model is not significant.\n")
      sim_res <- NULL
      p_value <- NULL
    }
    color_rates[[paste0("PC", i)]] <- list(obs = model_res, sim = sim_res, p.val = p_value)
  }
}

saveRDS(color_rates, paste0("./rdata/", ifelse(convergent_only, "conv", "all"), "_colorPCs_rate_models.rds"))
sig_col <- sapply(color_rates, function(i) !is.null(i$sim))
# plot
pdf("./new_figures/l1ou/colorPCs_multi_rate.pdf", width = 14, height = 6)
par(mfrow = c(1, length(color_rates[sig_col])))
lapply(names(color_rates[sig_col]), function(model) {
  plot.model_rates(color_rates[[model]], cex = c(1.1,1.1,1.3,1.9), legend.position = "topright", title = paste0("Multi-rate Evolution of Color Convergence\n(", model, ")"))
})
dev.off()

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
  
  pdf(paste0("./new_figures/l1ou/morphoPC_", i, "_independent_shifts_", criterion, "_new.pdf"), width = 20, height = 25)
  ew <- rep(1,nEdges) # to set default edge width of 1
  ew[fit_ind$shift.configuration] <- 3 # to widen edges with a shift
  plot(fit_ind, cex=0.5, label.offset=0.02, edge.width=ew)
  dev.off()
  
  # save model
  saveRDS(fit_ind, paste0("./rdata/morphoPC_", i, "_damsel_l1ou_estimated_ind_shifts_", criterion, "_new.rds"))
  
  ## then detect which of these shifts are convergent:
  fit_conv <- estimate_convergent_regimes(fit_ind, criterion="pBIC", nCores = 5)
  print(configuration_ic(trait.data$tree, fit_conv$Y, fit_conv$shift.configuration, criterion=criterion))
  
  pdf(paste0("./new_figures/l1ou/morphoPC_", i, "_convergent_shifts_", criterion, "_new.pdf"), width = 20, height = 25)
  ew <- rep(1,nEdges) # to set default edge width of 1
  ew[fit_conv$shift.configuration] <- 3 # to widen edges with a shift
  plot(fit_conv, cex=0.5, label.offset=0.02, edge.width=ew)
  dev.off()
  
  # save model
  saveRDS(fit_conv, paste0("./rdata/morphoPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, "_new.rds"))
  
}

# Custom plots:
###############
pcs <- 1:8 # first plot all, then select those with convergent regimes
pcs <- c(1,3,4,7)
convergent_only = FALSE
criterion = "pBIC" 

pdf(paste0("./new_figures/l1ou/l1ou_morpho_", ifelse(convergent_only, "conv", "all"), "_pBic_PC1-8.pdf"), width = 4 * length(pcs), height = 10)

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
for (i in pcs) {
  fit_conv <- readRDS(paste0("./rdata/morphoPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds"))
  
  if (i == 1) { 
    par(mar = c(2, 1, 4, 2))  # smaller left margin for first PC plot
  } else {
    par(mar = c(2, 4, 4, 2))  # regular margins for the rest
  }
  
  plot_shift_optima(fit_conv, 
                    color.by = "regime", 
                    convergent.only = convergent_only,
                    shift.edge.width = 5, 
                    cex.value = 1.1, 
                    cex = 0.5, 
                    show.tip.label = FALSE, 
                    title = paste0("morpho PC", i), 
                    legend.position = "topleft",
                    cex.legend = 1)
}

dev.off()

# test differential trait evolutionary rate across l1ou regimes
###############################################################
pcs <- 1:8 
pcs <- c(1,3,4,7)
init_iter <- 2000
max_iter <- 10000
convergent_only = FALSE
num_cores <- 5
n_simulations <- 100
morpho_rates <- list()
for (i in pcs) {
  
  maxit <- init_iter
  
  message("Processing morpho PC", i, "...")
          
  fit <- readRDS(paste0("./rdata/morphoPC_", i, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds"))
  
  mapped_tree <- l1ou2simmap(fit, convergent.only = convergent_only, plot = FALSE, randomize = FALSE)

  trait <- scale(morphoPCs[mapped_tree$tip.label,i])
  
  model_res <- brownie.lite(sim_tree, trait, maxit = maxit)

  model_res$convergence <- model_res$convergence == "Optimization has converged."
  
  repeat {
    
    if (model_res$convergence || maxit >= max_iter) {
      break
    }
    
    maxit <- maxit * 2
    cat("Model did not converge. Increasing max iterations to", maxit, "\n")
    
    if (maxit > max_iter) {
      maxit <- max_iter
    }

    model_res <- brownie.lite(mapped_tree, trait, maxit = maxit)
    model_res$convergence <- model_res$convergence == "Optimization has converged."
  }
  
  # Final convergence status
  if (!model_res$convergence) {
    cat("Warning: Model did not converge even after reaching max iterations.\n")
  } else {
    if (model_res$P.chisq < 0.05) {
      cat("The multi-rate model is significant! Performing simulations...\n")
      
      sim_res <- parallel::mclapply(1:n_simulations, function(j) {
       
        sim_tree <- l1ou2simmap(fit, convergent.only = convergent_only, plot = FALSE, randomize = TRUE)
          
        res <- brownie.lite(sim_tree, trait)
        
        return(res)
      }, mc.cores = num_cores)
        
      # Remove failed runs
      sim_res <- Filter(Negate(is.null), sim_res)
      
      observed_logLik <- model_res$logL.multiple
      simulated_logLik <- sapply(sim_res, function(result) result$logL.multiple)
      
      p_value <- sum(simulated_logLik <= observed_logLik) / length(simulated_logLik)
  
      if (p_value < 0.05) {
        cat("The observed model is significantly better than the simulated models.\n")
      } else {
        cat("The observed model is not significantly better than the simulated models.\n")
      }
      
    } else {
      cat("The multi-rate model is not significant.\n")
      sim_res <- NULL
      p_value <- NULL
    }
  
    morpho_rates[[paste0("PC", i)]] <- list(obs = model_res, sim = sim_res, p.val = p_value)
  }
}

saveRDS(morpho_rates, paste0("./rdata/", ifelse(convergent_only, "conv", "all"), "_morphoPCs_rate_models.rds"))

sig_morpho <- sapply(morpho_rates, function(i) !is.null(i$sim))
# plot
pdf("./new_figures/l1ou/morphoPCs_multi_rate.pdf", width = 12, height = 4)
par(mfrow = c(1, length(morpho_rates[sig_morpho])))
lapply(names(morpho_rates[sig_morpho]), function(model) {
  plot.model_rates(morpho_rates[[model]], cex = c(1, 1, 1.2, 1.8), legend.position = "topright", title = paste0("Multi-rate Evolution of Morphological Convergence\n(", model, ")"))
})
dev.off()


morpho_rates[[which(sapply(morpho_rates, function(i) (i$p.val <= 0.05)))]]$obs

