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
dir.create("new_results/multivariate")
dir.create("new_figures/multivariate")

col.palette <- sample(unique(as.character(sapply(c("Zissou1", "Darjeeling1", "Darjeeling2", "FantasticFox1"), function(x) wes_palette(x, 5)))))
ncores <- 5
color_analysis = TRUE
morpho_analysis = TRUE

# load input data and custom functions
source("./Rscripts/A0.Input_data.R")

## Multivariate analyses testing colPCA and morphoPCA association to a trait
#################################################################################
for (trait in discrete_traits){
  if (file.exists(paste0("./rdata/ASR_", trait, ".rda"))){
    ## Distance matrices
    load(paste0("./rdata/ASR_", trait, ".rda"))
    
    trait.vec = setNames(eco_traits_subset[, trait], rownames(eco_traits_subset))[tree$tip.label]
    
    cat("Working on trait:", trait, "...\n")
    tree <- fit.marginal$phy
    tree$node.label <- levels(trait.vec)[tree$node.label]
    n <- Ntip(tree)
  
    # color
    if (color_analysis) {
      ###########################
      ## Color vs discrete traits
      ###########################
      cat("Running Multivariate Phylogenetic Generalized Least Squares on color PCA, testing ", trait, "...\n")
      
      ## MVPGLS
      dat <- list(color = as.matrix(colorPCs), trait = trait.vec)
      
      cat("BM...\n")
      model_BM1 <- mvgls(color ~ 1, data = dat, tree = tree, model = "BM")
      model_BM_trait <- mvgls(color~trait, data = dat, tree = tree, model = "BM")
      
      cat("Lambda...\n")
      model_lambda1 <- mvgls(color ~ 1, data = dat, tree = tree, model = "lambda")
      model_lambda_trait <- mvgls(color~trait, data = dat, tree = tree, model = "lambda")
      
      cat("OU...\n")
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
      model_EB1 <- mvgls(color~1, data = dat, tree = tree, model = "EB")
      model_EB_trait <- mvgls(color~trait, data = dat, tree = tree, model = "EB")
      
      trait_col.models <- list(
        BM1 = model_BM1, 
        lambda1 = model_lambda1,
        OU1 = model_OU1, 
        EB1 = model_EB1,
        BM_trait = model_BM_trait,
        lambda_trait = model_lambda_trait,
        OU_trait = model_OU_trait,
        EB_trait = model_EB_trait
        )
      
      res <- lapply(trait_col.models, EIC, nbcores = ncores)
      
      mvgls.col.restab <- data.frame(model = names(res),
                                     logLik = sapply(res, function(i) i$LogLikelihood),
                                     k = sapply(res, function(i) i$p),
                                     se = sapply(res, function(i) i$se),
                                     EIC = eic <- sapply(res, function(i) i$EIC),
                                     delta.EIC = eic-min(eic),
                                     weight = unclass(aic.w(eic)))
      
      eic_values <- sapply(res, function(i) i$EIC)
      best_model_name <- names(res)[which.min(eic_values)]
      best_model.col <- get(paste0("model_", names(res)[which.min(sapply(res, function(i) i$EIC))]))
      
      # Run MANOVA test
      if (grepl("trait", best_model_name)) {
        rownames(best_model.col$coefficients) <- levels(trait.vec)
        manova_test <- manova.gls(best_model.col, nperm = 999, test = "Pillai", verbose = TRUE)
        mvgls.col.restab[mvgls.col.restab$model == best_model_name, "Pillai"] <- manova_test$stat
        mvgls.col.restab[mvgls.col.restab$model == best_model_name, "p-value"] <- manova_test$pvalue
      }
      
      # write table to hard drive
      cat("Writting discrete color models result table on hard disk...\n")
      write.table(mvgls.col.restab, file = file.path(res_folder, paste0("multivariate/colorPCs.vs.", trait, "_mvgls.models.fit.restable.txt")), 
                  quote = FALSE, sep = "\t", row.names = T, col.names = T)
      
      saveRDS(trait_col.models, file = paste0("./rdata/colorPCs.vs.", trait, "_mvgls.models.rds"))
    }
    
    # clean model objects
    rm(list = ls(pattern = "^model_"))
    
    # morpho
    if (morpho_analysis) {
      ################################
      ## Morphology vs discrete traits
      ################################
      cat("Running Multivariate Phylogenetic Generalized Least Squares on morpho PCA, testing ", trait, "...\n")
      ## MVPGLS
      dat <- list(morpho = as.matrix(morphoPCs), trait = trait.vec)
      
      cat("BM...\n")
      model_BM1 <- mvgls(morpho ~ 1, data = dat, tree = tree, model = "BM")
      model_BM_trait <- mvgls(morpho~trait, data = dat, tree = tree, model = "BM")
      
      cat("Lambda...\n")
      model_lambda1 <- mvgls(morpho ~ 1, data = dat, tree = tree, model = "lambda")
      model_lambda_trait <- mvgls(morpho~trait, data = dat, tree = tree, model = "lambda")
      
      cat("OU...\n")
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
      model_EB1 <- mvgls(morpho~1, data = dat, tree = tree, model = "EB")
      model_EB_trait <- mvgls(morpho~trait, data = dat, tree = tree, model = "EB")
      
      trait_morph.models <- list(
        BM1 = model_BM1, 
        lambda1 = model_lambda1,
        OU1 = model_OU1, 
        EB1 = model_EB1,
        BM_trait = model_BM_trait,
        lambda_trait = model_lambda_trait,
        OU_trait = model_OU_trait,
        EB_trait = model_EB_trait
      )
      
      res <- lapply(trait_morph.models, EIC, nbcores = ncores)
      
      mvgls.morph.restab <- data.frame(model = names(res),
                                       logLik = sapply(res, function(i) i$LogLikelihood),
                                       k = sapply(res, function(i) i$p),
                                       se = sapply(res, function(i) i$se),
                                       EIC = eic <- sapply(res, function(i) i$EIC),
                                       delta.EIC = eic-min(eic),
                                       weight = unclass(aic.w(eic)))
      
      eic_values <- sapply(res, function(i) i$EIC)
      best_model_name <- names(res)[which.min(eic_values)]
      best_model.morph <- get(paste0("model_", names(res)[which.min(eic_values)]))
      
      # Run MANOVA test
      if (grepl("trait", best_model_name)) {
        rownames(best_model.morph$coefficients) <- levels(trait.vec)
        manova_test <- manova.gls(best_model.morph, nperm = 999, test = "Pillai", verbose = TRUE)
        mvgls.morph.restab[mvgls.morph.restab$model == best_model_name, "Pillai"] <- manova_test$stat
        mvgls.morph.restab[mvgls.morph.restab$model == best_model_name, "p-value"] <- manova_test$pvalue
      }
      
      # write table to hard drive
      cat("Writting discrete morpho models result table on hard disk...\n")
      write.table(mvgls.morph.restab, file = file.path(res_folder, paste0("multivariate/morphoPCs.vs.", trait, "_mvgls.models.fit.restable.txt")), 
                  quote = FALSE, sep = "\t", row.names = T, col.names = T)
      
      saveRDS(trait_morph.models, file = paste0("./rdata/morphoPCs.vs.", trait, "_mvgls.models.rds"))
      
      # clean model objects
      rm(list = ls(pattern = "^model_"))
    }
  } else {
    warning(trait, "rda file not found!! Skipping analysis...\n")
  }
}

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
    trait_col.models <- readRDS(file.path("rdata", paste0("colorPCs.vs.", trait, "_mvgls.models.rds")))
    trait_morph.models <- readRDS(file.path("rdata", paste0("morphoPCs.vs.", trait, "_mvgls.models.rds")))
    
    ## load res tables
    col_restable <- read.table( file.path(res_folder, paste0("colorPCs.vs.", trait, "_mvgls.models.fit.restable.txt")))
    morph_restable <- read.table( file.path(res_folder, paste0("morphoPCs.vs.", trait, "_mvgls.models.fit.restable.txt")))
    
    ## get best models
    col_best_model_name <- col_restable$model[which.min(col_restable$EIC)]
    morpho_best_model_name <- morph_restable$model[which.min(morph_restable$EIC)]
    
    best_model.col <- trait_col.models[[col_best_model_name]]
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
    
    pdf(paste0("./new_figures/multivariate/col-morph.vs.", trait, ".pdf"), height = 15, width = 10)
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
    
    pdf(paste0("./new_figures/multivariate/ColorPCA1-4.vs.", trait, ".pdf"), width = 20, height = 12)
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
    
    pdf(paste0("./new_figures/multivariate/ColorPCA5-8.vs.", trait, ".pdf"), width = 20, height = 12)
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
    
    pdf(paste0("./new_figures/multivariate/MorphPCA1-4.vs.", trait, ".pdf"), width = 20, height = 12)
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
    
    pdf(paste0("./new_figures/multivariate/MorphPCA5-8.vs.", trait, ".pdf"), width = 20, height = 12)
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
sapply(col_traitList_best_models[c("DietEcotype", "Farming", "habitat1")], function(x) x$param)
sapply(morpho_traitList_best_models[c("DietEcotype", "Farming", "habitat1")], function(x) x$param)
