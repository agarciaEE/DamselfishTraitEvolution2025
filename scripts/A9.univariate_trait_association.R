lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", 
         "lme4", "mgcv", "phytools", "corHMM",
         "lmtest", "ColorAR", "castor", "tibble",
         "ggplot2", "ggpubr", "gridExtra",
         "knitr", "kableExtra")
sapply(lib, library, character.only = T)

#=============================#
#   Settings
#=============================#
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/")

res_folder <- "new_results"
dir.create(res_folder)
dir.create("new_figures/univariate/")
dir.create("new_results/univariate/")

nsim <- 1000

# load input data and custom functions
source("./Rscripts/A0.Input_data.R")

#=============================#
# Function: Run pGLS models
#=============================#
run_gls_models <- function(PC_matrix, traits_df, tree, prefix) {
  result <- list()
  for (i in 1:ncol(PC_matrix)) {
    response <- scale(PC_matrix[, i])
    model_data <- data.frame(response, traits_df)
    rownames(model_data) <- rownames(traits_df)
    
    gls_model <- gls(response ~ DietEcotype + habitat1,
                     correlation = corBrownian(phy = tree, form = ~ species),
                     data = model_data)
    result[[paste0(prefix, i)]] <- gls_model
  }
  return(result)
}

# Run for both color and morpho PCs
col_pgls_models <- run_gls_models(colorPCs, eco_traits_subset, tree, "colorPC")
morpho_pgls_models <- run_gls_models(morphoPCs, eco_traits_subset, tree, "morphoPC")

#=============================#
# Summarize GLS Results
#=============================#
summarize_gls <- function(model_list) {
  res <- lapply(model_list, function(model) {
    s <- summary(model)
    coefs <- s$tTable[,"Value"]
    pvals <- s$tTable[, "p-value"]
    return(data.frame(coefs, pvals))
  }) |> do.call("rbind", args = _)
  res$PC <- sapply(rownames(res), function(i) strsplit(i, "\\.")[[1]][1])
  res$trait <- unlist(sapply(rownames(res), function(i) {
    idx <- sapply(discrete_traits, grepl, strsplit(i, "\\.")[[1]][2])
    ifelse(any(idx), discrete_traits[idx], "(Intercept)")
    }))
  res$category <- unlist(sapply(rownames(res), function(i) {
    name <- strsplit(i, "\\.")[[1]][2]
    idx <- sapply(discrete_traits, grepl, name)
    ifelse(any(idx), sub(discrete_traits[idx], "", name), "Benthic-Sea anemone")
  }))
  res
}

color_results <- summarize_gls(col_pgls_models)
morpho_results <- summarize_gls(morpho_pgls_models)

write.csv(color_results, file = file.path(res_folder, "/univariate/colorPC_trait_gls.csv"), row.names = FALSE)
write.csv(morpho_results, file = file.path(res_folder, "/univariate/morphoPC_trait_gls.csv"), row.names = FALSE)

#=============================#
# Phylogenetic ANOVA
#=============================#
#color
traits_list <- colnames(colorPCs)

color_phylANOVA_results <- expand.grid(trait = traits_list, eco_trait = discrete_traits)
color_phylANOVA_results$F <- NA
color_phylANOVA_results$pvals <- NA

for (i in seq_len(nrow(color_phylANOVA_results))) {
  trait <- color_phylANOVA_results$trait[i]
  eco <- color_phylANOVA_results$eco_trait[i]
  
  y <- colorPCs[tree$tip.label, trait]
  group <- eco_traits_subset[tree$tip.label, eco]
  
  res <- tryCatch(phylANOVA(tree, group, y, nsim = nsim), error = function(e) NULL)
  if (!is.null(res)) {
    color_phylANOVA_results$F[i] <- res$F
    color_phylANOVA_results$pvals[i] <- res$Pf
  }
}

# morpho
traits_list <- colnames(morphoPCs)

morpho_phylANOVA_results <- expand.grid(trait = traits_list, eco_trait = discrete_traits)
morpho_phylANOVA_results$F <- NA
morpho_phylANOVA_results$pvals <- NA

for (i in seq_len(nrow(morpho_phylANOVA_results))) {
  trait <- morpho_phylANOVA_results$trait[i]
  eco <- morpho_phylANOVA_results$eco_trait[i]
  
  y <- morphoPCs[tree$tip.label, trait]
  group <- eco_traits_subset[tree$tip.label, eco]
  
  res <- tryCatch(phylANOVA(tree, group, y, nsim = nsim), error = function(e) NULL)
  if (!is.null(res)) {
    morpho_phylANOVA_results$F[i] <- res$F
    morpho_phylANOVA_results$pvals[i] <- res$Pf
  }
}

write.csv(color_phylANOVA_results, file = file.path(res_folder, "/univariate/colorPC_trait_phylANOVA.csv"), row.names = FALSE)
write.csv(morpho_phylANOVA_results, file = file.path(res_folder, "/univariate/morphoPC_trait_phylANOVA.csv"), row.names = FALSE)

#=============================#
# Output Summary
#=============================#
display_result <- function(df, p.value = 0.05, title = "") {
  # Filter only significant (p < 0.05) results
  df_filtered <- data.frame(df) %>% dplyr::filter(pvals < p.value)
  # Display 
  kable(df_filtered, digits = 4, row.names = FALSE, caption = title) %>%
    kable_styling(full_width = FALSE)
}

print("pGLS colorPC associations:")
print(color_results)
display_result(color_results, title = "Significant pGLS colorPC ~ trait  associations")

print("pGLS morphoPC associations:")
print(morpho_results)
display_result(morpho_results, title = "Significant pGLS morphoPC ~ trait associations")

print("Phylogenetic ANOVA colorPC summary:")
print(color_phylANOVA_results)
display_result(color_phylANOVA_results, title = "Significant Phylogenetic ANOVA colorPC ~ trait  associations")

print("Phylogenetic ANOVA morphoPC summary:")
print(morpho_phylANOVA_results)
display_result(morpho_phylANOVA_results, title = "Significant Phylogenetic ANOVA colorPC ~ trait  associations")
