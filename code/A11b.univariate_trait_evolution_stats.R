## ============================
##  Setup
## ============================
libs <- c("dplyr","tidyr","purrr","stringr","ggplot2","readr",
          "forcats","kableExtra","glue", "ggpmisc", "ggtext")
invisible(lapply(libs, require, character.only = TRUE))

setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")
source("./scripts/A1.Input_data.R")

dir.create("figures", showWarnings = FALSE)

## ============================
## Helper Functions
## ============================

# --------------------------
## Helper to standardize model names for readability
# --------------------------
standardize_model_names <- function(x) {
  x <- as.character(x)
  x[x == "white"]            <- "WN"
  x[x == "BM"]               <- "BM"
  x[x == "EB"]               <- "EB"
  x[x == "OU"]               <- "OU"
  x[grepl("^OUM_Diet", x)]   <- "OUM_Diet"
  x[grepl("^OUM_habitat",x)] <- "OUM_Habitat"
  x[grepl("^OUM_Farming",x)] <- "OUM_Symbiosis"
  x[x == "MC"]               <- "MC"
  x[x == "DDlin"]            <- "DDlin"
  x[x == "DDexp"]            <- "DDexp"
  x
}

# --------------------------
## Compact "family" label 
# --------------------------
model_family <- function(x) {
  x <- standardize_model_names(x)
  ifelse(x %in% c("OU","OUM_Diet","OUM_Habitat","OUM_Symbiosis"), "OU-type", 
         ifelse(x %in% c("DDlin", "DDexp", "MC"),"DD-type", x))
}

# --------------------------
## Combine different models list (base models + DD models)
# --------------------------
combine_modelsList <- function(base_models, dd_models) {
  # Ensure both inputs are lists
  stopifnot(is.list(base_models), is.list(dd_models))
  
  # Get common PC (or trait) names
  common_pcs <- intersect(names(base_models), names(dd_models))
  
  # Combine corresponding model sublists
  combined_list <- lapply(common_pcs, function(pc) {
    merged <- c(base_models[[pc]], dd_models[[pc]])
    # remove NULL elements if any
    merged[!sapply(merged, is.null)]
  })
  names(combined_list) <- common_pcs
  
  return(combined_list)
}

# --------------------------
## Extract AIC, k, n_obs, AICc from a model object
# --------------------------
extract_model_info <- function(model, tree, name = NA_character_) {
  if (is.null(model)) return(NULL)
  
  # --- Try to detect AIC ---
  AIC_val <- NA
  if (!is.null(model$aic)) AIC_val <- model$aic
  else if (!is.null(model$AIC)) AIC_val <- model$AIC
  else if (is.list(model$opt) && !is.null(model$opt$aic)) AIC_val <- model$opt$aic
  
  # --- Try to detect k ---
  k_val <- NA
  if (!is.null(model$opt$k)) k_val <- length(model$opt$k)
  else if (!is.null(model$param.count)) k_val <- model$param.count
  else if (!is.null(model$free.parameters)) k_val <- model$free.parameters
  
  # --- Detect n_obs ---
  n_val <- ape::Ntip(tree)
  
  # --- Try to detect AICc ---
  AICc_val <- NA
  if (!is.null(model$aicc)) AIC_val <- model$aicc
  else if (!is.null(model$AICc)) AIC_val <- model$AICc
  else if (is.list(model$opt) && !is.null(model$opt$aicc)) AIC_val <- model$opt$aicc
  
  # --- If not, compute AICc ---
  if (any(is.na(AICc_val))) {
    idx = is.na(AICc_val)
    AICc_val[idx] =  AIC_val[idx] + (2 * k_val[idx] * (k_val[idx] + 1)) / (n_val - k_val[idx] - 1)
  }
  
  tibble(model = name, AIC = AIC_val, k = k_val, n_obs = n_val, AICc = AICc_val)
}

# --------------------------
##  Extract AIC, k, n_obs, and AICc directly from models
# --------------------------
make_models_table <- function(trait_models_list, tree, 
                              trait_label = c("Color","Morphology"),
                              regional = FALSE, region_names = NULL) {
  trait_label <- match.arg(trait_label)
  pcs <- names(trait_models_list)
  
  purrr::map_dfr(pcs, function(pcname) {
    message("Processing ", pcname)
    
    # Extract all models for this PC
    all_models <- trait_models_list[[pcname]]
    all_models <- all_models[!vapply(all_models, is.null, logical(1))]
    
    if (length(all_models) == 0) return(NULL)
    
    if (regional) {
      ## Regional models
      purrr::imap_dfr(all_models, function(models, region_code) {
        if (is.null(models) || !is.list(models)) return(NULL)
        models <- models[!vapply(models, is.null, logical(1))]
        if (length(models) == 0) return(NULL)
        
        # extract parameters from all models
        df <- purrr::imap_dfr(models, ~extract_model_info(.x, tree = tree, name = .y))
        if (nrow(df) == 0) return(NULL)
        
        # compute AICc-based weights
        df <- df %>%
          mutate(
            model = standardize_model_names(model),
            PC = pcname,
            region = if (!is.null(region_names)) region_names[region_code] else region_code,
            TraitSet = trait_label
          ) %>%
          group_by(PC, region) %>%
          mutate(
            DeltaAICc = AICc - min(AICc, na.rm = TRUE),
            Weight = exp(-0.5 * DeltaAICc),
            Weight = Weight / sum(Weight, na.rm = TRUE)
          ) %>%
          ungroup()
        
        df
      })
    } else{
      ## Global models
      df <- purrr::imap_dfr(all_models, ~extract_model_info(.x, tree = tree, name = .y))
      
      df <- df %>%
        mutate(model = standardize_model_names(model),
               PC = pcname,
               TraitSet = trait_label) %>%
        group_by(PC) %>%
        mutate(DeltaAICc = AICc - min(AICc, na.rm = TRUE),
               Weight = exp(-0.5 * DeltaAICc),
               Weight = Weight / sum(Weight, na.rm = TRUE)) %>%
        ungroup()
      
      df
    }
    
  })
}

# --------------------------
##  Overall support patterns – shows fraction of PCs where each model is best
# --------------------------
summarize_best_counts_all <- function(df) {
  df %>%
    group_by(TraitSet, PC) %>%
    slice_max(order_by = Weight, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Model = factor(standardize_model_names(model), levels = model_levels)) %>%
    count(TraitSet, Model, name = "n_PC") %>%
    group_by(TraitSet) %>%
    mutate(frac = n_PC / sum(n_PC)) %>%
    ungroup() %>%
    # ensure every model appears per TraitSet, even if 0
    tidyr::complete(TraitSet, Model, fill = list(n_PC = 0, frac = 0)) %>%
    mutate(Model = factor(Model, levels = model_levels))
}

# --------------------------
##  Prepare data set for lollipop plot
# --------------------------
make_lolli_df <- function(df, trait_set_label = "Color", regional = FALSE, weight = "Weight") {
  
  df %>%
    dplyr::mutate(
      TraitSet = trait_set_label,
      Model = standardize_model_names(model)
    ) %>%
    {
      # Group by region if regional, otherwise only by PC
      if (regional && "region" %in% colnames(.)) {
        group_by(., PC, region)
      } else {
        group_by(., PC)
      }
    } %>%
    # Recompute AICc differences and weights (for safety)
    dplyr::mutate(
      best_model = Model[which.max(.data[[weight]])],
      ER = .data[[weight]][Model == best_model] / .data[[weight]],
      ER = ifelse(Model == best_model, Inf, ER),
      cat = dplyr::case_when(
        Model == best_model ~ "Best Model",
        ER < 3   ~ "ER < 3 (indistinguishable)",
        ER < 20  ~ "3-20 (positive)",
        ER < 150 ~ "20-150 (strong)",
        TRUE     ~ "> 150 (very strong)"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      PC = sub("color|morpho", "", PC),
      # Order factors for plotting
      Model = factor(Model, levels = rev(model_levels)),
      Region = if (regional && "region" %in% colnames(.))
        factor(region, levels = region_names)
      else
        NULL
    )
}

# --------------------------
# Plotting: lollipop per PC
# --------------------------
plot_lolli <- function(df_prepped, regional = FALSE, weight = "Weight") {
  
  # ensure consistent factor levels
  df_prepped$cat <- factor(
    df_prepped$cat,
    levels = c(
      "Best Model",
      "ER < 3 (indistinguishable)",
      "3-20 (positive)",
      "20-150 (strong)",
      "> 150 (very strong)"
    )
  )
  
  df_prepped$Model <- factor(
    df_prepped$Model,
    levels = c(
      "MC", "DDexp", "DDlin",
      "OUM_Symbiosis", "OUM_Diet", "OUM_Habitat",
      "OU", "EB", "BM", "WN"
    )
  )
  
  # background rectangles for model group blocks
  rect_df <- tibble::tibble(
    xmin = c(0.5, 3.5, 6.5),
    xmax = c(3.5, 6.5, 9.5),
    fill = c("#deebf7", "#fee6ce", "#f7f7f7"),
    group = c("DD models", "multi-OU models", "simple models")
  )
  
  p <- ggplot(df_prepped, aes(x = Model, y = .data[[weight]], color = cat)) +
    geom_rect(
      data = rect_df,
      inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = group),
      alpha = 0.5, colour = NA
    ) +
    geom_segment(aes(xend = Model, y = 0, yend = .data[[weight]]), linewidth = 0.6) +
    geom_point(aes(size = scales::rescale(.data[[weight]], to = c(1, 1.5))), shape = 16) +
    coord_flip() +
    scale_size(guide = "none") +
    geom_vline(xintercept = c(3.5, 6.5), linetype = "dashed", color = "grey70", linewidth = 0.4) +
    scale_fill_manual(
      name = NULL,
      values = setNames(rect_df$fill, rect_df$group),
      guide = "none"
    ) +
    scale_color_manual(
      name = "Support vs. best model",
      values = c(
        "Best Model"                   = "#FFD700",
        "ER < 3 (indistinguishable)"   = "#b2df8a",
        "3-20 (positive)"              = "#a6cee3",
        "20-150 (strong)"              = "#1f78b4",
        "> 150 (very strong)"          = "#08306b"
      ),
      drop = FALSE
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = ifelse(weight == "Weight", "Akaike weight", "Normalized Akaike weigth")) +
    theme_minimal(base_size = ifelse(regional, 18, 13)) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
  
  # dynamic faceting
  if (regional && "region" %in% colnames(df_prepped)) {
    p <- p + facet_grid(PC ~ region)
  } else {
    p <- p + facet_wrap(~ PC, ncol = 2)
  }
  
  return(p)
}

# --------------------------
# Plot to show how model weights vary across structural metrics
# --------------------------
plot_struct_bias <- function(df, metric_labels, family_colors, type = "all") {
  
  df$ModelFamily <- factor(df$ModelFamily, levels = model_families)
  df_all <- df
  if (type == "best") {
    df <- df %>% filter(DeltaAICc == 0)
  }
  
  model_with_n <- df %>%
    group_by(ModelFamily, metric) %>%
    count() %>%
    filter(n > 1) %>%
    pull(ModelFamily) %>% unique()
  
  df <- df %>% filter(ModelFamily %in% model_with_n) %>%
    mutate(ModelFamily = factor(ModelFamily, levels = model_families))
  
  ggplot(df_all, aes(x = value, y = Weight, color = ModelFamily, fill = ModelFamily)) +
    geom_point(alpha = 0) +
    geom_point(data = df, alpha = 0.6, size = 2.5) +
    geom_smooth(data = df, method = "lm", se = TRUE, alpha = 0.25, linewidth = 0.8, show.legend = FALSE) +

    # Equation per model class 
    stat_poly_eq(data = df,
      aes(
        label = paste(model_with_n, "~",
          after_stat(eq.label),
          after_stat(rr.label),
          after_stat(p.value.label),
          sep = "~~~"
        ),
        group = ModelFamily
      ),
      formula = y ~ x,
      parse = TRUE,
      size = 3,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    
    scale_color_manual(values = family_colors, labels = names(family_colors), name = "Model class", drop = FALSE) +
    scale_fill_manual(values = family_colors, labels = names(family_colors), guide = "none", drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    coord_cartesian(ylim = c(0, 1.3)) +
    guides(
      color = guide_legend(
        override.aes = list(size = 4, alpha = 1),
        title = "Model class"
      )
    ) +
    facet_wrap(
      ~ metric,
      scales = "free_x",
      labeller = labeller(metric = metric_labels)
    ) +
    
    labs(
      x = "Metric value",
      y = paste(ifelse(type == "best", "Best", "Relative"), "model support (Akaike weight)")
    ) +
    
    theme_light(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(size = 13, face = "bold", color = "grey20"),
      plot.margin = margin(30, 20, 25, 20),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 13),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      legend.key.size = unit(0.8, "lines"),
      panel.spacing = unit(1, "lines")
    )
}

# --------------------------
# Get model parameters
# --------------------------
.model_param <- function(m, name = NULL) {

  if ("theta" %in% names(m)) {
    setNames(m$theta[,1], rep( ifelse(!is.null(name), name, "OUM"),
                                      length(m$theta[,1])))
  } else if ("S" %in% names(m)) {
    setNames(m$S,  ifelse(!is.null(name), name, "MC"))
  } else if ("b" %in% names(m)) {
    setNames(m$b,  ifelse(!is.null(name), name, "DDlin"))
  } else if ("r" %in% names(m)) {
    setNames(m$r,  ifelse(!is.null(name), name, "DDexp"))
  } else if ("a" %in% names(m$opt)) {
    setNames(m$opt$a, ifelse(!is.null(name), name, "EB"))
  } else if ("alpha" %in% names(m$opt)) {
    setNames(m$opt$alpha,  ifelse(!is.null(name), name, "OU"))
  }  else {
    NULL
  }
}

## ============================
##  Colors & Levels
## ============================

# Define color palette and labels for ALL models
model_levels <- c("WN", "BM", "EB", "OU", "OUM_Diet", "OUM_Habitat", "OUM_Symbiosis", "DDlin", "DDexp", "MC")

## Assign family-based colors
model_colors <- c(
  WN          = "lightgrey",
  BM          = "black",
  EB          = "green3",
  OU          = "#e41a1c",       # strong red
  OUM_Diet    = "#fb6a4a",       # lighter red
  OUM_Habitat = "#a50f15",       # dark red
  OUM_Symbiosis = "#fcae91",     # pale red
  DDlin       = "#6baed6",       # light blue
  DDexp       = "#2171b5",       # medium blue
  MC          = "#08306b"        # dark blue
)

model_families <- c("WN", "BM", "EB", "OU-type", "DD-type")

family_colors <- c(
  "WN" = "grey60",
  "BM" = "grey20",
  "EB" = "green3",
  "OU-type" = "#e6550d",
  "DD-type" = "#3182bd"
)

## Region name map (match your earlier code)
region_names <- setNames(c("IO","IAA","CPO","EPO","AO"), LETTERS[1:5])

metric_labels <- c(
  imbalance = "Tree imbalance (Colless Index)",
  K         = "Phylogenetic signal (Blomberg's K)",
  MNTD      = "Mean Nearest Taxon Distance",
  MPD       = "Mean Pairwise Distance",
  PD        = "Phylogenetic Diversity",
  trait_var = "Trait variance"
)

### ================================= ### 
###           Review results          ###
### ================================= ### 

## ============================
##  GLOBAL MODELS
## ============================

##  Load saved full model fits
## ============================
# === Color === 
color_global_models     <- readRDS("./rdata/uni_colorPC_global_traitevo_models.rds")

sapply(color_global_models, function(pc) unlist(lapply(names(pc), function(m) .model_param(pc[[m]], name = m))), simplify = FALSE)
# === Morpho === 
morpho_global_models     <- readRDS("./rdata/uni_morphoPC_global_traitevo_models.rds")

sapply(morpho_global_models, function(pc) unlist(lapply(names(pc), function(m) .model_param(pc[[m]], name = m))), simplify = FALSE)[c(1,3,6,7)]

##  Make models table
## ============================
color_global_df  <- make_models_table(color_global_models, tree, "Color")
morpho_global_df <- make_models_table(morpho_global_models, tree, "Morphology")

# Merge color + morphology global model summaries
global_summary_full <- bind_rows(color_global_df, morpho_global_df) %>%
  dplyr::select(
    TraitSet, PC, model, AIC, k, n_obs, AICc, DeltaAICc, Weight
  )

# Save compact models table
readr::write_csv(global_summary_full, "./results/univariate/TraitEvo_global_model_AICc_full.csv")

##  Overall support patterns 
## ============================
best_counts_all <- bind_rows(color_global_df, morpho_global_df) %>% summarize_best_counts_all()

p_fig1_all <- ggplot(best_counts_all, aes(x = TraitSet, y = frac, fill = Model)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(aes(label = ifelse(frac > 0, scales::percent(frac, accuracy = 1L), "")),
            position = position_stack(vjust = 0.5), size = 6) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, .02))) +
  scale_fill_manual(values = model_colors, breaks = model_levels, limits = model_levels, drop = FALSE) +
  labs(x = NULL, y = "Proportion of PCs (best-supported model)",
       fill = "Model") +
  theme_minimal(base_size = 20) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "right")

# Make Lollipop Plots
## ============================
color_lolli_df  <- make_lolli_df(color_global_df,  trait_set_label = "Color",       regional = FALSE)
morpho_lolli_df <- make_lolli_df(morpho_global_df, trait_set_label = "Morphology",  regional = FALSE)

p_color_lolli  <- plot_lolli(color_lolli_df,  regional = FALSE)
p_morph_lolli  <- plot_lolli(morpho_lolli_df, regional = FALSE)

combined_lolli <- ggarrange(
  p_color_lolli + 
    ggtitle("(a) Color PCs") + 
    theme(plot.title = element_text(hjust = -0.25, size = 14, face = "bold")),
  
  p_morph_lolli + 
    ggtitle("(b) Morphology PCs") + 
    theme(plot.title = element_text(hjust = -0.25, size = 14, face = "bold")),
  
  ncol = 2, 
  common.legend = TRUE, 
  legend = "bottom"
)

# Save high-res 
## ============================
ggsave("./figures/univariate/TraitEvol_overall_model_support.pdf", p_fig1_all, width = 12, height = 9)
ggsave("./figures/univariate/TraitEvol_overall_model_support.png", p_fig1_all, width = 12, height = 9, dpi = 300)
ggsave("./figures/univariate/TraitEvo_LolliPanels.pdf", combined_lolli, width = 14, height = 10)

## ============================
##  REGIONAL MODELS
## ============================

col_reg_models2 <- list()
for (pc in names(col_reg_models)) {
  for (reg in names(col_reg_models[[pc]])) {
    col_reg_models2[[pc]][[reg]] <- c(col_reg_models[[pc]][[reg]], list(EB = col_reg_EBs[[1]][[pc]]$fits[[reg]]))
    col_reg_models2[[pc]][[reg]] <- col_reg_models2[[pc]][[reg]][sub("Diet", "DietEcotype", model_levels)]
  }
}
saveRDS(col_reg_models2, "./rdata/uni_colorPC_reg_traitevo_models.rds")

morpho_reg_models2 <- list()
for (pc in names(morpho_reg_models)) {
  for (reg in names(morpho_reg_models[[pc]])) {
    morpho_reg_models2[[pc]][[reg]] <- c(morpho_reg_models[[pc]][[reg]], list(EB = morpho_reg_EBs[[1]][[pc]]$fits[[reg]]))
    morpho_reg_models2[[pc]][[reg]] <- morpho_reg_models2[[pc]][[reg]][sub("Diet", "DietEcotype", model_levels)]
  }
}
saveRDS(morpho_reg_models2, "./rdata/uni_morphoPC_reg_traitevo_models.rds")

##  Load saved full model fits
## ============================
# === Color === 
col_reg_models        <- readRDS("./rdata/uni_colorPC_reg_traitevo_models.rds")

# === Morpho === 
morpho_reg_models     <- readRDS("./rdata/uni_morphoPC_reg_traitevo_models.rds")

##  Make models table
## ==================
color_reg_df     <- make_models_table(col_reg_models,      tree, "Color",     regional = TRUE,  region_names)
morpho_reg_df    <- make_models_table(morpho_reg_models,   tree, "Morphology", regional = TRUE, region_names)

# Save compact models table
readr::write_csv(bind_rows(color_reg_df, morpho_reg_df),
                 "./results/univariate/TraitEvo_regional_model_AICc_full.csv")

## ============================
##  Regional metrics & best-model support
## ============================

# Region structural stats were accumulated during fitting
col_region_stats_df    <- readRDS("./rdata/uni_colorPC_tree_region_stats.rds")
morpho_region_stats_df <- readRDS("./rdata/uni_morphoPC_tree_region_stats.rds")

# ================================================================
# Combine ALL model weights with regional metrics
# ================================================================

# (1) Join color model fits with structural metrics
color_weights_full <- color_reg_df %>%
  left_join(
    col_region_stats_df %>% dplyr::rename(PC = trait) %>% select(PC, region, imbalance, K, MNTD, MPD, PD, trait_var),
    by = c("PC", "region")
  ) %>%
  mutate(
    ModelFamily = model_family(model),
    ModelFamily = factor(ModelFamily, levels = c("WN", "BM", "EB", "OU-type", "DD-type")),
    TraitSet = "Color"
  )

# (2) Join morphology model fits with its region metrics
morpho_weights_full <- morpho_reg_df %>%
  left_join(
    morpho_region_stats_df %>% dplyr::rename(PC = trait) %>% select(PC, region, imbalance, K, MNTD, MPD, PD, trait_var),
    by = c("PC", "region")
  ) %>%
  mutate(
    ModelFamily = model_family(model),
    ModelFamily = factor(ModelFamily, levels = c("WN", "BM", "EB", "OU-type", "DD-type")),
    TraitSet = "Morphology"
  )

# ================================================================
# Pivot metrics and prepare for plotting
# ================================================================

color_weights_df <- color_weights_full %>%
  pivot_longer(
    cols = all_of(names(metric_labels)),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value), Weight > 0)

morpho_weights_df <- morpho_weights_full %>%
  pivot_longer(
    cols = all_of(names(metric_labels)),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value), Weight > 0)

# ================================================================
# Plot model weights variation across structural metrics
# ================================================================
combined_struct_support <- plot_struct_bias(bind_rows(color_weights_df, morpho_weights_df), metric_labels, family_colors, type = "best")
ggsave("./figures/univariate/Struct_bias_regional_model_support.pdf", plot = combined_struct_support, width = 15, height = 10)

##  Regional support patterns 
## ============================
regional_best <- bind_rows(color_reg_df, morpho_reg_df) %>%
  group_by(TraitSet, region, PC) %>%
  slice_max(order_by = Weight, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(ModelFamily = model_family(model)) %>%
  count(TraitSet, region, ModelFamily, name = "n_PC") %>%
  group_by(TraitSet, region) %>%
  mutate(frac = n_PC / sum(n_PC)) %>%
  ungroup() %>%
  tidyr::complete(
    TraitSet,
    region      = region_names,      # ensure all regions
    ModelFamily = model_families,    # ensure all model families
    fill = list(n_PC = 0, frac = 0)
  ) %>%
  mutate(
    region      = factor(region,      levels = region_names),
    ModelFamily = factor(ModelFamily, levels = rev(model_families))
  )

p_fig2 <- ggplot(regional_best, aes(x = region, y = ModelFamily, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = scales::percent(frac, accuracy = 1L)), size = 6) +
  scale_fill_gradient(
    low = "#f0f0f0", high = "#377eb8",
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    x = "Region",
    y = "Model class",
    fill = "Proportion of PCs\n"
  ) +
  facet_wrap(~ TraitSet, ncol = 1) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid   = element_blank(),
    strip.text   = element_text(face = "bold")
  )

ggsave("./figures/univariate/TraitEvo_regional_model_support_heatmap.pdf", p_fig2, width = 12, height = 9)
ggsave("./figures/univariate/TraitEvo_regional_model_support_heatmap.png", p_fig2, width = 12, height = 9, dpi = 300)

##  Prepare lollipop regional data
# ================================
color_reg_lolli_df  <- make_lolli_df(color_reg_df,  trait_set_label = "Color",       regional = TRUE)
morpho_reg_lolli_df <- make_lolli_df(morpho_reg_df, trait_set_label = "Morphology",  regional = TRUE)

# Make Regional plots
# ================================
p_color_reg  <- plot_lolli(color_reg_lolli_df,  regional = TRUE) + theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1.5))
p_morpho_reg <-  plot_lolli(morpho_reg_lolli_df,  regional = TRUE) + theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1.5))

combined_reg_lolli <- ggarrange(
  p_color_reg + 
    ggtitle("(a) Color PCs") + 
    theme(plot.title = element_text(hjust = -0.1, size = 26, face = "bold"),
          strip.text = element_text(size = 24, face = "bold"),
          legend.title = element_text(size = 24, face = "bold"),
          legend.text  = element_text(size = 20),
          legend.key.size = unit(2.5, "lines")),
  
  p_morpho_reg + 
    ggtitle("(b) Morphology PCs") + 
    theme(plot.title = element_text(hjust = -0.1, size = 26, face = "bold"),
          strip.text = element_text(size = 24, face = "bold"),
          legend.title = element_text(size = 24, face = "bold"),
          legend.text  = element_text(size = 20),
          legend.key.size = unit(2.5, "lines")),
  
  nrow = 2, 
  common.legend = TRUE, 
  legend = "right"
) + guides(fill = guide_legend(label.theme = element_text(size = 22)))

# Save high-res version
# ================================
ggsave("./figures/univariate/TraitEvo_Color_reg_LolliPanels.pdf", p_color_reg, width = 15, height = 17)
ggsave("./figures/univariate/TraitEvo_Morpho_reg_LolliPanels.pdf", p_morpho_reg, width = 15, height = 17)
ggsave("./figures/univariate/TraitEvo_reg_LolliPanels.pdf", combined_reg_lolli, width = 25, height = 35)
