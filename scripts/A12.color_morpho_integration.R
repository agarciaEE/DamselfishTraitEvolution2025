lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", "phyloTop",
         "lme4", "mgcv", "phytools", "corHMM", "ggpmisc",
         "lmtest", "ColorAR", "castor", "tibble", "dplyr", "purrr",
         "knitr", "kableExtra", "scales", "cowplot", "png", "grid",
         "ggplot2", "ggpubr", "gridExtra", "pheatmap")
sapply(lib, library, character.only = T)

#=============================#
#   Settings
#=============================#
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

# load custom functions
source("./scripts/custom_functions.R")

res_folder <- "results"
dir.create(res_folder)
dir.create("results/col_morpho/")
dir.create("figures/col_morpho/")

col.palette <- sample(unique(as.character(sapply(c("Zissou1", "Darjeeling1", "Darjeeling2", "FantasticFox1"), function(x) wes_palette(x, 5)))))
ncores <- 5
geo_model = FALSE

# load input data and custom functions
source("./scripts/A1.Input_data.R")

#=============================#
# Pgls morpho vs color
#=============================#

# Make sure species are in the same order across all inputs
common_species <- intersect(rownames(colorPCs), rownames(morphoPCs))
common_species <- intersect(common_species, tree$tip.label)

# Prune tree and matrices
tree_pruned <- keep.tip(tree, common_species)
tree_pruned <- reorder.phylo(tree_pruned, order = "postorder")
colorPCs <- colorPCs[common_species, ]
morphoPCs <- morphoPCs[common_species, ]

# Ensure matching order
colorPCs <- colorPCs[tree_pruned$tip.label, ]
morphoPCs <- morphoPCs[tree_pruned$tip.label, ]

# Create correlation structure
cor_struct <- corBrownian(phy = tree_pruned, form = ~ species)
cor_struct <- corPagel(value = 0.5, phy = tree_pruned, fixed = FALSE, form = ~ species)

# Univariate
library(nlme)
library(MuMIn)
library(dplyr)

# Helper: fit one model with alternative correlation structures
fit_best_gls <- function() {
  fits <- list()
  
  # BM
  fits$BM <- try(gls(Color ~ Morpho, data = df,
                     correlation = corBrownian(phy = tree_pruned, form = ~ species),
                     method = "REML"), silent = TRUE)
  
  # Pagel’s λ with different starts
  mod <- lapply(seq(0.1, 1.5, by = 0.1), function(st) {
    m <- try(gls(Color ~ Morpho, data = df,
                   correlation = corPagel(value = st, phy = tree_pruned, fixed = FALSE, form = ~ species),
                   method = "REML"), silent = TRUE)
    if (inherits(m, "try-error")) {
      m <- NULL
    }
    m
  })
  mod <- mod[sapply(mod, length) > 0]
  mod <- mod[[which.max(sapply(mod, function(m) m$logLik))]]
  fits$Pagel <- mod

  # OU with different starts
  mod <- lapply(seq(0.1, 1.5, by = 0.1), function(st) {
    m <- try(gls(Color ~ Morpho, data = df,
                 correlation = corMartins(value = st, phy = tree_pruned, fixed = FALSE, form = ~ species),
                 method = "REML"), silent = TRUE)
    if (inherits(m, "try-error")) {
      m <- NULL
    }
    m
  })
  mod <- mod[sapply(mod, length) > 0]
  mod <- mod[[which.max(sapply(mod, function(m) m$logLik))]]
  fits$OU <- mod
  
  fits <- fits[sapply(fits, inherits, "gls")]
  if (length(fits) == 0L) return(NULL)
  
  # Compare models by AICc
  aictab <- data.frame(t(data.frame(sapply(fits, function(i) c(i$dims$p + 1,
                                                    MuMIn::AICc(i))), row.names = c("df", "AICc"))))
  aictab$Model <- names(fits)
  aictab <- aictab[order(aictab$AICc), , drop = FALSE]
  aictab$delta  <- aictab$AICc - min(aictab$AICc)
  aictab$weight <- exp(-0.5 * aictab$delta)
  aictab$weight <- aictab$weight / sum(aictab$weight)
  
  best_name <- rownames(aictab)[which.min(aictab$AICc)]
  best_model <- fits[[best_name]]
  
  return(list(best_model = best_model,
              best_name = best_name,
              aictab = aictab))
}

# =====================
# Univariate color~morpho across PCs
# =====================
results <- list()

for (i in 1:8) {
  for (j in 1:8) {
    cat("Fitting PGLS for colorPC", i , "and morphoPC", j, "...\n")
    df <- data.frame(
      Color = scale(colorPCs[, i]),
      Morpho = scale(morphoPCs[, j]),
      species = rownames(colorPCs)
    )
    
    fit <- fit_best_gls()
    
    if (is.null(fit)) {
      results[[paste0("ColorPC", i, "_MorphoPC", j)]] <- NA
    } else {
      summary_out <- summary(fit$best_model)
      slope <- coef(fit$best_model)["Morpho"]
      pval <- summary_out$tTable["Morpho", "p-value"]
      r2   <- cor(df$Color, predict(fit$best_model))^2
      
      results[[paste0("ColorPC", i, "_MorphoPC", j)]] <- list(
        slope = slope,
        p_value = pval,
        r_squared = r2,
        best_model = fit$best_name
      )
    }
  }
}

# Extract results to data frame
pgls_table <- do.call(rbind, lapply(names(results), function(name) {
  res <- results[[name]]
  if (is.null(res) || all(is.na(res))) return(NULL)
  
  data.frame(
    Comparison = name,
    colorPC = strsplit(name, "_")[[1]][1],
    morphoPC = strsplit(name, "_")[[1]][2],
    Slope = res$slope,
    P = res$p_value,
    R2 = res$r_squared,
    BestModel = res$best_model
  )
}))

# FDR correction within each colorPC
pgls_table <- pgls_table %>%
  group_by(colorPC) %>%
  mutate(P_adj = p.adjust(P, method = "fdr")) %>%
  dplyr::ungroup()

# Sort by p-value
pgls_table <- pgls_table[order(pgls_table$P_adj), ]
pgls_table_sig <- pgls_table[pgls_table$P_adj < 0.05,]

pgls_table_mat <- tidyr::spread(pgls_table %>%
                                  dplyr::select(colorPC, morphoPC, Slope), 
                                morphoPC, Slope) %>%
  as.data.frame()
rownames(pgls_table_mat) <- pgls_table_mat[,1]
pgls_table_mat <- pgls_table_mat[,-1]

pgls_table_mat_sig <- tidyr::spread(pgls_table %>%
                                      dplyr::select(colorPC, morphoPC, P_adj), 
                                    morphoPC, P_adj) %>%
  as.data.frame()
rownames(pgls_table_mat_sig) <- pgls_table_mat_sig[,1]
pgls_table_mat_sig <- pgls_table_mat_sig[,-1]

pgls_table_mat_sig <- ifelse(pgls_table_mat_sig < 0.05, ifelse(pgls_table_mat_sig < 0.01, ifelse(pgls_table_mat_sig < 0.001, "***", "**"), "*"), "")

zlims <- range(pretty(max(abs(pgls_table_mat)) * c(-1,1)))
pdf("./figures/col_morpho/heatmap_univariate_pairwise_pgls_colvsmorph.pdf", width = 8, height = 7)
pheatmap(pgls_table_mat,  angle_col = 45, legend_breaks = seq(zlims[1], zlims[2], length.out = 5),
         cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = pgls_table_mat_sig, 
         fontsize_number = 20, number_color = "black")
dev.off()

write.csv(pgls_table, "./results/col_morpho/color_morph_pgls.tex.csv")

out_file <- "results/col_morpho/color_morph_pgls.tex"

pgls_table <- pgls_table %>%
  mutate(
    colorPC_num  = as.integer(gsub("\\D+", "", colorPC)),
    morphoPC_num = as.integer(gsub("\\D+", "", morphoPC))
  ) %>%
  arrange(colorPC_num, morphoPC_num)
    
sig_idx <- which(pgls_table$P_adj < 0.05)

latex_tab <- pgls_table %>%
  mutate(
    colorPC = sub("Color", "", colorPC),
    morphoPC = sub("Morpho", "", morphoPC),
    Slope = sprintf("%.3f", Slope),
    P     = pvalue(P, accuracy = 0.0001),
    R2    = sprintf("%.3f", R2),
    `p(FDR)` = pvalue(P_adj, accuracy = 0.0001)  # rename for display
  ) %>%
  dplyr::select(colorPC, morphoPC, BestModel, R2, Slope, P, `p(FDR)`) %>% 
  kable(
    format   = "latex",
    booktabs = TRUE,
    align    = c("c","c","c","c","c","c","c"),
    caption  = "Phylogenetic generalized least squares (PGLS) associations between colour and morphological principal components (PCs). Shown are the PC pairing, slope, raw \\emph{p}-value, coefficient of determination ($R^2$), best-fitting correlation structure (Pagel or OU), and FDR-adjusted \\emph{p}-value.",
    col.names = c("Colour","Morphology","Best model","$R^2$","Slope","p","p(FDR)"),
    label    = "tab:color_morph_pgls",
    escape   = FALSE
  ) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  row_spec(c(0,sig_idx), bold = TRUE)

cat(latex_tab, file = out_file)

pgls_sig <- pgls_table %>%
  mutate(
    colorPC_num  = as.integer(gsub("\\D+", "", colorPC)),
    morphoPC_num = as.integer(gsub("\\D+", "", morphoPC)),
    sign         = ifelse(Slope > 0, "positive", "negative"),
    abs_slope    = abs(Slope)
  ) 

# Node positions ---------------------------------------------------------
n_color  <- length(unique(pgls_table$colorPC))
n_morpho <- length(unique(pgls_table$morphoPC))

color_nodes <- data.frame(
  name = sort(unique(pgls_table$colorPC)),
  x    = 0,
  y    = seq(from = n_color, to = 1, length.out = n_color),
  type = "Colour"
)

morpho_nodes <- data.frame(
  name = sort(unique(pgls_table$morphoPC)),
  x    = 1,
  y    = seq(from = n_morpho, to = 1, length.out = n_morpho),
  type = "Morphology"
)

nodes <- bind_rows(color_nodes, morpho_nodes)

# Helper to get node coordinates by PC name
get_node_coords <- function(pc_name) {
  nodes[nodes$name == pc_name, c("x", "y")]
}

# Edge table for significant links --------------------------------------
edges <- pgls_sig %>%
  rowwise() %>%
  mutate(
    x_color = get_node_coords(colorPC)$x,
    y_color = get_node_coords(colorPC)$y,
    x_morph = get_node_coords(morphoPC)$x,
    y_morph = get_node_coords(morphoPC)$y,
    signif_lab  = ifelse(P_adj < 0.05, "Significant", "Non-significant"),
    mid_x = (x_color + x_morph) / 2,
    mid_y = (y_color + y_morph) / 2 + ifelse(sign == "positive", 0.15, -0.15),
    beta_label = sprintf("%.3f", Slope),
    angle_raw = atan2(y_morph - y_color, x_morph - x_color) * 180 / pi,
    angle_label = round(ifelse(angle_raw > 90 | angle_raw < -90,
                         angle_raw + 180, 
                         ifelse(angle_raw < 0, angle_raw + 37,
                                ifelse(angle_raw > 45, angle_raw - 45,
                                       angle_raw - 37))))
  ) %>%
  ungroup()

# Simple two-colour palette for sign
sign_cols <- c(non_sig = "grey90",
               negative = "#D7301F",  # reddish
               positive = "#0571B0")  # bluish


p_bipartite <- ggplot() +
  # Non-significant edges
  geom_segment(
    data = edges %>% dplyr::filter(signif_lab == "Non-significant"),
    aes(x = x_color, xend = x_morph,
        y = y_color, yend = y_morph,
        colour = "non_sig"),
    size = 0.75,
    alpha = 0.8
  ) +
  # Significant edges
  geom_segment(
    data = edges %>% dplyr::filter(signif_lab == "Significant"),
    aes(x = x_color, xend = x_morph,
        y = y_color, yend = y_morph,
        colour = sign),
    size = 1.5,
    alpha = 0.8
  ) +
  geom_text(
    data = edges %>% dplyr::filter(signif_lab == "Significant"),
    aes(x = mid_x, y = mid_y, label = beta_label, 
        color = sign, angle = angle_label),
    size = 5, 
    fontface = "bold", 
    show.legend = FALSE 
  ) + 
  # Nodes
  geom_point(
    data = nodes,
    aes(x = x, y = y),
    shape = 21,
    size = 4,
    fill = "white",
    colour = "black",
    stroke = 1
  ) +
  geom_text(
    data = nodes,
    aes(x = x + ifelse(type == "Colour", -0.04, 0.04),
        y = y,
        label = name),
    hjust = ifelse(nodes$type == "Colour", 1, 0),
    size = 4.5
  ) +
  # Size = |beta|
  scale_size_continuous(
    name  = expression("|" * beta * "|"),
    range = c(0.2, 3)
  ) +
  
  # Colours = effect direction
  scale_colour_manual(
    values = sign_cols,
    breaks = c("positive", "non_sig", "negative"),
    labels = c("Positive", "Non-significant", "Negative")
  ) +
  coord_cartesian(
    xlim   = c(-0.2, 1.2),
    ylim   = c(min(nodes$y) - 0.5, max(nodes$y) + 0.5),
    expand = FALSE
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    axis.title  = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank(),
    legend.position = "top",
    legend.box      = "horizontal",
    legend.box.just = "left",
    plot.margin     = margin(0, 0, 0, 0),
    panel.spacing   = unit(0, "pt")
  ) +
  guides(
    # 3-row legend for effect direction
    colour = guide_legend(
      title = "Effect direction:  ",
      nrow  = 3,
      byrow = TRUE,
      override.aes = list(linewidth = 1.5),
      order = 1
    ),
    # 3-row legend for |beta|
    size = guide_legend(
      title = expression("|" * beta * "|"),
      nrow  = 3,
      byrow = TRUE,
      order = 2
    )
  )

p_bipartite

# Helper to read and turn a list of PC legend PNGs into a row strip -------
make_pc_strip <- function(pc_names, type = c("color", "morpho")) {
  type <- match.arg(type)
  
  p_list <- lapply(pc_names, function(pc) {
    pc_num <- gsub("\\D+", "", pc)
    file   <- sprintf("./figures/pca_legends/%sPC%s_legend.png",
                      ifelse(type == "color", "color", "morpho"),
                      pc_num)
    
    img <- png::readPNG(file)
    g   <- rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
    
    ggplot() +
      annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      coord_cartesian(expand = FALSE, clip = "off") +   # <- important
      theme_void() +
      theme(
        plot.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(0, "pt")
      )
  })
  
  # vertical strip where each panel stretches fully horizontally
  plot_grid(plotlist = p_list, ncol = 1, align = "v") +
    theme(plot.background = element_blank())
}

# Order PCs to match the node order we used -------------------------------
color_order  <- color_nodes$name
morpho_order <- morpho_nodes$name

color_strip  <- make_pc_strip(color_order,  type = "color")
morpho_strip <- make_pc_strip(morpho_order, type = "morpho")

# Final multi-panel figure: legends + network -----------------------------
combined <- plot_grid(
  color_strip,
  p_bipartite,
  morpho_strip,
  ncol = 3,
  rel_widths = c(0.9, 2.5, 0.9),   # tweak to taste
  align = "h",
  axis  = "tb"                     # align top & bottom
)

final_figure <- ggdraw(combined) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# Save
ggsave("./figures/col_morpho/fig_color_morpho_integration.png",
       final_figure, width = 15, height = 8, dpi = 300)

# =====================
# Multivariate color~morpho
# =====================
# Data prep 
Y <- as.matrix(colorPCs[tree_pruned$tip.label, 1:8])
X <- as.matrix(morphoPCs[tree_pruned$tip.label, 1:8])
dat <- list(Y = Y, X = X)

# --- fit candidates ---
# Use PL-LOOCV for stability with multivariate responses
fit_list <- list()

# BM
fit_list$BM_null    <- mvgls(Y ~ 1, data = dat, tree = tree_pruned, model = "BM",     method = "PL-LOOCV")
fit_list$BM_morpho  <- mvgls(Y ~ X, data = dat, tree = tree_pruned, model = "BM",     method = "PL-LOOCV")

# lambda (try a couple of starts implicitly handled by mvgls; it's robust for lambda)
fit_list$LAM_null   <- mvgls(Y ~ 1, data = dat, tree = tree_pruned, model = "lambda", method = "PL-LOOCV")
fit_list$LAM_morpho <- mvgls(Y ~ X, data = dat, tree = tree_pruned, model = "lambda", method = "PL-LOOCV")

# OU (set an 'upper' bound; escalate if needed)
ou_fit <- function(form) {
  u <- as.numeric(fit_list$LAM_null$start_values[2])
  max_attempts <- 10
  attempt <- 1
  repeat {
    model <- tryCatch(
      mvgls(form, data = dat, tree = tree_pruned, model = "OU", upper = u, method = "PL-LOOCV"),
      warning = function(w) "warning"
    )
    
    if (!is.character(model)) break  # Model fit successful, exit loop
    
    if (attempt >= max_attempts) {
      stop("Model failed after ", max_attempts, " attempts; last upper=", u)
    }
    
    cat("Warning caught: increasing 'upper' to", u, "\n")
    u <- u * 1.5  # Increase 'upper' by 50%
    attempt <- attempt + 1
  }
  return(model)
}

u <- as.numeric(fit_list$LAM_null$start_values[2])
max_attempts <- 10
attempt <- 1
repeat {
  model <- tryCatch(
    mvgls(Y ~ 1, data = dat, tree = tree_pruned, model = "OU", upper = u, method = "PL-LOOCV"),
    warning = function(w) "warning"
  )
  
  if (!is.character(model)) break  # Model fit successful, exit loop
  
  if (attempt >= max_attempts) {
    stop("Model failed after ", max_attempts, " attempts; last upper=", u)
  }
  
  cat("Warning caught: increasing 'upper' to", u, "\n")
  u <- u * 1.5  # Increase 'upper' by 50%
  attempt <- attempt + 1
}
fit_list$OU_null    <- model
attempt <- 1
repeat {
  model <- tryCatch(
    mvgls(Y ~ X, data = dat, tree = tree_pruned, model = "OU", upper = u, method = "PL-LOOCV"),
    warning = function(w) "warning"
  )
  
  if (!is.character(model)) break  # Model fit successful, exit loop
  
  if (attempt >= max_attempts) {
    stop("Model failed after ", max_attempts, " attempts; last upper=", u)
  }
  
  cat("Warning caught: increasing 'upper' to", u, "\n")
  u <- u * 1.5  # Increase 'upper' by 50%
  attempt <- attempt + 1
}
fit_list$OU_morpho  <- model

# EB
fit_list$EB_null    <- mvgls(Y ~ 1, data = dat, tree = tree_pruned, model = "EB", method = "PL-LOOCV")
fit_list$EB_morpho  <- mvgls(Y ~ X, data = dat, tree = tree_pruned, model = "EB", method = "PL-LOOCV")

eics <- lapply(fit_list, EIC, nbcores = ncores)

# --- compare by EIC ---
model_names <- names(fit_list)
tab  <- data.frame(
  model = names(eics),
  logLik = sapply(eics, `[[`, "LogLikelihood"),
  k      = sapply(eics, `[[`, "p"),
  se     = sapply(eics, `[[`, "se"),
  EIC    = sapply(eics, `[[`, "EIC"),
  row.names = NULL
)
tab$delta  <- tab$EIC - min(tab$EIC, na.rm = TRUE)
w          <- exp(-0.5 * tab$delta)
tab$weight <- w / sum(w, na.rm = TRUE)
tab$param <- ifelse(grepl("^LAM", model_names), "lambda", 
                    ifelse(grepl("^OU", model_names), "alpha", 
                           ifelse(grepl("^EB", model_names), "a", NA)))
tab$param_value <- sapply(fit_list, `[[`, "param")

# select model
idx <- which.min(tab$EIC)
best_name <- tab$model[idx]
best_model <- fit_list[[best_name]]

# permutation MANOVA only on the winning model
man <- manova.gls(best_model, nperm = 999, test = "Pillai", verbose = TRUE)

tab$pillai[tab$model == best_name] <-  man$stat
tab$pvalue[tab$model == best_name] <-  man$pvalue

# Write to csv
write.csv(tab, "./results/col_morpho/mv_col_morpho.csv")
