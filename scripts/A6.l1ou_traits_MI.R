#=============================#
# Libraries
#=============================#
lib <- c("ape", "picante", "sp", "geiger", "wesanderson",
         "hypervolume", "car", "bayou", "nlme", "l1ou", 
         "parallel", "genlasso", "mvMORPH", "pheatmap",
         "lme4", "mgcv", "phytools", "corHMM",
         "lmtest", "ColorAR", "castor", "mvMORPH",
         "RPANDA","OUwie", "ggplot2", "ggpubr", "gridExtra")
sapply(lib, library, character.only = T)

#=============================#
# Settings
#=============================#
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/")

col.palette <- sample(unique(as.character(sapply(c("Zissou1", "Darjeeling1", "Darjeeling2", "FantasticFox1"), function(x) wes_palette(x, 5)))))
ncores <- 5

# load input data and custom functions
source("./Rscripts/A0.Input_data.R")

dir.create("new_figures/l1ou")

# load Ancestral trait reconstructions
for (trait in discrete_traits){
  if (file.exists(paste0("./rdata/ASR_", trait, ".rda"))){
    ## load trait ASR
    load(paste0("./rdata/ASR_", trait, ".rda"))
    
    assign(paste0("SCM_", trait), trait_SCM_trees)
    
    assign(paste0("simmap_", trait), make.consensus.simmap(trait_SCM_trees))
  }  
}

# Testing if trait convergence is explained by ecological traits
################################################################
## Look at data
selected_pcs <- 1:8
convergent.only <- TRUE
### Color PCA, regimes and traits
for (trait in discrete_traits){
  
  trait.vec = setNames(eco_traits_subset[, trait], rownames(eco_traits_subset))[rownames(colorPCs)]
  
  trait.cols <- setNames(col.palette[1:length(levels(trait.vec))], levels(trait.vec))
  legend_text <- names(trait.cols)
  
  # plot color pc vs trait (raw)
  pdf(paste0("./new_figures/", trait, ".vs.colPC", pc, ".pdf"), width = 10, height = 12)
  layout(cbind(matrix(1:8, nrow = 4, ncol = 2, byrow = TRUE), 9), widths = c(1,1,0.4))
  par(mar = c(2.1,5.1,2.1,2.1))
  lapply(selected_pcs, function(i) {
    plot(colorPCs[,i]~trait.vec, pch = 21, col = trait.cols, 
         bty = "n", xlab = "", xaxt = "n", ylab = paste("color PC", i), 
         cex = 1.2, cex.axis = 1.2, cex.lab = 1.5)
  })
  par(mar = c(0,0,0,0))
  plot(NULL, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "", bty = "n")
  legend("left", legend_text, bty = "n", pch = 21, pt.cex = 1.5, 
         pt.bg = trait.cols, cex = 1.2)
  dev.off()
  
  # plot morpho pc vs trait (raw)
  pdf(paste0("./new_figures/", trait, ".vs.morphoPC", pc, ".pdf"), width = 10, height = 12)
  layout(cbind(matrix(1:8, nrow = 4, ncol = 2, byrow = TRUE), 9), widths = c(1,1,0.4))
  par(mar = c(2.1,5.1,2.1,2.1))
  lapply(selected_pcs, function(i) {
    plot(morphoPCs[,i]~trait.vec, pch = 21, col = trait.cols, 
         bty = "n", xlab = "", xaxt = "n", ylab = paste("morpho PC", i), 
         cex = 1.2, cex.axis = 1.2, cex.lab = 1.5)
  })
  par(mar = c(0,0,0,0))
  plot(NULL, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "", bty = "n")
  legend("left", legend_text, bty = "n", pch = 21, pt.cex = 1.5, 
         pt.bg = trait.cols, cex = 1.2)
  dev.off()
  
  
  for (pc in selected_pcs){
    
    if (file.exists(paste0("./rdata/colPC_", pc, "_damsel_l1ou_estimated_conv_shifts_pBIC_new.rds"))){
      
      ## COLOR
      
      ## load l1ou model
      l1ou_model <- readRDS(paste0("./rdata/colPC_", pc, "_damsel_l1ou_estimated_conv_shifts_pBIC_new.rds"))
      
      PCdata <- setNames(colorPCs[,pc], rownames(colorPCs))
      
      # extract l1ou shifts and plot versus trait
      shift.nodes <- l1ou_model$shift.configuration
      
      regimes <- names(shift.nodes)
      
      if (convergent.only) {
        is_convergent <- duplicated(regimes) | duplicated(regimes, fromLast = TRUE)
        shift.nodes <- shift.nodes[is_convergent] 
        regimes <- regimes[is_convergent]
      } else {
        is_convergent <- rep(TRUE, length(shift.nodes))
      }
      
      if (length(regimes) > 0) {
        # Mark nodes with shift names
        for (n in seq_along(shift.nodes)) {
          l1ou_model$tree$node.label[shift.nodes[[n]]] <- names(shift.nodes)[n]
        }
        
        # Identify tips under each shift
        shift.tips <- lapply(shift.nodes, function(nodes) {
          unique(unlist(sapply(nodes, function(n) {
            if (n > Ntip(l1ou_model$tree)) tips(l1ou_model$tree, n)
            else l1ou_model$tree$tip.label[n]
          })))
        })
        
        # Create a labeling vector
        shift.vec <- setNames(rep("Background", Ntip(l1ou_model$tree)), l1ou_model$tree$tip.label)
        for (i in seq_along(shift.tips)) {
          shift.vec[shift.tips[[i]]] <- paste("Regime", names(shift.tips)[i])
        }
        
        # Prepare plotting data
        data <- data.frame(
          species = names(PCdata),
          PC = PCdata,
          Trait = trait.vec,
          Regime = factor(shift.vec[names(PCdata)])
        )
        
        # Plot
        p <- ggplot(data, aes(x = Trait, y = PC, fill = Regime)) +
          geom_boxplot(outlier.size = 2, linewidth = 1.2, alpha = 0.8, position = position_dodge(width = 0.8)) +
          scale_fill_manual(values = c("Background" = "grey80", 
                                       setNames(colorRampPalette(c("#FEC194", "#FA824C"))(length(unique(data$Regime))-1),
                                                levels(data$Regime)[levels(data$Regime) != "Background"]))) +
          labs(x = "", y = paste0("Color PC", pc)) +
          facet_wrap(~ Regime, nrow = 1) +
          theme_light(base_size = 16) +
          coord_flip() +
          theme(
            plot.margin = margin(20, 20, 20, 20),
            axis.title.x = element_text(vjust = -1),
            axis.title.y = element_text(vjust = +1),
            legend.position = "none",
            strip.text = element_text(size = 14, face = "bold"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey85")
          )
  
        ggsave(filename = paste0("./new_figures/l1ou/colorPCs", pc, "_Trait_", trait, "_L1OU_regimes.pdf"), plot = p, width = 12, height = 6)
        
        ## MORPHO
        
        ## load l1ou model
        l1ou_model <- readRDS(paste0("./rdata/morphoPC_", pc, "_damsel_l1ou_estimated_conv_shifts_pBIC_new.rds"))
        
        # extract l1ou shifts and plot versus trait
        shift.nodes <- l1ou_model$shift.configuration
        
        regimes <- names(shift.nodes)
        
        if (convergent.only) {
          is_convergent <- duplicated(regimes) | duplicated(regimes, fromLast = TRUE)
          shift.nodes <- shift.nodes[is_convergent] 
          regimes <- regimes[is_convergent]
        } else {
          is_convergent <- rep(TRUE, length(shift.nodes))
        }
        
        if (length(regimes) > 0) {
          # Mark nodes with shift names
          for (n in seq_along(shift.nodes)) {
            l1ou_model$tree$node.label[shift.nodes[[n]]] <- names(shift.nodes)[n]
          }
          
          # Identify tips under each shift
          shift.tips <- lapply(shift.nodes, function(nodes) {
            unique(unlist(sapply(nodes, function(n) {
              if (n > Ntip(l1ou_model$tree)) tips(l1ou_model$tree, n)
              else l1ou_model$tree$tip.label[n]
            })))
          })
          
          # Create a labeling vector
          shift.vec <- setNames(rep("Background", Ntip(l1ou_model$tree)), l1ou_model$tree$tip.label)
          for (i in seq_along(shift.tips)) {
            shift.vec[shift.tips[[i]]] <- paste("Regime", names(shift.tips)[i])
          }
          
          # Prepare plotting data
          data <- data.frame(
            species = names(PCdata),
            PC = PCdata,
            Trait = trait.vec,
            Regime = factor(shift.vec[names(PCdata)])
          )
          
          # Plot
          p <- ggplot(data, aes(x = Trait, y = PC, fill = Regime)) +
            geom_boxplot(outlier.size = 2, linewidth = 1.2, alpha = 0.8, position = position_dodge(width = 0.8)) +
            scale_fill_manual(values = c("Background" = "grey80", 
                                         setNames(colorRampPalette(c("#FEC194", "#FA824C"))(length(unique(data$Regime))-1),
                                                  levels(data$Regime)[levels(data$Regime) != "Background"]))) +
            labs(x = "", y = paste0("Color PC", pc)) +
            facet_wrap(~ Regime, nrow = 1) +
            theme_light(base_size = 16) +
            coord_flip() +
            theme(
              plot.margin = margin(20, 20, 20, 20),
              axis.title.x = element_text(vjust = -1),
              axis.title.y = element_text(vjust = +1),
              legend.position = "none",
              strip.text = element_text(size = 14, face = "bold"),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_line(color = "grey85")
            )
          
          print(p)
          ggsave(filename = paste0("./new_figures/l1ou/morphoPCs", pc, "_Trait_", trait, "_L1OU_regiomes.pdf"), plot = p, width = 12, height = 6)
        }
      }
    }
  }
}

#################################################
## Color
########
convergent_pcs <- 1:8
convergent_pcs <- c(1,2,6)
criterion <- "pBIC"
convergent.only <- FALSE

# List of ecological traits to compare with col PCs
traits <- list(diet = SCM_DietEcotype, habitat = SCM_habitat1, symbiosis = SCM_Farming)
col_MI_result <- list()
col_p_values_matrix <- col_MI_matrix <- matrix(NA, nrow = length(convergent_pcs), ncol = length(traits), dimnames = list(paste0("PC", convergent_pcs), names(traits)))

# Loop through each convergent PC
for (i in 1:length(convergent_pcs)) {
  
  pc <- convergent_pcs[i]
  
  cat("Analyzing color PC", pc, "...\n")
  
  # Read the model results for the convergent PC
  fit <- readRDS(paste0("./rdata/colPC_", pc, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds"))
  
  # Generate mapped tree for the PC
  mapped_tree <- l1ou2simmap(fit, convergent.only = convergent.only, plot = FALSE, randomize = FALSE)

  # Loop through each ecological trait to test MI
  for (trait_name in names(traits)) {
    trait_states <- traits[[trait_name]]  

    # Run the MI test
    MI <- MI.test(trait_states, mapped_tree, nsim = 100, plot = FALSE, normalize = TRUE) 

    # Store results
    col_MI_result[paste0("PC", pc, "-", trait_name)] <- list(MI)

    # Save mean MI value in the matrix
    col_MI_matrix[i, which(names(traits) == trait_name)] <- mean(MI$obs)
    
    # Save p-value in the matrix
    col_p_values_matrix[i, which(names(traits) == trait_name)] <- mean(MI$p.val)
    
    # Print the results for the current PC
    cat("Results for", trait_name, ":\n")
    print(setNames(c(mean(MI$obs), mean(MI$p.val)), c("MI", "p-value")))
  }
}
saveRDS(col_MI_result, "./rdata/l1ou_allShifts_MI_traits.rds")
# plot
pdf("./new_figures/l1ou/colorPCs_traits_MI_all.pdf", width = 15, height = 25)
par(mfrow = c(length(convergent_pcs), length(traits)))
lapply(names(col_MI_result), function(mi) {
  title <- strsplit(mi, "-")[[1]]
  plot.MI(col_MI_result[[mi]], cex = c(0.8, 1), title = paste0(title[2], " and color ", title[1], "\nMutual Information"))
  })
dev.off()

# Plot the heatmap
heatmap_data_asterisks <- col_p_values_matrix
heatmap_data_asterisks[heatmap_data_asterisks < 0.05] <- "*"
heatmap_data_asterisks[heatmap_data_asterisks >= 0.05] <- ""
heatmap_data_asterisks <- matrix(paste(round(col_MI_matrix, 2), heatmap_data_asterisks), nrow = nrow(col_MI_matrix), ncol = ncol(col_MI_matrix))

rownames(col_MI_matrix) <- paste("color", rownames(col_MI_matrix))
# Plot
pdf(paste0("./new_figures/l1ou/colorPCs_traits_MI_heatmap_all.pdf"), width = 15, height = 6)
pheatmap(
  t(col_MI_matrix),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "#FFD6E7", "#FF99CC", "#FF3399", "#CC0066"))(100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  display_numbers = t(heatmap_data_asterisks),
  number_color = "black",           # Numbers in black for contrast
  fontsize_number = 14,             # Bigger numbers
  fontsize_row = 10,                # Bigger row labels
  fontsize_col = 10,                # Bigger col labels
  angle_col = 0,                    # Horizontal column names
  main = "Mutual Information\nbetween Color Pattern and Ecological Traits",
  legend = TRUE,
  cex = 1.2,
  legend_breaks = seq(0, max(pretty(col_MI_matrix, na.rm = TRUE)), length.out = 5),
  legend_labels = round(seq(0, max(pretty(col_MI_matrix)), length.out = 5), 2)
)
dev.off()

## Morpho
##########
convergent_pcs <- 1:8
convergent_pcs <- c(1,3,4,7)
criterion <- "pBIC"
convergent.only <- TRUE

# List of ecological traits to compare with morpho PCs
traits <- list(diet = SCM_DietEcotype, habitat = SCM_habitat1, symbiosis = SCM_Farming)
morpho_MI_result <- list()
morpho_p_values_matrix <- morpho_MI_matrix <- matrix(NA, nrow = length(convergent_pcs), ncol = length(traits), dimnames = list(paste0("PC", convergent_pcs), names(traits)))

# Loop through each convergent PC
for (i in 1:length(convergent_pcs)) {
  
  pc <- convergent_pcs[i]
  
  cat("Analyzing morpho PC", pc, "...\n")
  
  # Read the model results for the convergent PC
  fit <- readRDS(paste0("./rdata/morphoPC_", pc, "_damsel_l1ou_estimated_conv_shifts_", criterion, ".rds"))
  
  # Generate mapped tree for the PC
  mapped_tree <- l1ou2simmap(fit, convergent.only = convergent.only, plot = FALSE, randomize = FALSE)

  # Loop through each ecological trait to test MI
  for (trait_name in names(traits)) {
    trait_states <- traits[[trait_name]]  
    
    # Run the MI test
    MI <- MI.test(trait_states, mapped_tree, nsim = 100, plot = FALSE, normalize = TRUE) 
    
    # Store results
    morpho_MI_result[paste0("PC", pc, "-", trait_name)] <- list(MI)
    
    # Save mean MI value in the matrix
    morpho_MI_matrix[i, which(names(traits) == trait_name)] <- mean(MI$obs)
    
    # Save p-value in the matrix
    morpho_p_values_matrix[i, which(names(traits) == trait_name)] <- mean(MI$p.val)
    
    # Print the results for the current PC
    cat("Results for", trait_name, ":\n")
    print(setNames(c(mean(MI$obs), mean(MI$p.val)), c("MI", "p-value")))
  }
}
saveRDS(morpho_MI_result, "./rdata/l1ou_convShifts_MI_traits.rds")

# plot
pdf("./new_figures/l1ou/morphoPCs_traits_MI.pdf", width = 15, height = 15)
par(mfrow = c(length(convergent_pcs),length(traits)))
lapply(names(morpho_MI_result), function(mi) {
  title <- strsplit(mi, "-")[[1]]
  plot.MI(morpho_MI_result[[mi]], cex = c(0.8, 1), title = paste0(title[2], " and morpho ", title[1], "\nMutual Information"))
})
dev.off()

# Plot the heatmap
heatmap_data_asterisks <- morpho_p_values_matrix
heatmap_data_asterisks[heatmap_data_asterisks < 0.05] <- "*"
heatmap_data_asterisks[heatmap_data_asterisks >= 0.05] <- ""
heatmap_data_asterisks <- matrix(paste(round(morpho_MI_matrix, 2), heatmap_data_asterisks), 
                                 nrow = nrow(morpho_MI_matrix), ncol = ncol(morpho_MI_matrix))

rownames(morpho_MI_matrix) <- paste("morpho", rownames(morpho_MI_matrix))
# Plot
pdf(paste0("./new_figures/l1ou/morphoPCs_traits_MI_heatmap.pdf"), width = 10, height = 6)
par(mar = c(8,8,8,8))
pheatmap(
  t(morpho_MI_matrix),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "#FFE5B4", "#FFC266", "#FF9900", "#CC6600"))(100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  display_numbers = t(heatmap_data_asterisks),
  number_color = "black",           # Numbers in black for contrast
  fontsize_number = 14,             # Bigger numbers
  fontsize_row = 10,                # Bigger row labels
  fontsize_col = 10,                # Bigger col labels
  angle_col = 0,                    # Horizontal column names
  main = "Mutual Information\nbetween Morphology and Ecological Traits",
  legend = TRUE,
  cex = 1.2,
  legend_breaks = seq(0, max(pretty(morpho_MI_matrix, na.rm = TRUE)), length.out = 5),
  legend_labels = round(seq(0, max(pretty(morpho_MI_matrix, na.rm = TRUE)), length.out = 5), 2)
)
dev.off()

