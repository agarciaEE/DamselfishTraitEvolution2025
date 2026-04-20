################################################################################
################################################################################
####----------------------------------------------------------------------######
####----------------------------------------------------------------------######
###             Damselfish ancestral state reconstructions                ######
####----------------------------------------------------------------------######
####----------------------------------------------------------------------######
################################################################################
################################################################################

library(phytools)
library(corHMM)
library(castor)
library(lmtest)
library(ColorAR)
library(wesanderson)

setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

# load input data and custom functions
source("./scripts/A1.Input_data.R")

# create output directory
dir.create("./results/AncRec")

# settings
width = 15 # figure width
height = 15 # figure height
nsim = 100
HSM = FALSE # Don't run Hidden States Models for now
mcmc = FALSE # Run SCM on fitted models
ncores = 5

################################################################################
####                                                                      ######
#### --------------   Reconstructing ancestral states    ---------------- ######
####                                                                      ######
################################################################################

################################################################################
###  DISCRETE TRAIT  ####
################################################################################

tree_rescaled <- geiger::rescaleTree(tree, 1.0)
k <- max(ape::branching.times(tree))

# Specific transition matrices 
# ============================

## DietEcotype
diet_states <- c("Benthic","Intermediate","Pelagic")
Q_diet_conservative <- matrix(0, 3, 3, dimnames = list(diet_states, diet_states))
Q_diet_conservative["Benthic","Intermediate"]    <- 1  # a1
Q_diet_conservative["Intermediate","Benthic"]    <- 2  # a2
Q_diet_conservative["Intermediate","Pelagic"]    <- 3  # b1
Q_diet_conservative["Pelagic","Intermediate"]    <- 4  # b2
# B<->P remain 0 (forbidden)

Q_diet_relaxed <- Q_diet_conservative
Q_diet_relaxed["Benthic","Pelagic"] <- 5  # c
Q_diet_relaxed["Pelagic","Benthic"] <- 5  # c (symmetric rare jump)

## Symbiosis
sym_states <- c("Free-living","Comensalistic","Mutualistic")
Q_sym_conservative <- matrix(0, 3, 3, dimnames = list(sym_states, sym_states))
Q_sym_conservative["Free-living","Comensalistic"]  <- 1  # g1
Q_sym_conservative["Comensalistic","Mutualistic"]  <- 2  # g2
Q_sym_conservative["Mutualistic","Comensalistic"]  <- 3  # l2
Q_sym_conservative["Comensalistic","Free-living"]  <- 4  # l1
# Free<->Mutualistic remain 0

Q_sym_relaxed <- Q_sym_conservative
Q_sym_relaxed["Free-living","Mutualistic"] <- 5  # r
Q_sym_relaxed["Mutualistic","Free-living"] <- 5  # r

## Habitat
hab_states <- c("sea anemone","freshwater","non-reef","rocky-reef","coral-reef")
Q_hab_conservative <- matrix(0, 5, 5, dimnames = list(hab_states, hab_states))

# Coral <-> Anemone
Q_hab_conservative["coral-reef","sea anemone"] <- 1  # m (gain specialization)

# Marine stepwise
Q_hab_conservative["non-reef","rocky-reef"]    <- 2  # n
Q_hab_conservative["rocky-reef","non-reef"]    <- 3  # p
Q_hab_conservative["rocky-reef","coral-reef"]  <- 4  # q
Q_hab_conservative["coral-reef","rocky-reef"]  <- 5  # r

# Freshwater isolated via non-reef only
Q_hab_conservative["non-reef","freshwater"]    <- 6  # s

# All other off-diagonals remain 0 (forbidden)

Q_hab_relaxed <- Q_hab_conservative

# Rare non-reef <-> coral-reef jump
Q_hab_relaxed["non-reef","coral-reef"] <- 7  # t
Q_hab_relaxed["coral-reef","non-reef"] <- 8  # t

# Rocky-reef -> rocky-reef
Q_hab_relaxed["rocky-reef","freshwater"] <- 9 # w

# Rocky-reef -> Sea anemone
Q_hab_relaxed["rocky-reef","sea anemone"] <- 10 # x

# Non-reef -> Sea anemone
Q_hab_relaxed["reef-reef","sea anemone"] <- 11 # y

Q_List<- list(DietEcotype = list(CONSERVATIVE = Q_diet_conservative,
                                 RELAXED = Q_diet_relaxed),
              Symbiosis = list(CONSERVATIVE = Q_sym_conservative,
                               RELAXED = Q_sym_relaxed),
              Habitat = list(CONSERVATIVE = Q_hab_conservative,
                             RELAXED = Q_hab_relaxed))

for (trait in discrete_traits) {
  
  message("Working on ", trait, " trait...")
  trait.vec <- setNames(eco_traits_subset[, trait], rownames(eco_traits_subset))
  # make sure the order of the trait matches the tree
  trait.vec <- trait.vec[tree$tip.label]
  states <- levels(trait.vec)
  nstates <- length(states)
  trait_title = ifelse(trait == "DietEcotype", "Diet ecotype", trait)
  
  # sample color palette
  c_i <- runif(1)
  if (c_i < 0.33){
    cols<-setNames(viridis::viridis(nstates), states)
  } else if (c_i < 0.66) {
    cols<-setNames(as.character(wes_palette(sample(c("Zissou1", "Darjeeling1", 
                                                     "Darjeeling2", "FantasticFox1"), 1), 
                                            nstates)), states)
  } else {
    cols<-setNames(ggsci::pal_jco("default")(nstates), states)
  }
  barplot(1:length(cols), col = cols, names.arg = states)

  # Description of some parameters
  cat("\nNumber of tips:", length(trait.vec))
  cat("\nNumber of states:", nstates)
  cat("\nStates:", states)
  cat("\n")

  message("Plotting phylogeny with trait on tips...")
  pdf(paste0("./results/AncRec/phylo_", trait, ".pdf"), width = width, height = height)
  # plot trait on tips
  plotTree(tree,show.tip.label = TRUE, type = "fan", part = 0.88,
           fsize = 0.5, offset = 10, lwd = 3)
  xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(tree)]
  yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(tree)]
  # add ages lines
  obj<-axis(1,pos=-2,at=seq(10,50,by=10),cex.axis=0.5,labels=FALSE)
  text(obj,rep(-5,length(obj)), rev(obj),cex=0.6)
  text(mean(obj),-8,"time (ma)",cex=0.8)
  for(i in 1:(length(obj)-1)){
    a1<-atan(-2/obj[i])
    a2<-0.88*2*pi
    plotrix::draw.arc(0,0,radius=obj[i],a1,a2,lwd=1.5,
             col=make.transparent("grey40",0.25))
  }
  ## add traits
  points(xx*1.025, yy*1.025, pch = 21, bg = cols[trait.vec[tree$tip.label]])
  # add legend
  legend(x="topleft", legend=levels(trait.vec),
         pt.cex=2, pch=16, col=cols, title = trait_title, bty = "n")
  dev.off()
  
  ################
  ## Fit Mk models
  ################
  message("Carrying out Markov k-state (Mk) Models...")
  Mk_models <- c("ER", "ARD", "SYM", "DIR1", "DIR2", "DIR3", "CONSERVATIVE", "RELAXED")

  trait_Mk_fit <- lapply(Mk_models, function(i) {
    
    if (i %in% c("CONSERVATIVE", "RELAXED")) {
      Q <- Q_List[[trait]][[i]]
    } else {
      # make Q matrix
      Q <- make_Q_matrix(i, nstates, states)
      # forbid transitions from freshwater and sea anemone as they are one directional specializations (only transitions to allowed)
      if (trait == "Habitat") {
        Q["freshwater","sea anemone"] = 0
        Q["sea anemone","freshwater"] = 0
      }
    }
    print(Q)
    # fit Mk
    message("Fitting ", i, " model...")
    fitMk(tree_rescaled, trait.vec, model = Q, pi = "fitzjohn")
  })
  names(trait_Mk_fit) <- Mk_models

  message("Plotting Mk Models transitions...")
  pdf(paste0("./results/AncRec/", trait, "_Mk.fitted_Q.transitions.pdf"), width = width+5, height = height+5)
  par(mfrow = c(ceiling(sqrt(length(Mk_models))), ceiling(sqrt(length(Mk_models)))))
  lapply(names(trait_Mk_fit), function(i) plot(trait_Mk_fit[[i]], 
                                               main  = paste0(i, "\nLogLik = ", round(trait_Mk_fit[[i]]$logLik,3),
                                                             ";k = ", attr(logLik(trait_Mk_fit[[i]]), "df")),
                                               offset = 0.05, lwd = 2, show.zeros = F, spacer = 0.25,
                                               cex.main = 2, cex.traits = 2, cex.rates = 1.5))
  dev.off()
  
  aic_vec <- sapply(trait_Mk_fit, AIC)
  k_vec   <- sapply(trait_Mk_fit, function(m) attr(logLik(m),"df"))
  logL    <- sapply(trait_Mk_fit, logLik)
  n_obs   <- length(trait_vec)
  AICc    <- aic_vec + (2*k_vec*(k_vec+1))/(n_obs - k_vec - 1)
  wts     <- geiger::aicw(AICc)$w
  
  Mk_res_table <- tibble::tibble(model = Mk_models, logLik = as.numeric(logL), k = as.numeric(k_vec), AIC = as.numeric(aic_vec),
                                 AICc = as.numeric(AICc)) %>%
    mutate(delta.AICc = AICc - min(AICc, na.rm=TRUE), weight = wts) %>% arrange(AICc)
  
  if (HSM) {
    
    ####################
    ## Fit Hidden models
    ####################
    message("Carrying out Hidden State Models...")
    
    # Hidden Rates models:
    ## nstates models in which only one state has an extra hidden state
    ## nstates * (nstates - 1) / 2 models in which two states have an extra hidden state
    ## 1 model with all states having and extra hidden state
    ### !!! Combinatorial hidden states models may be too computationally complex and not necessary, specially when number of states > 3.
    
    # get the best ER, SYM or ARD model
    best_simple <- Mk_res_table %>%
      filter(model %in% c("ER", "SYM", "ARD")) %>%
      arrange(AIC) %>%
      slice(1)
    
    nhrm <- nstates + nstates*(nstates-1)/2 + 1
    HRM_matrix <- matrix(1, nrow = nstates, ncol = nhrm, 
                         dimnames = list(states, paste0("H", 1:nhrm)))
    
    hrm <- c(lapply(1:nstates, function(i) i), combn(nstates, 2, simplify = F), list(c(1:nstates)))
    for (i in 1:length(hrm)) {
      HRM_matrix[hrm[[i]], i] = 2
    }
    HRM_matrix <- t(unique(t(HRM_matrix)))
    nhrm <- ncol(HRM_matrix)
    HR_models <-  paste0("H", 1:nhrm)
    i
    # using phytools:
    trait_HR.fit <- lapply(1:length(HR_models), function(i) {
      cat("Fitting Hidden Rates Model with", HRM_matrix[,i])
      fitHRM(tree_rescaled, trait.vec, 
             model = best_simple$model, ncores = ncores, 
             ncat = HRM_matrix[,i], pi="fitzjohn")
      })

    names(trait_HR.fit) <- HR_models
    
    message("Plotting Hidden State Models transitions...")
    pdf(paste0("./results/AncRec/", trait, "_HR.fitted_Q.transitions.pdf"), width = width-5, height = height-5)
    par(mfrow = c(2, 3))
    lapply(names(trait_HR.fit), function(i) plot(as.Qmatrix(trait_HR.fit[[i]]), main  = i, offset = 0.1, spacer = 0.25, lwd = 2, show.zeros = T,
                                                 cex.main = 2, cex.traits = 1, cex.rates = 1))
    dev.off()
    
    HR_res_table <- tibble::tibble(model = HR_models,
                                   logLik = sapply(HR_fit, function(i) i$loglik),
                                   k      = sapply(HR_fit, function(i) i$DF),
                                   AIC    = sapply(HR_fit, function(i) i$AIC)) %>%
      mutate(AICc = AIC + (2*k*(k+1))/(n_obs - k - 1),
             delta.AICc = AICc - min(AICc, na.rm=TRUE),
             weight = geiger::aicw(AICc)$w) %>% arrange(AICc)
    
  } else {
    trait_HR.fit <- NULL
    HR_models <- NULL
    HR.res.table <- NULL
  }
  
  ####
  # build a table with the results of all Mk and HR models
  ####
  trait.discrete.fit <- c(trait_Mk_fit, trait_HR.fit)
  trait.discrete.models <- c(Mk_models, HR_models)
  
  res.table <- rbind(Mk_res_table, HR.res.table)

  # write table to hard drive
  message("Writting discrete models result table on hard disk...")
  write.table(res.table, file = paste0("./results/AncRec/", trait, "_models.fit.restable.txt"), 
                                       quote = FALSE, sep = "\t", row.names = T, col.names = T)
  
  # select best discrete model
  best_discrete_model <- res.table$model[which.min(res.table$AIC)]

  if (best_discrete_model %in% Mk_models) { 
    best_Qmat <- trait.discrete.fit[[best_discrete_model]]$index.matrix
    best_Qmat[is.na(best_Qmat)] = 0
    rate.cat <-  1
  } else {
    best_Qmat <- trait.discrete.fit[[best_discrete_model]]$index.mat
    best_Qmat[is.na(best_Qmat)] = 0
  }
  
  best_Qrates <- best_Qmat
  for (i in unique(as.vector(best_Qrates))) {
    if (i == 0) next
    rate = trait.discrete.fit[[best_discrete_model]]$rates[i]
    # back-transformed scale
    rate <- rate / k
    best_Qrates[best_Qrates == i] = rate
  }
  
  #######################################
  # Joint ancestral state reconstructions
  #######################################

  fit.joint <- corHMM(tree, 
                      data.frame(species = names(trait.vec), 
                                 trait = trait.vec), 
                      node.states = "joint", 
                      rate.mat = best_Qmat,
                      rate.cat = rate.cat, 
                      root.p = "fitzjohn", 
                      nstarts = 10, n.cores = 5, 
                      get.tip.states = FALSE)

  # plot
  message("Plotting Joint ancestral state reconstruction...")
  pdf(paste0("./results/AncRec/Joint.ASR_model.fit_", best_discrete_model, "_", trait, ".pdf"), width = width, height = height)
  plotTree(tree,show.tip.label = TRUE, type = "fan", part = 0.88,
           fsize = 0.5, offset = 10, lwd = 3)
  xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(tree)]
  yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(tree)]
  # add ages lines
  obj<-axis(1,pos=-2,at=seq(10,50,by=10),cex.axis=0.5,labels=FALSE)
  text(obj,rep(-5,length(obj)), rev(obj),cex=0.9)
  text(mean(obj),-8,"time (ma)",cex=0.9)
  for(i in 1:(length(obj)-1)){
    a1<-atan(-2/obj[i])
    a2<-0.88*2*pi
    plotrix::draw.arc(0,0,radius=obj[i],a1,a2,lwd=1.5,
                      col=make.transparent("grey40",0.25))
  }
  ## add traits
  points(xx*1.025, yy*1.025, pch = 21, bg = cols[trait.vec[tree$tip.label]])
  # add legend
  colLegend = unique(data.frame(id = as.numeric(fit.joint$data.legend[,2]), 
                                trait = as.character(trait.vec[fit.joint$data.legend[,1]])))
  colLegend <- sort(setNames(colLegend[,1], colLegend[,2]))
  legend(x="topleft", legend=names(colLegend),
         pt.cex=2, pch=16, col=cols[names(colLegend)], title = trait_title, bty = "n")
  # add states at nodes
  nodelabels(pie = to.matrix(names(colLegend)[fit.joint$phy$node.label],
                             names(colLegend)), piecol = cols[names(colLegend)], cex = 0.2)
  dev.off()

  ##########################################
  # Marginal ancestral state reconstructions
  ##########################################
  
  fit.marginal <- corHMM(tree, 
                         data.frame(species = names(trait.vec), 
                                    trait = trait.vec), 
                      node.states = "marginal", 
                      rate.mat = best_Qmat,
                      rate.cat = rate.cat, 
                      root.p = "fitzjohn", 
                      nstarts = 10, n.cores = 5, 
                      get.tip.states = FALSE)
  
  # plot
  message("Plotting Marginal ancestral state reconstruction...")
  pdf(paste0("./results/AncRec/Marginal.ASR_model.fit_", best_discrete_model, "_", trait, ".pdf"), width = width, height = height)
  plotTree(tree,show.tip.label = TRUE, type = "fan", part = 0.88,
           fsize = 0.5, offset = 10, lwd = 3)
  xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(tree)]
  yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(tree)]
  # add ages lines
  obj<-axis(1,pos=-2,at=seq(10,50,by=10),cex.axis=0.5,labels=FALSE)
  text(obj,rep(-5,length(obj)), rev(obj),cex=0.9)
  text(mean(obj),-8,"time (ma)",cex=0.9)
  for(i in 1:(length(obj)-1)){
    a1<-atan(-2/obj[i])
    a2<-0.88*2*pi
    plotrix::draw.arc(0,0,radius=obj[i],a1,a2,lwd=1.5,
                      col=make.transparent("grey40",0.25))
  }
  ## add traits
  points(xx*1.025, yy*1.025, pch = 21, bg = cols[trait.vec[tree$tip.label]])
  # add legend
  colLegend = unique(data.frame(id = as.numeric(fit.marginal$data.legend[,2]), 
                                trait = as.character(trait.vec[fit.marginal$data.legend[,1]])))
  colLegend <- sort(setNames(colLegend[,1], colLegend[,2]))
  legend(x="topleft", legend=names(colLegend),
         pt.cex=2, pch=16, col=cols[names(colLegend)], title = trait_title, bty = "n")
  # add states at nodes
  nodelabels(pie = to.matrix(names(colLegend)[fit.marginal$phy$node.label],
                             names(colLegend)), piecol = cols[names(colLegend)], cex = 0.2)
  dev.off()

  ##############################
  # Stochastic character mapping
  ##############################
  if (mcmc) {
    message("Generating 1000 stochastic character maps in which the transition rate is sampled from its posterior distribution through MCMC...")
    # generate 1000 stochastic character maps in which the transition rate is sampled 
    # from its posterior distribution
    trait_SCM_trees <- phytools::make.simmap(tree, trait.vec,
                                   model = best_Qmat,
                                   nsim = nsim, pi="fitzjohn",
                                   Q = "mcmc", vQ = 0.01,
                                   prior = list(use.empirical = TRUE), samplefreq = 10)
    
    prefix = "MCMC_"
  } else {
    message("Generating 1000 stochastic character maps using the specified transition matrix from the ", best_discrete_model, " model...")
    ## Only if best model is not a Hidden Rate model
    if (best_discrete_model %in% Mk_models){
      # simulate nsim stochastic models with the best discrete evolutionary model
      trait_SCM_trees <- phytools::make.simmap(tree, trait.vec,
                                     model = best_Qmat,
                                     nsim = nsim)
    } else { 
      trait_SCM_trees = NULL
    }
    prefix = ""
  }
  
  # write stochastic maps
  message("Writting stochastic maps...")
  lapply(trait_SCM_trees, write.simmap, 
         file= paste0("./results/AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_n", nsim, ".map"), 
         append = TRUE) 

  # plot posterior density from stochastic mapping
  message("Plotting posterior density from stochastic mapping...")
  pdf(paste0("./results/AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_posterior.density.distribution.pdf"), width = width, height = height)
  par(mar = c(5.1, 4.1, 2.1, 2.1))
  plot(d <- density(sapply(trait_SCM_trees, function(x) x$Q[1,2]), bw = 0.005),
       bty = "n", main = paste("msmc SCM on", trait), xlab = "q", xlim = c(0, 0.5), ylab = "Posterior density from MCMC",
       las = 1, cex.axis = 0.8)
  polygon(d, col = make.transparent("blue", .25))
  abline(v = fit.marginal$solution[1,2])
  text(x = fit.marginal$solution[1,2], y = max(d$y), "MLE (q)", pos = 4)
  dev.off()

  # plot 100 stochastic maps from the 1000
  message("Plotting 100 out of ", nsim, " stochastic maps...")
  pdf(paste0("./results/AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_n100maps.pdf"), width = width, height = height)
  par(mfrow = c(10, 10))
  null <- sapply(trait_SCM_trees[sample(1:100, 100)], plot, colors = cols, lwd = 1, ftype = "off")
  dev.off()
  par(mfrow = c(1,1))
  
  # plot transitions densities
  message("Plotting transitions densities...")
  pdf(paste0("./results/AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_transition_densities.pdf"), width = width, height = height)
  dd<-density(trait_SCM_trees)
  plot(dd)
  dev.off()

  # compute posterior probabilities at nodes
  SCM_pbn <- summary(trait_SCM_trees)
  
  # plot fan
  message("Plotting posterior probabilities at nodes...")
  pdf(paste0("./results/AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_posterior_probabilities_at_nodes.pdf"), width = width, height = height)
  plot(SCM_pbn, colors = cols[colnames(SCM_pbn$ace)], type = "fan", fsize = 0.5, part = 0.88, tip.cols = "white", lwd = 2, cex = c(0.3,0.1), offset = 6)
  xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(tree)]
  yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(tree)]
  # add ages lines
  obj<-axis(1,pos=-2,at=seq(10,50,by=10),cex.axis=0.5,labels=FALSE)
  text(obj,rep(-5,length(obj)), rev(obj),cex=0.9)
  text(mean(obj),-8,"time (ma)",cex=0.9)
  for(i in 1:(length(obj)-1)){
    a1<-atan(-2/obj[i])
    a2<-0.88*2*pi
    plotrix::draw.arc(0,0,radius=obj[i],a1,a2,lwd=1.5,
                      col=make.transparent("grey40",0.25))
  }
  plot(SCM_pbn, colors = cols[colnames(SCM_pbn$ace)], type = "fan", fsize = 0.5, part = 0.88, lwd = 2, cex = c(0.3,0.1), offset = 6, add = T)
  # add legend
  legend(x="topleft", legend=colnames(SCM_pbn$ace),
         pt.cex=2, pch=16, col=cols[colnames(SCM_pbn$ace)], title = trait_title, cex = 1.5,  bty = "n")
  dev.off()

  # Marginal Likelihoods vs Posterior Probabilities
  message("Plotting Marginal Likelihoods vs Posterior Probabilities...")
  pdf(paste0("./results/AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_MargLik_vs_PostProb.pdf"), width = width, height = height)
  par(mar = c(5.1,4.1,2.1,2.1))
  plot(as.vector(fit.marginal$states), as.vector(SCM_pbn$ace[1:tree$Nnode,]), pch = 21,
       cex = 1.2, bg = "grey", xlab = "Marginal scaled likelihoods", 
       ylab = "Posterior probabilities", bty = "n", las = 1, cex.axis = 0.8)
  lines(c(0,1), c(0,1), col = alpha("blue", 0.5), lwd = 2)
  dev.off()
  
  # create densityMap object to visualize the posterior probability of being in 
  # each state across all the edges and nodes of the tree
  message("Plotting density map to visualize the posterior probability of being in each state along the tree...")
  if (nstates == 2){
    # For a binary trait:
    ####################
    pdf(paste0("./results/AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_densityMap.pdf"), width = width, height = height)
    trait.densityMap <- densityMap(trait_SCM_trees, states = levels(trait.vec), plot = FALSE)
    trait.densityMap <- setMap(trait.densityMap, cols)
    plot(trait.densityMap, fsize = c(0.3, 0.7), lwd = c(3, 4))
    dev.off()
    
  } else if (nstates == 3) {
    # For a 3-states trait:
    ####################
    pdf(paste0("./results/AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_densityMap.pdf"), width = width, height = height)
    densityTree(trait_SCM_trees, method="plotSimmap", fix.depth = TRUE, ftype = "off",
                fsize = 36*par()$pin[2]/par()$pin[1]/Ntip(tree)[1],
                lwd=8, alpha = NULL, nodes="intermediate",
                colors=cols,compute.consensus=FALSE, show.axis = FALSE)
    rect(-0.01, par()$usr[3], par()$usr[2], par()$usr[4], col = "white", border = "white")
    h<-max(nodeHeights(tree))
    labs <- seq(0, h, by=5)
    axis(1, pos=-0.025*h, at=labs, labels=labs, tick = TRUE,  cex.axis=2, lwd=2, lend=2)
    text(mean(par()$usr[1:2]), -0.175*h, "Million Years Ago", cex = 3)
    ylims <- sort(quantile(par()$usr[3:4], c(0.95, 0.85)))
    r <- diff(ylims)/diff(par()$usr[3:4]) * 1.5
    a <- diff(par()$usr[1:2]) * r
    xlims <- quantile(par()$usr[1:2], c(0.95, 0.8))
    xlims <- c(xlims[1], xlims[1] + a)
    make.3state.legend(cols, names = levels(trait.vec), fill.center = "grey90", 
                       xlims = xlims, ylims = ylims, size = 0.9, pc = 1, b = 100,
                       border.col = "transparent", offset = 0, cex = 5, add = T)
    dev.off()
  }
  
  ################
  ## SAVE OBJECT
  ################
  save(trait.discrete.fit, fit.joint, fit.marginal, trait_SCM_trees, SCM_pbn, tree_rescaled,
       file = paste0("./rdata/ASR_", trait, ".rda"))
}
