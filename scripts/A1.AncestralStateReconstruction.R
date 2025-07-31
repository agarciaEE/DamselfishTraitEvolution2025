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

setwd("~/Unil/Research/5_Damselfish_evo/")

# load custom functions
source("./Rscripts/custom_functions.R")

# load tree
damsel_tree_subset <- read.tree("./data/damsel_subset.tre")
# load traits
eco_traits_subset <- read.csv("./data/eco_traits_subset_final.csv", row.names = 1)

# get discrete and continuous traits
discrete_traits <- colnames(eco_traits_subset)[sapply(1:ncol(eco_traits_subset), function(i) is.character(eco_traits_subset[,i]))]
discrete_traits <- discrete_traits[discrete_traits != "species"] # remove species as trait
discrete_traits = c("DietEcotype", "Farming", 'habitat1')

continuous_traits <- colnames(eco_traits_subset)[sapply(1:ncol(eco_traits_subset), function(i) is.numeric(eco_traits_subset[,i]))]

# discrete traits as factors
eco_traits_subset$DietEcotype <- factor(eco_traits_subset$DietEcotype, levels = c("B", "I", "P"))
eco_traits_subset$Farming <- factor(eco_traits_subset$Farming, levels = c("No", "Farming", "Mutualism"))
eco_traits_subset$Symbiosis <- factor(as.character(eco_traits_subset$Symbiosis), levels = c("No", "Yes"))
eco_traits_subset$habitat1 <- factor(eco_traits_subset$habitat1, levels = c("rocky-reef", "coral-reef", "seagrass-rubble-sand", "sea anemone", "freshwater"))
eco_traits_subset$waterColumn <- factor(eco_traits_subset$waterColumn, levels = c("benthic/site-attached", "pelagic/non-site-attached"))
eco_traits_subset$type <- factor(eco_traits_subset$type, levels = c("inshore-reef", "both", "offshore-reef"))
eco_traits_subset$DepthRange <- factor(eco_traits_subset$DepthRange, levels =c("shallow", "mid-deep", "deep", "bathyal"))
eco_traits_subset$BodyShape <- factor(eco_traits_subset$BodyShape, levels = unique(eco_traits_subset$BodyShape))
eco_traits_subset$SIZE.3state <- factor(eco_traits_subset$SIZE.3state, levels = c("S", "M", "L"))
eco_traits_subset$SIZE.5state <- factor(eco_traits_subset$SIZE.5state, levels = c("XS", "S", "M", "L", "XL"))
eco_traits_subset$Habitat <- factor(eco_traits_subset$Habitat, levels = unique(eco_traits_subset$Habitat[order(eco_traits_subset$habitat1)]))
eco_traits_subset$ReproMode <- factor(eco_traits_subset$ReproMode, levels = c("dioecism", "protogyny", "protandry"))

################################################################################
####                                                                      ######
#### --------------   Reconstructing ancestral states    ---------------- ######
####                                                                      ######
################################################################################

################################################################################
###  DISCRETE TRAIT  ####
################################################################################
width = 15 # figure width
height = 15 # figure height
nsim = 100
HSM = FALSE # Don't run Hidden States Models for now
mcmc = FALSE # Run SCM on fitted models
cat("trait", "ntips", "nstates", "states", "model", "Lambda", "p value", "\n", file = "./AncRec/traits_parameters.txt")

for (trait in discrete_traits[3]) {
  
  message("Working on ", trait, " trait...")
  trait.vec <- setNames(eco_traits_subset[, trait], rownames(eco_traits_subset))
  # make sure the order of the trait matches the tree
  trait.vec <- trait.vec[damsel_tree_subset$tip.label]
  states <- levels(trait.vec)
  nstates <- length(states)
  trait_title = ifelse(trait == "DietEcotype", "Diet", ifelse(trait == "Farming", "Symbiosis", "Habitat"))
  
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
  barplot(1:length(cols), col = cols)

  # Description of some parameters
  cat("\nNumber of tips:", length(trait.vec))
  cat("\nNumber of states:", nstates)
  cat("\nStates:", states)
  cat("\n")

  message("Plotting phylogeny with trait on tips...")
  pdf(paste0("./Traits/phylo_", trait, ".pdf"), width = width, height = height)
  # plot trait on tips
  plotTree(damsel_tree_subset,show.tip.label = TRUE, type = "fan", part = 0.88,
           fsize = 0.5, offset = 10, lwd = 3)
  xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(damsel_tree_subset)]
  yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(damsel_tree_subset)]
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
  points(xx*1.025, yy*1.025, pch = 21, bg = cols[trait.vec[damsel_tree_subset$tip.label]])
  # add legend
  legend(x="topleft", legend=levels(trait.vec),
         pt.cex=2, pch=16, col=cols, title = trait_title, bty = "n")
  dev.off()
  
  ################
  ## Fit Mk models
  ################
  message("Carrying out Markov k-state (Mk) Models...")
  
  Mk_models <- c("ER", "ARD", "SYM", "DIR1", "DIR2", "DIR3")

  trait_Mk.fit <- lapply(Mk_models, function(i) {
    # make Q matrix
    Q <- make_Q_matrix(i, nstates, states)
    # forbid transitions from freshwater (only transitions to freshwater allowed)
    if (trait == "habitat1") {
      Q[states == "freshwater",] = 0
    }
    # fit Mk
    message("Fitting ", i, " model...")
    fitMk(damsel_tree_subset, trait.vec, model = Q)
  })
  names(trait_Mk.fit) <- Mk_models

  message("Plotting Mk Models transitions...")
  pdf(paste0("./AncRec/", trait, "_Mk.fitted_Q.transitions.pdf"), width = width-5, height = height-5)
  par(mfrow = c(2, 3))
  lapply(names(trait_Mk.fit), function(i) plot(trait_Mk.fit[[i]], main  = i, offset = 0.1, spacer = 0.25, lwd = 2, show.zeros = T,
                                               cex.main = 2, cex.traits = 2, cex.rates = 2))
  dev.off()

  Mk.res.table <- data.frame(model = Mk_models,
                          logLik = sapply(trait_Mk.fit, logLik),
                          k = sapply(trait_Mk.fit, function(m)  attr(logLik(m), "df")),
                          AIC = aic <- sapply(trait_Mk.fit, AIC),
                          delta.AIC = aic-min(aic),
                          weight = unclass(aic.w(aic)))
  
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
    
    # using phytools:
    #trait_HR.fit <- lapply(1:length(HR_models), function(i) fitHRM(damsel_tree_subset, trait.vec, ncat = HRM_matrix[,i], pi="fitzjohn"))
    
    # using corHMM:
    trait_HR.fit <- lapply(1:length(HR_models), function(i) corHMM(damsel_tree_subset, cbind(names(trait.vec), trait.vec), 
                                                                   rate.cat = sum(HRM_matrix[,i])-1, nstarts = 10, n.cores = 5, root.p = "maddfitz"))
    
    names(trait_HR.fit) <- HR_models
    
    message("Plotting Hidden State Models transitions...")
    pdf(paste0("./AncRec/", trait, "_HR.fitted_Q.transitions.pdf"), width = width-5, height = height-5)
    par(mfrow = c(2, 3))
    lapply(names(trait_HR.fit), function(i) plot(as.Qmatrix(trait_HR.fit[[i]]), main  = i, offset = 0.1, spacer = 0.25, lwd = 2, show.zeros = T,
                                                 cex.main = 2, cex.traits = 1, cex.rates = 1))
    dev.off()
    
    HR.res.table <- data.frame(model = HR_models,
                            logLik = sapply(trait_HR.fit, function(i) i$loglik),
                            k = sapply(trait_HR.fit, function(i) i$rate.cat),
                            AIC = aic <- sapply(trait_HR.fit, function(i) i$AIC),
                            delta.AIC = aic-min(aic),
                            weight = unclass(aic.w(aic)))
    
  } else {
    
    trait_HR.fit <- NULL
    HR_models <- NULL
    HR.res.table <- NULL
    
  }
  
  ####
  # build a table with the results of all Mk and HR models
  ####
  trait.discrete.fit <- c(trait_Mk.fit, trait_HR.fit)
  trait.discrete.models <- c(Mk_models, HR_models)
  
  res.table <- rbind(Mk.res.table, HR.res.table)

  # write table to hard drive
  message("Writting discrete models result table on hard disk...")
  write.table(res.table, file = paste0("./AncRec/", trait, "_models.fit.restable.txt"), 
                                       quote = FALSE, sep = "\t", row.names = T, col.names = T)
  
  # Likelihood-Ratio tests
  #sapply(trait.discrete.models[trait.discrete.models %in% best_discrete_model], lrtest, trait.discrete.fit[[best_discrete_model]])
  
  # select best discrete model
  best_discrete_model <- res.table$model[which.min(res.table$AIC)]

  if (best_discrete_model %in% Mk_models) { 
    best_Qmat <- trait.discrete.fit[[best_discrete_model]]$index.matrix
    best_Qmat[is.na(best_Qmat)] = 0
    rate.cat <-  1
  } else {
    best_Qmat <- trait.discrete.fit[[best_discrete_model]]$index.mat
    best_Qmat[is.na(best_Qmat)] = 0
    
    # order Qmat if coming from phytools hrm 
    #   idx <- which(HR_models == best_discrete_model)
    #   rate.cat <- sum(HRM_matrix[,idx])
    #   nc <- ncol(best_Qmat)
    #   hs <- which(HRM_matrix[,idx] == 2)
    #   h_i <- which(1:ns %% 2 == 0)[1]
    #   k <- (1:nc)[!1:nc %in% i]
    #   best_Qmat <- cbind(best_Qmat[,k], best_Qmat[,h_i], matrix(0, nrow = ns, ncol = nstates-length(hs)))
    #   best_Qmat <- rbind(best_Qmat[k,], best_Qmat[h_i,], matrix(0, nrow = nstates-length(hs), ncol = nstates*2))
  }
  
  ###########################################
  # Compute discrete trait phylgenetic signal
  ###########################################
  message("Computing phygenetic signal on ", trait, "...")
  # From https://github.com/mrborges23/delta_statistic DOI:10.1093/bioinformatics/bty800
  # phylogenetic analog of the Shannon entropy for measuring the degree of phylogenetic signal between a categorical trait and a phylogeny
  
  # trait_delta <- delta(trait.vec, damsel_tree_subset, 0.1, 0.0589, 10000, 10, 100)
  # 
  # random_delta <- rep(NA,100)
  # for (i in 1:100){
  #   rtrait <- sample(trait.vec)
  #   random_delta[i] <- delta(rtrait,damsel_tree_subset, 0.1, 0.0589, 10000, 10, 100)
  # }
  # p_value <- (sum(random_delta>deltaA) + 1)/ (length(random_delta) + 1)
  # x <- hist(random_delta)
  # 
  # message("Plotting phygenetic signal test on the trait...")
  # pdf(paste0("./AncRec/phylosig_", trait, ".pdf"), width = width, height = height)
  # par(mfrow = c(1,1))
  # plot(x, xlim = range(c(x$breaks, trait_delta)), ylim = c(0,max(x$counts)*1.1), main = paste0("Significance test of phylogenetic signal \ndiscrete trait: ", trait), xlab = "Delta")
  # abline(v=trait_delta, col="blue")
  # text(trait_delta, max(x$counts)*1.1, paste("p-value =", sprintf("%.2e", p_value)), pos = 4)
  # dev.off()
  if (FALSE){
  # Testing phylogenetic signal on discrete traits following phytools.org (http://blog.phytools.org/2018/02/how-to-fit-tree-transformation-for.html)
  #######################################################################
  # phytools.org function to get logLik of Pagel's lambda on discrete trait
  lk.lambda<-function(lambda,tree,x,...) {
    -logLik(fitMk(phytools:::lambdaTree(tree,lambda),
                  x,...))
  }
  # convert trait vector to matrix
  X <- to.matrix(levels(trait.vec)[trait.vec],
                 levels(trait.vec))
  rownames(X) <- names(trait.vec)
  
  # get optimal lambda estimate using ML
  opt<-optimize(lk.lambda,c(0,phytools:::maxLambda(damsel_tree_subset)),tree=damsel_tree_subset,
                x=X,model=best_Qmat)
  
  lam<-opt$minimum # get estimated lambda
  
  # explore lokLik among a sequence of lambdas to compare with observed
  lambda<-seq(0,1,by=0.01)
  lik<--sapply(lambda,lk.lambda,tree=damsel_tree_subset,x=X,model=best_Qmat)

  # fit best lambda and lambda 0
  fit.lambda<-fitMk(phytools:::lambdaTree(damsel_tree_subset,lam),X,model=best_Qmat)
  fit.h0<-fitMk(phytools:::lambdaTree(damsel_tree_subset,0),X,model=best_Qmat)
  
  # get chi sqr statistic testing against a lambda 0 for phylogenetic signal
  LR<--2*(logLik(fit.h0)-logLik(fit.lambda))
  P.chisq<-as.numeric(pchisq(LR,df=1,lower.tail=FALSE))
  
  # write trait parameters to file
  cat(trait, length(trait.vec), nstates, paste(states, collapse = ","), best_discrete_model, lam, P.chisq, "\n", file = "./AncRec/traits_parameters.txt", append = TRUE)
  
  # plot result
  message("Plotting phygenetic signal test on the trait...")
  pdf(paste0("./AncRec/phylosig_", trait, ".pdf"), width = width, height = height)
  plot(lambda,lik,type="l",xlab=expression(lambda),ylab="log(likelihood)",
       lwd=2,col="darkgrey", main = bquote(paste("Likelihood Surface and ML of ", lambda, " for ", .(trait))))
  abline(v=lam,lty="dashed")
  text(x=lam,y=-1.005*opt$objective,expression(paste("ML(",lambda,")")), pos=3)
  text(x=lam, y= min(lik), labels = bquote(paste(chi^2 , " p-value = " , .(format(P.chisq, scientific = T)), sep="")), pos = 2)
  dev.off()
  }
  #######################################
  # Joint ancestral state reconstructions
  #######################################
  
  fit.joint <- corHMM(damsel_tree_subset, 
                      cbind(species = names(trait.vec), 
                            as.data.frame(trait.vec)), 
                      node.states = "joint", 
                      rate.mat = best_Qmat,
                      rate.cat = rate.cat, 
                      root.p = "maddfitz", 
                      nstarts = 10, n.cores = 5, 
                      get.tip.states = FALSE)

  # plot
  message("Plotting Joint ancestral state reconstruction...")
  pdf(paste0("./AncRec/Joint.ASR_model.fit_", best_discrete_model, "_", trait, ".pdf"), width = width, height = height)
  plotTree(damsel_tree_subset,show.tip.label = TRUE, type = "fan", part = 0.88,
           fsize = 0.5, offset = 10, lwd = 3)
  xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(damsel_tree_subset)]
  yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(damsel_tree_subset)]
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
  points(xx*1.025, yy*1.025, pch = 21, bg = cols[trait.vec[damsel_tree_subset$tip.label]])
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
  
  fit.marginal <- corHMM(damsel_tree_subset, 
                      cbind(species = names(trait.vec), 
                            as.data.frame(trait.vec)), 
                      node.states = "marginal", 
                      rate.mat = best_Qmat,
                      rate.cat = rate.cat, 
                      root.p = "maddfitz", 
                      nstarts = 5, n.cores = 5, 
                      get.tip.states = FALSE)
  
  # plot
  message("Plotting Marginal ancestral state reconstruction...")
  pdf(paste0("./AncRec/Marginal.ASR_model.fit_", best_discrete_model, "_", trait, ".pdf"), width = width, height = height)
  plotTree(damsel_tree_subset,show.tip.label = TRUE, type = "fan", part = 0.88,
           fsize = 0.5, offset = 10, lwd = 3)
  xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(damsel_tree_subset)]
  yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(damsel_tree_subset)]
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
  points(xx*1.025, yy*1.025, pch = 21, bg = cols[trait.vec[damsel_tree_subset$tip.label]])
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
  trait_SCM_trees <- make.simmap(damsel_tree_subset, trait.vec,
                                      model = best_Qmat,
                                      nsim = nsim,
                                      Q = "mcmc", vQ = 0.01,
                                      prior = list(use.empirical = TRUE), samplefreq = 10)
  
  prefix = "MCMC_"
  } else {
    message("Generating 1000 stochastic character maps using the specified transition matrix from the ", best_discrete_model, " model...")
    ## Only if best model is not a Hidden Rate model
    if (best_discrete_model %in% Mk_models){
      # simulate nsim stochastic models with the best discrete evolutionary model
      trait_SCM_trees <- make.simmap(damsel_tree_subset, trait.vec,
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
         file= paste0("./AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_n", nsim, ".map"), 
         append = TRUE) 
  
  # plot posterior density from stochastic mapping
  message("Plotting posterior density from stochastic mapping...")
  pdf(paste0("./AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_posterior.density.distribution.pdf"), width = width, height = height)
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
  pdf(paste0("./AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_n100maps.pdf"), width = width, height = height)
  par(mfrow = c(10, 10))
  null <- sapply(trait_SCM_trees[sample(1:100, 100)], plot, colors = cols, lwd = 1, ftype = "off")
  dev.off()
  par(mfrow = c(1,1))
  
  # plot transitions densities
  message("Plotting transitions densities...")
  pdf(paste0("./AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_transition_densities.pdf"), width = width, height = height)
  dd<-density(trait_SCM_trees)
  plot(dd)
  dev.off()
  states
  # compute posterior probabilities at nodes
  SCM_pbn <- summary(trait_SCM_trees)
  
  # plot fan
  message("Plotting posterior probabilities at nodes...")
  pdf(paste0("./AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_posterior_probabilities_at_nodes.pdf"), width = width, height = height)
  plot(SCM_pbn, colors = cols[colnames(SCM_pbn$ace)], type = "fan", fsize = 0.5, part = 0.88, tip.cols = "white", lwd = 2, cex = c(0.3,0.1), offset = 6)
  xx<-get("last_plot.phylo",envir=.PlotPhyloEnv)$xx[1:Ntip(damsel_tree_subset)]
  yy<-get("last_plot.phylo",envir=.PlotPhyloEnv)$yy[1:Ntip(damsel_tree_subset)]
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
  pdf(paste0("./AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_MargLik_vs_PostProb.pdf"), width = width, height = height)
  par(mar = c(5.1,4.1,2.1,2.1))
  plot(as.vector(fit.marginal$states), as.vector(SCM_pbn$ace[1:damsel_tree_subset$Nnode,]), pch = 21,
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
    pdf(paste0("./AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_densityMap.pdf"), width = width, height = height)
    trait.densityMap <- densityMap(trait_SCM_trees, states = levels(trait.vec), plot = FALSE)
    trait.densityMap <- setMap(trait.densityMap, cols)
    plot(trait.densityMap, fsize = c(0.3, 0.7), lwd = c(3, 4))
    dev.off()
    
  } else if (nstates == 3) {
    # For a 3-states trait:
    ####################
    pdf(paste0("./AncRec/SCM_", prefix, "model.fit_", best_discrete_model,"_", trait, "_densityMap.pdf"), width = width, height = height)
    densityTree(trait_SCM_trees, method="plotSimmap", fix.depth = TRUE, ftype = "off",
                fsize = 36*par()$pin[2]/par()$pin[1]/Ntip(damsel_tree_subset)[1],
                lwd=8, alpha = NULL, nodes="intermediate",
                colors=cols,compute.consensus=FALSE, show.axis = FALSE)
    rect(-0.01, par()$usr[3], par()$usr[2], par()$usr[4], col = "white", border = "white")
    h<-max(nodeHeights(damsel_tree_subset))
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
  save(trait.discrete.fit, fit.joint, fit.marginal, trait_SCM_trees, SCM_pbn,
       file = paste0("./rdata/ASR_", trait, ".rda"))
}

################################################################################
### CONTINUOUS TRAIT ####
################################################################################

## !! check if log-transform of continuous variable is needed!
hist(trait.vec)
log = FALSE
if (log) {
  trait.vec <- log(trait.vec) # appropriate distribution for BM (Normal - like)
}

####
# Check phylogenetic signal
#//////////////////////////

# Using K
#########
K_tl <- phylosig(damsel_tree_subset, MaxTL, test = TRUE, nsim = 1000) # although significant

# plot
par(cex=0.8, mar=c(5.1,4.1,2.1,2.1), mfrow = c(1,1))
plot(K_tl, las=1, cx.axis=0.9)

# simulate null (BM) evolution
nullX <- fastBM(damsel_tree_subset, nsim = 10000)
# test phylosig for each
nullK <- apply(nullX, 2, phylosig, tree = damsel_tree_subset)
# Calculate p-values
Pval_bd <- mean(nullK<=K_bd$K)

# plot
par(mfrow=c(1,1))
hist(c(nullK, K_bd$K), breaks = 30, col="lightgray",
     border="lightgray", main="", xlab="K", las=1,
     cex.axis=0.7, cex.lab=0.9, ylim=c(0,4000))
arrows(x0=K_bd$K, y0=par()$usr[4], y1=0, length=0.12,
       col=make.transparent("blue", 0.5), lwd=2)
text(K_bd$K, 0.96*par()$usr[4],
     paste("observed value of K (P = ",
           round(Pval_bd, 4), ")", sep=""), pos = 4, cex= 0.8)
mtext(paste("Phylogenetic sygnal test of", trait), line=1, adj=0)

# Using lambda
##############

lambda_bd <- phylosig(damsel_tree_subset, BodyDepth, method = "lambda", test = TRUE)
lambda_bd

# Plot likelihood surfaces
par(mfrow=c(3,1), mar=c(5.1,4.1,2.1,2.1), cex=0.8)
plot(lambda_bd, las=1, cex.axis=0.9, bty="n", xlim=c(0,1.1))
mtext(paste("Likelihood surfaces of phylogenetic signal of", trait), line=1, adj=0)

# Check statistical significance
LR_bd <- -2*(lambda_bd$lik(1) - lambda_bd$logL)
LR_bd
Pval_lambda_bd <- pchisq(LR_bd, df=1, lower.tail = FALSE)
Pval_lambda_bd

#########################
# Fit evolutionary models
#########################

# Brownian Motion
#################
fitBM <- geiger::fitContinuous(damsel_tree_subset, trait.vec)

# Early Burst
#############
fitEB_bd <- fitContinuous(damsel_tree_subset, BodyDepth, model = "EB")
fitEB_bd

# Ornstein-Uhlenbeck
####################
fitOU_bd <- fitContinuous(damsel_tree_subset, BodyDepth, model = "OU", 
                          bounds = list(alpha=c(0,1e+7)))
fitOU_bd

#### 
# Create results table
####
aic_bd <- setNames(c(AIC(fitBM_bd), 
                     AIC(fitEB_bd),
                     AIC(fitOU_bd)),
                   c("BM", "EB", "OU"))
aic_bd

aic.w(aic_bd)


################################################################################
####                                                                      ######
#### --------------   Multi-regime evolutionary models  ----------------- ######
####                                                                      ######
################################################################################


# load simmap tree from discrete trait
for (d_trait in discrete_traits){
  load("./AncRec/ASR_", trait, ".rda")
  
  
  # plot
  tips <- getStates(damsel_tree_subset_simmap_diet, "tips")
  tip.cols <- cols[tips]
  plotTree.barplot(damsel_tree_subset_simmap_diet, colpc[,4], 
                   args.plotTree = list(fsize = 0.4), args.barplot = list(col = tip.cols,
                                                                          xlab = "PC1", cex.lab = 0.8))
  legend("topleft", levels(diet), pch = 22, pt.bg = cols, pt.cex = 1.5, cex = 0.9)
  
  
  # Brownian Motion Single-rate (one sigma)
  ###############################################
  fitBM_owie <- OUwie(damsel_tree_subset_simmap_diet, owie.data, model = "BM1", simmap.tree = TRUE,
                      ub=1e+10)
  fitBM_owie
  
  # Brownian Motion Multi-rate (multiple sigma)
  ###############################################
  fitBMS_owie <- OUwie(damsel_tree_subset_simmap_diet, owie.data, model = "BMS", 
                       simmap.tree = TRUE, root.station = FALSE,
                       ub=1e+10)
  fitBMS_owie
  
  # with phytools:
  brownie.lite(damsel_tree_subset_simmap_diet, scores(col.pca)[,1])
  
  # or: also from phytools: 
  fitMV.all <- evolvcv.lite(damsel_tree_subset_simmap_diet, scores.colpca[,1:2])
  fitMV.all
  
  # Ornstein-Uhlenbeck (OU) multi-regime
  ######################################
  fitOUM_owie <- OUwie(damsel_tree_subset_simmap_diet, owie.data, model = "OUM", 
                       simmap.tree = TRUE, root.station = FALSE,
                       ub=1e+10)
  fitOUM_owie
  
  # Ornstein-Uhlenbeck (OU) multi-regime and multiple sigma
  #########################################################
  fitOUM_owie <- OUwie(damsel_tree_subset_simmap_diet, owie.data, model = "OUMV", 
                       simmap.tree = TRUE, root.station = FALSE,
                       ub=1e+10)
  fitOUM_owie
  
  # Ornstein-Uhlenbeck (OU) multi-regime and multiple alpha
  #########################################################
  fitOUM_owie <- OUwie(damsel_tree_subset_simmap_diet, owie.data, model = "OUMA", 
                       simmap.tree = TRUE, root.station = FALSE,
                       ub=1e+10)
  fitOUM_owie
  
  # Ornstein-Uhlenbeck (OU) multi-regime and multiple sigma and alpha
  ###################################################################
  fitOUM_owie <- OUwie(damsel_tree_subset_simmap_diet, owie.data, model = "OUMVA", 
                       simmap.tree = TRUE, root.station = FALSE,
                       ub=1e+10)
  fitOUM_owie
  
  #### 
  # Create results table
  ####
  aic_pc1 <- setNames(c(fitBM_owie$AIC, 
                        fitBMS_owie$AIC,
                        fitOUM_owie$AIC),
                      c("BM1", "BMS", "OUM"))
  aic_pc1
  aic.w(aic_pc1)
  
  
}

################################################################################
# 5.4.1 Testing for temporal shifts in the rate of evolution
################################################################################

# fit single-rate model (no rate shift)
fit1 <- rateshift(damsel_tree_subset, colpc[,1]) # try different starting values or computer seeds to converge

# fit two-rate model (one rate shift)
fit2 <- rateshift(damsel_tree_subset, colpc[,1], nrates = 2) # try different starting values or computer seeds to converge

# fit three-rate model (two rate shift)
fit3 <- rateshift(damsel_tree_subset, colpc[,1], nrates = 3) # try different starting values or computer seeds to converge

# fit EB using geiger::fitContinuous
fitEB <- fitContinuous(damsel_tree_subset, colpc[,1], model = "EB")

fits <- list(fit1, fitEB, fit2, fit3)
res.fits <- data.frame(model = (c("BM", "EB", "two-rate", "three-rate")),
                       logL = sapply(fits, logLik),
                       k = sapply(fits, function(x) attr(logLik(c), "df")),
                       AIC = sapply(fits, AIC),
                       weight = unclass(aic.w(sapply(fits, AIC))))

# plot each fitted model both graphed onto the tree and as a line plot illustrating the change in rate through time
## compute the total heigth of our tree
h <-  max(nodeHeights(damsel_tree_subset))

# split plot into 8 panels
par(mfrow = c(4,2))

# panel a) single-rate model
plot(fit1, mar = c(1.1, 4.1, 2.1, 0.1), ftype = "i", fsize =0.5, col = "gray")
mtext("(a)", adj = 0, line = 0)
# panel b) line graph of the single-rate model
par(mar = c(4.1, 4.1, 2.1, 1.1))
plot(NA, xlim = c(0,h), ylim = c(0,fit1$sig2*2), xlab = "time", ylab= expression(sigma^2), bty = "n")
lines(c(0,h), rep(fit1$sig2,2), lwd = 3, col = "gray")
mtext("(b)", adj = 0, line = 0)

# panel c) EB model
## compute sigma^2 through time under the fitted model
s2 <- fitEB$opt$sigsq*exp(fitEB$opt$a*seq(h/200, h-h/200, length.out = 100))
s2.index <- round((s2-min(s2)) / diff(range(s2)) * 100) + 1
# use make.era.map to pain fitted RB model onto tree
tmp <- make.era.map(damsel_tree_subset, setNames(seq(0, h, length.out = 101), s2.index))
# set colors for graph
cols <- setNames(gray.colors(101, 0.9, 0), 1:101)
# plot tree
plot(tmp, cols, mar=c(1.1, 4.1, 2.1, 0.1), ftype = "i", ylim = c(-0.1*Ntip(damsel_tree_subset), Ntip(damsel_tree_subset)), fsize=0.5)
add.color.bar(leg=0,.5*h, cols = cols, prompt=FALSE, x=0, y=-0.05*Ntip(damsel_tree_subset), lims = round(range(s2), 3), title = expression(sigma^2))
mtext("(c)", adj = 0, line = 0)
# panel d) line graph for the EB model
par(mar = c(4.1, 4.1, 2.1, 1.1))
plot(NA, xlim = c(0,h), ylim = c(0,max(s2)), xlab = "time", ylab= expression(sigma^2), bty = "n")
lines(seq(0,h, length.out = 100), s2, lwd = 3, col = "gray")
mtext("(d)", adj = 0, line = 0)

# panel e) two-rate model
plot(fit2, mar = c(1.1, 4.1, 2.1, 0.1), ftype = "i", fsize =0.5, col = "gray")
mtext("(e)", adj = 0, line = 0)
# panel b) line graph of the two-rate model
par(mar = c(4.1, 4.1, 2.1, 1.1))
plot(NA, xlim = c(0,h), ylim = c(0,12), xlab = "time", ylab= expression(sigma^2), bty = "n")
lines(c(0, fit2$shift, h), c(fit2$sig2, fit2$sig2[2]), type = "s", lwd = 3, col = "gray")
mtext("(f)", adj = 0, line = 0)

# panel g) three-rate model
plot(fit3, mar = c(1.1, 4.1, 2.1, 0.1), ftype = "i", fsize =0.5, col = "gray")
mtext("(e)", adj = 0, line = 0)
# panel h) line graph of the three-rate model
par(mar = c(4.1, 4.1, 2.1, 1.1))
plot(NA, xlim = c(0,h), ylim = c(0,12), xlab = "time", ylab= expression(sigma^2), bty = "n")
lines(c(0, fit3$shift, h), c(fit3$sig2, fit3$sig2[3]), type = "s", lwd = 3, col = "gray")
mtext("(f)", adj = 0, line = 0)

################################################################################
# 5.4.2. Exploring heterogeneity across branches and clades
################################################################################

if(!require('bayou')) devtools::install_github("uyedaj/bayou", auth_token = "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")
library(bayou)

# define priors
par(bty="n")
priorOU <- make.prior(damsel_tree_subset, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"),
                      param=list(dalpha=list(scale=0.1), 
                                 dsig2=list(scale=0.1), 
                                 dk=list(lambda=10, kmax=50), # number of regimes shifts
                                 dsb=list(bmax=1, prob=1), 
                                 dtheta=list(mean=mean(BodyDepth), sd=1.5*sd(BodyDepth))), 
                      plot.prior = TRUE)

startpars <- priorSim(priorOU, damsel_tree_subset, plot = FALSE)$pars[[1]]
startpars

priorOU(startpars)

mcmcOU <- bayou.makeMCMC(damsel_tree_subset, BodyDepth, prior=priorOU, 
                         plot.freq = NULL, file.dir = NULL, ticker.freq = 1000000)

damsel.bd.rjMCMC <- mcmcOU$run(100000)

damsel.bd.rjMCMC <- set.burnin(damsel.bd.rjMCMC, 0.3)

summary(damsel.bd.rjMCMC)

damsel.bd.rjmcmc.result <- summary(damsel.bd.rjMCMC)
damsel.bd.rjmcmc.result$statistics


par(mar= c(1.1,1.1,3.1,0.1))
plotSimmap.mcmc(damsel.bd.rjMCMC, edge.type="regimes", 
                lwd=2, pp.cutoff=0.25,cex=0.6)
mtext("(a)", adj=0, line=1)
plotSimmap.mcmc(damsel.bd.rjMCMC, edge.type="theta", 
                lwd=2, pp.cutoff=0.25,cex=0.6, 
                legend_settings=list(x=0.2*max(nodeHeights(damsel_tree_subset)),
                                     y=0.7*Ntip(damsel_tree_subset)))
mtext("(b)", adj=0, line=1)

##################################################################################################

# plot phenogram to visualize the best evoluationary model to the evolution of the trait
par(mar = c(5.1, 4.1, 2.1, 2.1))
pehnogram(damsel_tree_subset, MaxSizeTL, fsize = 0.6, 
          color = make.transparent("blue", 0.5),
          spread.cost = c(1,0), cex.axis = 0.8, las = 1)

# Estimate Ancestral States using ML and getting 95% confidence intervals
fit.maxTL <- fastAnc(damsel_tree_subset, MaxSizeTL, vars = TRUE, CI = TRUE)
print(fit.maxTL, printlen = 10)

# plot tree
plotTree(damsel_tree_subset, ftype = "i", fsize = 0.5, lwd = 1)
labelnodes(1:damsel_tree_subset$Nnode+Ntip(damsel_tree_subset),
           1:damsel_tree_subset$Nnode+Ntip(damsel_tree_subset),
           interactive = FALSE, cex = 0.5)

# plotting recontructed ancestral states
damsel.contMap <- contMap(damsel_tree_subset, MaxSizeTL, plot = FALSE, lims = c(2.7, 5.8))
# custom colors
damsel.contMap <- setMap(damsel.contMap, c("white", "orange", "black"))

plot(damsel.contMap, sig = 2, fsize = c(0.4, 0.7), lwd = c(2,3),
     leg.txt = "Max Size Total Length (cm)")

# identify tips 
node <- 1
tips <- extract.clade(damsel_tree_subset, node)$tip.label
tips

# prune contmap to retain specific tips
pruned.contMap <- keep.tip.contMap(damsel.contMap, tips)
# plot
plot(pruned.contMap, xlim = c(-2, 90), lwd = c(3,4), fsize = c(0.7, 0.9))
errorbar.contMap(pruned.contMap, lwd = 8)

phenogram(eel.trees[[1]],bsize,lwd=3,colors=
            setNames(c("blue","red"),c("suction","bite")),
          spread.labels=TRUE,spread.cost=c(1,0),fsize=0.6,
          ftype="i")
add.simmap.legend(colors=setNames(c("blue","red"),c("suction","bite")),
                  prompt=FALSE,shape="circle",x=0,y=250)
obj<-summary(eel.trees)
nodelabels(pie=obj$ace,piecol=setNames(c("red","blue"),colnames(obj$ace)),
           cex=0.6)
tiplabels(pie=to.matrix(fmode[eel.tree$tip.label],colnames(obj$ace)),
          piecol=setNames(c("red","blue"),colnames(obj$ace)),cex=0.4)

################################################################################
####                                                                      ######
#### -------- 7. Other Models of discrete character evolution ----------- ######
####                                                                      ######
################################################################################

# Correlated BINARY traits
#//////////////////////////
object <- plotTree.datamatrix(damsel_tree_subset, data.frame(diet = DietEcotype, faming = FarmingEcotype),
                              fsize = 0.6, yexp = 1, header = FALSE, xexp = 1.45, palettes = c("YlOrRd", "PuBuGn"))

leg <- legend(x = "topright", names(object$colors$diet), cex = 0.7, pch = 22, pt.bg = object$colors$diet,
              pt.cex = 1.5, bty = "n", title = "diet")
leg <- legend(x = leg$rect$left * 0.995, y = leg$rect$top-leg$rect$h, names(object$colors$faming), 
              cex = 0.7, pch = 22, pt.bg = object$colors$faming,
              pt.cex = 1.5, bty = "n", title = "faming")

DietEcotype_bin <- DietEcotype
DietEcotype_bin[DietEcotype_bin == "I"] = "P"
DietEcotype_bin <- factor(DietEcotype_bin, levels = c("B", "P"))

FarmingEcotype_bin <- FarmingEcotype
FarmingEcotype_bin[FarmingEcotype_bin == 2] = 1
FarmingEcotype_bin <- factor(FarmingEcotype_bin, levels = 0:1)

diet.fit <- fitPagel(damsel_tree_subset, DietEcotype_bin, FarmingEcotype_bin)
print(diet.fit)
plot(diet.fit, signif = 4, cex.main = 1, cex.sub = 0.8, cex.taits = 0.7, cex.rates = 0.7, lwd = 1)

## !! CAUTION Madison & FitzJohn (2105): 
# unique or singular evolutionary events could lead to significant model fit 
# of the dependent model compared to the independent model
# RECOMMENDATION: ALWAYS use with caution and visualize the characters on the tree to identify potential
# unique evlutionary events.

# Multi-regime Mk model
#//////////////////////
library(lmtest)
damsel_tree_subset_simmap_diet <- make.simmap(damsel_tree_subset, DietEcotype)
damsel_tree_subset_simmap_habitat <- make.simmap(damsel_tree_subset, HabitatEcotype)
damsel_tree_subset_simmap_farming <- make.simmap(damsel_tree_subset, FarmingEcotype)
damsel_tree_subset_simmap_bodyshape <- make.simmap(damsel_tree_subset, BodyShapeEcotype)

fit.multi <- fitmultiMk(damsel_tree_subset_simmap_diet, FarmingEcotype)
print(fit.multi)

fit.single <- fitmultiMk(damsel_tree_subset_simmap_diet, FarmingEcotype, model = "ER")
print(fit.single, digits = 2)

lrtest(fit.single, fit.multi)

H <- numDeriv::hessian(fit.multi$lik, fit.multi$rates)
v <- diag(solve(-H))
se <- sqrt(v)
se
#plot
par(bty = "n")
stripchart(fit.multi$rates~fit.multi$regimes, vertical = TRUE, 
           bty = "n", ylim = c(-0.5, 0.5), pch = 21, cex = 1.2, bg = "gray", ylab = "Estimated rate, q", 
           cex.lab =0.8, ce.axis = 0.7)
abline(h= 0, lty = "dotted")
for (i in 1:length(se)){
  lines(x = rep(i,3), y = c(fit.multi$rates[1] - se[i], fit.multi$rates[i], fit.multi$rates[i] + se[i]))
  points(i, fit.multi$rates[i], pch = 21, cex = 1.2, bg = "gray")
  lines(c(i-0.05, i+0.05), rep(fit.multi$rates[i]-se[i], 2))
  lines(c(i-0.05, i+0.05), rep(fit.multi$rates[i]+se[i], 2))
  
}

# 7.3. Modelling rate variation using hidden-rate models
################################################################################
# SIMULATE BINARY CHARACTER WIHTOUT HIDDEN RATE
set.seed(7)
# create a transition matrix between states under a simple Mk model
Q.mk <- matrix(c(-1,1,1,-1), 2, 2, dimnames = list(0:1, 0:1))
# simulate characteri history under constan rate model
mk.tree <- sim.history(tree <- pbtree(n=100, scale = 2), Q.mk, anc ="0")
# SIMULATE THREE STATES CHARACTER WITH ONE OF THEM HIDDEN 
# create a hidden-rate transision matrix
# -- this matrix has two different values for character 1: 1 and 1* --
Q.hrm <- matrix(c(-1,1,0,1,-1.5, 0.5, 0, 0.1, -0.1), 3, 3, byrow = TRUE, dimnames = list(c(0:1, "1*"), c(0:1, "1*")))
# simulate under the hidden rate model
hrm.tree <- sim.history(tree, Q.hrm, anc ="0", message = FALSE)

# visualize results
par(mfrow = c(1, 3))
# mk model
cols <- setNames(c("lightgray", "black"), 0:1)
plot(mk.tree, colors = cols, ftype = "off", mar = c(1.1, 2.1, 3.1, 0.1))
legend("bottomleft", names(cols), pch = 15, col = cols, pt.cex = 2, bty = "n")
mtext("(a)", line = 0, adj = 0)
# hidden model
cols <- setNames(c("lightgray", "black", "slategray"), c(0:1, "1*"))
plot(hrm.tree, colors = cols, ftype = "off", mar = c(1.1, 2.1, 3.1, 0.1))
legend("bottomleft", names(cols), pch = 15, col = cols, pt.cex = 2, bty = "n")
mtext("(b)", line = 0, adj = 0)
# plot hidden model but with two 1 states merged
cols <- setNames(c("lightgray", "black"), 0:1)
plot(tree <- mergeMappedStates(hrm.tree, c("1", "1*"), "1"), colors = cols, 
     ftype = "off", mar = c(1.1, 2.1, 3.1, 0.1))
legend("bottomleft", c("0", "1/1*"), pch = 15, col = cols, pt.cex = 2, bty = "n")
mtext("(c)", line = 0, adj = 0)

# HRM with one hidden states for each category
diet.hrm1 <- fitHRM(damsel_tree_subset, DietEcotype_bin, ncat = 2, model= "ARD", 
                    umbral = TRUE, pi = "fitzjohn", niter = 5, opt.method = "nlminb")
print(diet.hrm1, digits = 4)
# HRM with one hidden states for only the second category
diet.hrm2 <- fitHRM(damsel_tree_subset, DietEcotype_bin, ncat = c(1, 2), model= "ARD", 
                    umbral = TRUE, pi = "fitzjohn", niter = 5, opt.method = "nlminb")
print(diet.hrm2, digits = 4)
# HRM with one hidden states for only the first category
diet.hrm3 <- fitHRM(damsel_tree_subset, DietEcotype_bin, ncat = c(2, 1), model= "ARD", 
                    umbral = TRUE, pi = "fitzjohn", niter = 5, opt.method = "nlminb")
print(diet.hrm3, digits = 4)
# HRM with fitting standard Mk model
diet.hrm4 <- fitHRM(damsel_tree_subset, DietEcotype_bin, ncat = 1, model= "ARD", 
                    umbral = TRUE, pi = "fitzjohn", niter = 1, opt.method = "nlminb")
print(diet.hrm4, digits = 4)

# plot models
par(mfrow = c(2,2))
plot(diet.hrm1, spacer = 0.25, mar = c(0.1, 1.1, 2.1, 0.1))
mtext("(a)", line = 0, adj = 0)
plot(diet.hrm2, spacer = 0.25, mar = c(0.1, 1.1, 2.1, 0.1))
mtext("(b)", line = 0, adj = 0)
plot(diet.hrm3, spacer = 0.25, mar = c(0.1, 1.1, 2.1, 0.1))
mtext("(c)", line = 0, adj = 0)
plot(diet.hrm4, spacer = 0.25, mar = c(0.1, 1.1, 2.1, 0.1))
mtext("(d)", line = 0, adj = 0)

# store models stats
hrm_res <- data.frame(model = c("4-state HRM", "Pelagic hidden", "Benthic hidden", "Mk model"),
                      logL = sapply(list(diet.hrm1, diet.hrm2, diet.hrm3, diet.hrm4), logLik),
                      k = sapply(list(diet.hrm1, diet.hrm2, diet.hrm3, diet.hrm4), function(x) length(x$rates)),
                      AIC = sapply(list(diet.hrm1, diet.hrm2, diet.hrm3, diet.hrm4), AIC))
hrm_res
library(corHMM)
fit.diet <- corHMM(damsel_tree_subset, data.frame(Genus.species = names(DietEcotype_bin), 
                                                  diet = as.integer(DietEcotype_bin)), 
                   rate.cat = 3, nstarts = 10, root.p = "maddfitz")

plotMKmodel(fit.diet, display = "square", text.scale = 0.5, vertex.scale = 0.6, arrow.scale = 0.5)
plot(as.Qmatrix(fit.diet), show.zeros = FALSE, lwd = 1, cex.traits = 0.7)

## plot hidden states onto tree
# create matrix containing the tip and internal node states
states <- rbind(fit.diet$tip.states[damsel_tree_subset$tip.label,],
                fit.diet$states)

rownames(states) <- 1:max(damsel_tree_subset$edge)
# normalize each row to sum to 1.0
states <- t(apply(states, 1, function(x) x/(sum(x))))
# set colors
reds <- c("#ec9488", "#eb5a46", "#933b27")
blues <- c("#8bbdd0", "#0079bf", "#094c72")
COLS <- c(reds[1], blues[1], reds[2], blues[2], reds[3], blues[3])
dev.off()
for (i in 1:ncol(states)){
  tree <- damsel_tree_subset
  edge.col <- rep(NA, nrow(tree$edge))
  for (j in 1:nrow(tree$edge)){
    edge.col[j] <- make.transparent(COLS[i], mean(states[tree$edge[j,],i]))
    tree <- paintBranches(tree, tree$edge[j,2],
                          as.character(j))
  }
  cols <- setNames(edge.col, 1:nrow(tree$edge))
  plot(tree, type = "fan", colors = cols, ftype = "off", lwd = 1, add = (i!=1))
}
# add node labels
par(fg = "transparent")
nodelabels(pie = fit.diet$states, piecol = COLS, cex = 0.2)
par(fg = "black")
legend("topleft", rownames(fit.diet$solution)[c(1,3,5,2,4,6)], 
       pch =15, col = COLS[c(1,3,5,2,4,6)],
       pt.cex = 2, bty = "n")

# 7.4.1. Transient models
################################################################################
# polymotphism is an inherently less stable condition than monomorphism
# Assumes polymorphic state is acquired at one (constant) rate and lost at another (presumably faster) constant rate.
dev.off()
plotTree(damsel_tree_subset, ftype = "i", fsize = 0.8, offset = 0.5)
diet_poly <- diet
diet_poly <- as.character(diet_poly)
diet_poly[diet_poly == "I"] = "B+P"
diet_poly <- factor(setNames(diet_poly, names(diet)), levels = c("B", "B+P", "P"))
xx <- strsplit(as.character(diet_poly), split = "+", fixed = TRUE)
pp <- matrix(0,length(diet_poly), 2, dimnames = list(names(diet_poly), c("B", "P")))
for (i in 1:nrow(pp)) pp[i,xx[[i]]] <- 1/length(xx[[i]])
tiplabels(pie = pp, piecol = c("white", "black"), cex = 0.2) 
legend("topleft", c("B", "P"), pch = 21, pt.bg = c("white", "black"), pt.cex = 2, bty = "n")
diet_poly.ER <- fitpolyMk(damsel_tree_subset, diet_poly, model = "ER", quiet = TRUE)
diet_poly.SYM <- fitpolyMk(damsel_tree_subset, diet_poly, model = "SYM", quiet = TRUE)
diet_poly.ARD <- fitpolyMk(damsel_tree_subset, diet_poly, model = "ARD", quiet = TRUE)
diet_poly.transient <- fitpolyMk(damsel_tree_subset, diet_poly, model = "transient", quiet = TRUE)

diet_poly_res <- data.frame(model = c("ER", "SYM", "ARD", "transient"),
                            logLik = c(logLik(diet_poly.ER), logLik(diet_poly.SYM), logLik(diet_poly.ARD), logLik(diet_poly.transient)),
                            k = c(attr(logLik(diet_poly.ER), "df"), attr(logLik(diet_poly.SYM), "df"), 
                                  attr(logLik(diet_poly.ARD), "df"), attr(logLik(diet_poly.transient), "df")),
                            AIC = aic <- c(AIC(diet_poly.ER), AIC(diet_poly.SYM), AIC(diet_poly.ARD), AIC(diet_poly.transient)),
                            weight = unclass(aic.w(aic)))

plot(diet_poly.ARD, signif = 2, mar = rep(1.1,4), cex.trait = 0.6, cex.rates = 0.6)

# 7.5. Threshold models
################################################################################
# Discrete character state is determined by the value of an unobserved continuous trait called 
# 'liability'. Whenever 'liability' crosses the threshold, the discrete character changes 
# state.
# set number of generations for the MCMC
ngen <- 4e6
# run the MCMC in threshBayes
library(phytools)
mcmc.damsel <- threshBayes(damsel_tree_subset, eco_traits_subset[,2:3], type = c("disc", "disc"),
                           ngen = ngen, plot = FALSE, control = list(print.interval=ngen/10))
mcmc.damsel
plot(mcmc.damsel)
par(mar = c(5.1, 4.1, 2.1, 2.1))
plot(density(mcmc.damsel), cex.lab = 0.8, cex.axis = 0.7)

#load coda package
library(coda)
## extract post burn-in sample for r
r.mcmc <- tail(mcmc.damsel$par$r,
               0.8*nrow(mcmc.damsel$par))
class(r.mcmc) <- "mcmc"
hpd.r <- HPDinterval(r.mcmc)
hpd.r
print(mcmc.damsel)

plot(density(mcmc.damsel), cex.lab = 0.8, cex.axis = 0.7)
h <- 0-par()$usr[3]
lines(x = hpd.r, y = rep(-h/2,2))
lines(x = rep(hpd.r[1],2), y = c(-0.3, -0.7)*h)
lines(x = rep(hpd.r[2],2), y = c(-0.3, -0.7)*h)

dev.off()