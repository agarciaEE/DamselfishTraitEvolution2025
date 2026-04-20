lib <- c("ape","picante","sp","geiger","wesanderson",
         "hypervolume","car","bayou","nlme","l1ou",
         "parallel","genlasso","mvMORPH",
         "lme4","mgcv","phytools","corHMM",
         "lmtest","ColorAR","castor","tibble",
         "ggplot2","ggpubr","gridExtra","knitr","kableExtra",
         "MuMIn","dplyr", "tidyr", "stringr")   
sapply(lib, library, character.only = TRUE)

#==============#
#   Settings
#==============#
# set working directory
setwd("~/Unil/Research/5_Damselfish_evo/DamselTraitEvol2025/")

res_folder <- "results"
dir.create(res_folder)
dir.create("figures/univariate/", recursive = TRUE)
dir.create("results/univariate/", recursive = TRUE)

nsim <- 1000

# load input data and custom functions
source("./scripts/A1.Input_data.R")

#===========================================#
#   PGLS with correlation model selection   #
#===========================================#
fit_one_pgls <- function(formula, dat, tree) {
  fit_with_starts <- function(struct = c("Pagel","OU"),
                              starts = c(0.1, 0.5, 0.9)) {
    struct <- match.arg(struct)
    fits <- list(); lls <- rep(-Inf, length(starts))
    for (i in seq_along(starts)) {
      st <- starts[i]
      mod <- try(
        gls(formula, data = dat,
            correlation = switch(struct,
                                 "Pagel" = corPagel(value = st, phy = tree, fixed = FALSE, form = ~ species),
                                 "OU"    = corMartins(value = st, phy = tree, fixed = FALSE, form = ~ species)
            ),
            method = "REML"),
        silent = TRUE
      )
      if (inherits(mod, "gls")) {
        fits[[i]] <- mod
        lls[i]    <- as.numeric(logLik(mod))
      } else {
        fits[[i]] <- NULL
      }
    }
    if (all(!sapply(fits, inherits, "gls"))) return(NULL)
    fits[[which.max(lls)]]
  }
  
  m_bm <- try(
    gls(formula, data = dat,
        correlation = corBrownian(phy = tree, form = ~ species),
        method = "REML"),
    silent = TRUE
  )
  if (!inherits(m_bm, "gls")) m_bm <- NULL
  
  m_pag <- fit_with_starts("Pagel", starts = c(0.1, 0.5, 0.9))
  m_ou  <- fit_with_starts("OU",    starts = c(0.1, 0.5, 0.9))
  
  fits <- list(BM = m_bm, Pagel = m_pag, OU = m_ou)
  fits <- fits[sapply(fits, inherits, "gls")]
  if (length(fits) == 0L) stop("All correlation structures failed for this model.")
  
  # AICc model selection
  if (length(fits) == 2) {
    aictab <- MuMIn::AICc(fits[[1]], fits[[2]])
  } else if (length(fits) == 3) {
    aictab <- MuMIn::AICc(fits[[1]], fits[[2]], fits[[3]])
  } else {
    df1 <- tryCatch(fits[[1]]$dims$p + 1, error = function(e) NA_real_)
    aictab <- data.frame(df = df1, AICc = MuMIn::AICc(fits[[1]]))
    rownames(aictab) <- names(fits)[1]
  }
  aictab <- aictab[, c("df","AICc")]
  aictab$Model <- names(fits)
  rownames(aictab) <- NULL
  aictab <- aictab[order(aictab$AICc), , drop = FALSE]
  aictab$delta  <- aictab$AICc - min(aictab$AICc)
  aictab$weight <- exp(-0.5 * aictab$delta)
  aictab$weight <- aictab$weight / sum(aictab$weight)
  
  best_name <- aictab$Model[1]
  best_mod  <- fits[[best_name]]
  
  # extract lambda/alpha if applicable
  lambda <- NA_real_; alpha <- NA_real_
  if ("Pagel" %in% names(fits)) {
    lambda <- try(as.numeric(coef(fits[["Pagel"]]$modelStruct$corStruct, unconstrained = TRUE)),
                  silent = TRUE)
    if (inherits(lambda, "try-error")) lambda <- NA_real_
  } 
  if ("OU" %in% names(fits)) {
    alpha  <- try(as.numeric(coef(fits[["OU"]]$modelStruct$corStruct, unconstrained = TRUE)),
                  silent = TRUE)
    if (inherits(alpha, "try-error")) alpha <- NA_real_
  }
  
  aictab$param[aictab$Model == "Pagel"] = lambda
  aictab$param[aictab$Model == "OU"] = alpha
  
  # coefficients table
  s  <- summary(best_mod)
  tt <- as.data.frame(s$tTable)
  tt$term <- rownames(tt)
  colnames(tt)[1:4] <- c("estimate","SE","t","p")
  
  # 95% CIs for fixed effects
  ci <- try(intervals(best_mod, which = "coef"), silent = TRUE)
  if (!inherits(ci, "try-error") && !is.null(ci$coef)) {
    ci_df <- as.data.frame(ci$coef)
    ci_df$term <- rownames(ci_df)
    names(ci_df)[1:3] <- c("lower","estimate_ci","upper")
    tt <- dplyr::left_join(tt, ci_df[, c("term","lower","upper")], by = "term")
  } else {
    tt$lower <- NA_real_
    tt$upper <- NA_real_
  }
  
  list(
    best      = best_mod,
    best_name = best_name,
    aictab    = aictab[, c("Model", "param", "df","AICc","delta","weight")],
    coeftab   = tt
  )
}

run_pgls_models <- function(PC_matrix, traits_df, tree, add_FDR = TRUE) {
  
  candidate_rhs <- list(
    `NULL`            = "1",
    DietEcotype       = "DietEcotype",
    Habitat           = "Habitat",
    Symbiosis         = "Symbiosis",
    Habitat_Diet      = "Habitat + DietEcotype",
    Symbiosis_Diet    = "Symbiosis + DietEcotype",
    Habitat_Symbiosis = "Habitat + Symbiosis"
  )
  
  out_best_models <- vector("list", ncol(PC_matrix))
  out_best_coefs  <- list()
  out_pred_aic    <- list()
  out_corr_aic    <- list()
  
  for (i in seq_len(ncol(PC_matrix))) {
    resp_name <- colnames(PC_matrix)[i]
    cat("Fitting candidate PGLS for", resp_name, "...\n")
    
    response <- scale(PC_matrix[, i])
    model_data <- data.frame(response = as.numeric(response), traits_df, check.names = FALSE)
    rownames(model_data) <- rownames(traits_df)
    model_data$species   <- rownames(model_data)
    
    if ("DietEcotype" %in% names(model_data)) {
      model_data$DietEcotype <- factor(model_data$DietEcotype,
                                       levels = c("Intermediate","Benthic","Pelagic"))
    }
    if ("Habitat" %in% names(model_data)) {
      model_data$Habitat <- factor(model_data$Habitat,
                                    levels = c("non-reef","coral-reef","rocky-reef","sea anemone","freshwater"))
    }
    if ("Symbiosis" %in% names(model_data)) {
      model_data$Symbiosis <- factor(model_data$Symbiosis,
                                   levels = c("Comensalistic","Free-living","Mutualistic"))
    }
    
    cand_fits <- list()
    cand_aic_best <- list()
    
    for (nm in names(candidate_rhs)) {
      rhs <- candidate_rhs[[nm]]
      form <- as.formula(paste("response ~", rhs))
      fit  <- try(fit_one_pgls(form, model_data, tree), silent = TRUE)
      if (inherits(fit, "try-error")) next
      
      fit$aictab$predictor <- nm
      # for predictor model nm, take the min AICc across corStructs
      aic_best_row <- fit$aictab[which.min(fit$aictab$AICc), , drop = FALSE]
      
      cand_aic_best[[nm]] <- aic_best_row
      cand_fits[[nm]] <- fit
    }
    
    if (length(cand_fits) == 0L) {
      warning("No candidate model succeeded for ", resp_name, ". Skipping.")
      next
    }
    
    # select best predictor model by AICc
    aic_tbl <- do.call(rbind, cand_aic_best)
    rownames(aic_tbl) <- NULL
    aic_tbl$delta_pred  <- aic_tbl$AICc - min(aic_tbl$AICc)
    aic_tbl$weight_pred <- exp(-0.5 * aic_tbl$delta_pred)
    aic_tbl$weight_pred <- aic_tbl$weight_pred / sum(aic_tbl$weight_pred)
    
    best_pred <- aic_tbl$predictor[which.min(aic_tbl$AICc)]
    best_fit  <- cand_fits[[best_pred]]
    
    # store best model object
    out_best_models[[i]] <- best_fit$best
    
    # coefficients
    coef_i <- best_fit$coeftab
    coef_i$PC <- resp_name
    coef_i$predictor  <- best_pred
    coef_i$best_correlation <- best_fit$best_name
    coef_i$param <- best_fit$aictab$param[1]
    out_best_coefs[[i]] <- coef_i
    
    # keep correlation-structure AICc table for the best predictor
    corr_aic_i <- best_fit$aictab
    corr_aic_i$PC <- resp_name
    corr_aic_i$predictor <- best_pred
    out_corr_aic[[i]] <- corr_aic_i
    
    # keep across-predictor comparison 
    aic_tbl$PC <- resp_name
    out_pred_aic[[i]] <- aic_tbl
  }
  
  coef_df <- dplyr::bind_rows(out_best_coefs)
  pred_aic_df  <- dplyr::bind_rows(out_pred_aic)
  corr_aic_df <- dplyr::bind_rows(out_corr_aic)
  
  if (add_FDR && nrow(coef_df) > 0) {
    coef_df <- coef_df %>%
      dplyr::group_by(PC) %>%
      dplyr::mutate(p_FDR = p.adjust(replace(p, term == "(Intercept)", NA), method = "fdr")) %>%
      dplyr::ungroup()
  }
  
  list(
    models     = out_best_models,
    coefs      = coef_df,
    aic_pred   = pred_aic_df,  # across predictor comparisons
    aic_corr   = corr_aic_df   # corStruct AICc comparisons
  )
}

col_pgls <- run_pgls_models(colorPCs, eco_traits_subset, tree)
morph_pgls <- run_pgls_models(morphoPCs[,1:8], eco_traits_subset, tree)

# Save outputs
saveRDS(col_pgls, "./rdata/uni_col_pgls.rds")
saveRDS(morph_pgls, "./rdata/uni_morpho_pgls.rds")

# --- small formatter helpers ---
fmt_num <- function(x, digits = 3) ifelse(is.na(x), "---", formatC(x, format = "f", digits = digits))
fmt_p   <- function(p) ifelse(is.na(p), "---",
                              ifelse(p < 0.001, "\\textless{}0.001*", ifelse(p<0.05, 
                                                                  paste0(formatC(p, format = "f", digits = 3), "*"),
                                                                  formatC(p, format = "f", digits = 3))))
# ===========================
# COEFFICIENTS (best model per PC)
# ===========================
make_coef_tables <- function(tab, outdir = "tex_pgls", filename = "pgls_coefs.tex", caption = "") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  coef_df <- tab$coefs %>%
    mutate(
      lambda = ifelse(best_correlation == "Pagel", param, NA_real_),
      alpha  = ifelse(best_correlation == "OU",    param, NA_real_),
      term   = case_when(
        term == "(Intercept)" ~ "(Intercept)",
        TRUE ~ gsub("Symbiosis|DietEcotype|Habitat", "", term)
      )
    ) %>%
    select(PC, predictor, best_correlation, lambda, alpha,
           term, estimate, SE, t, p, lower, upper, p_FDR)

  # ombined table across PCs 
  combined <- coef_df %>%
    mutate(
      Model = paste0(sub("_", "+", predictor), " | ", best_correlation,
                     ifelse(!is.na(lambda), paste0(" ($\\lambda$=", fmt_num(lambda,3), ")"), ""),
                     ifelse(!is.na(alpha),  paste0(" ($\\alpha$=", fmt_num(alpha,3), ")"), ""))
    ) %>%
    select(PC, Model, term, estimate, SE, t, p, p_FDR, lower, upper) %>%
    mutate(
      Estimate = fmt_num(estimate, 3),
      SE = fmt_num(SE, 3),
      t  = fmt_num(t, 2),
      p  = fmt_p(p),
      FDR = fmt_p(p_FDR),
      CI = ifelse(is.na(lower) | is.na(upper), "---",
                        paste0("[", fmt_num(lower, 3), ", ", fmt_num(upper, 3), "]"))
    ) %>%
    select(PC, Model, Term = term, Estimate, SE, CI, t, p,FDR)
  colnames(combined) <- c("PC","Model","Term","Estimate","SE","$95\\%$ CI","t","p","p(FDR)")
  
  combined %>%
    kable(format = "latex", booktabs = TRUE, escape = FALSE,
          caption = caption) %>%
    kable_styling(latex_options = c("hold_position")) %>%
    save_kable(file.path(outdir, filename))
}

# ===========================
# REDICTOR COMPARISONS 
# ===========================
make_combined_pred_table <- function(pred,
                                     outdir = "tex_pgls",
                                     out_tex = "pred_model_comparison_combined.tex",
                                     caption = "Univariate PGLS predictor model comparison across colour PCs (best correlation structure per predictor).",
                                     label = "tab:pred_model_comparison_combined") {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  stopifnot(all(c("PC","predictor","Model","param","df","AICc","delta_pred","weight") %in% names(pred)))
  
  tab <- pred %>%
    mutate(
      # nice names
      Predictor  = sub("_", "+", predictor),
      CorStruct  = dplyr::case_when(
        Model == "Pagel" ~ "Pagel ($\\lambda$)",
        Model == "OU"    ~ "OU ($\\alpha$)",
        Model == "BM"    ~ "BM",
        TRUE ~ Model
      ),
      Param      = dplyr::case_when(
        Model == "BM" ~ "---",
        is.na(param)      ~ "---",
        TRUE              ~ sprintf("%.3f", param)
      ),
      # numeric formatting
      AICc       = sprintf("%.3f", AICc),
      DeltaAICc  = sprintf("%.3f", delta_pred),
      Weight     = sprintf("%.3f", weight_pred)
    ) %>%
    arrange(PC, as.numeric(DeltaAICc)) %>%
    group_by(PC) %>%
    mutate(is_best = as.numeric(DeltaAICc) == min(as.numeric(DeltaAICc))) %>%
    ungroup() %>%
    mutate(
      Predictor  = ifelse(is_best, paste0("\\textbf{", Predictor, "}"), Predictor),
      CorStruct  = ifelse(is_best, paste0("\\textbf{", CorStruct, "}"), CorStruct),
      Param      = ifelse(is_best, paste0("\\textbf{", Param, "}"), Param),
      AICc       = ifelse(is_best, paste0("\\textbf{", AICc, "}"), AICc),
      DeltaAICc  = ifelse(is_best, paste0("\\textbf{", DeltaAICc, "}"), DeltaAICc),
      Weight     = ifelse(is_best, paste0("\\textbf{", Weight, "}"), Weight)
    ) %>%
    select(
      PC, Predictor, CorStruct, Param, df,
      AICc, DeltaAICc, Weight
    )
  
  # build latex
  latex_tab <-
    kable(tab,
          format = "latex",
          booktabs = TRUE,
          escape = FALSE,
          caption = caption,
          label = label,
          col.names = c("PC","Predictor","Correlation","Param","df","AICc","$\\Delta$AICc","Weight"),
          align = c("l","l","l","c","c","r","r","r")) |>
    kable_styling(latex_options = c("hold_position")) |>
    collapse_rows(columns = 1, latex_hline = "major")
  
  save_kable(latex_tab, file.path(outdir, out_tex))
  invisible(latex_tab)
}

# ===========================
# CORRELATION-STRUCTURE COMPARISONS 
# ===========================
make_corr_cmp_tables <- function(pred, outdir = "tex_pgls",
                                 out_tex = "corr_model_comparison_combined.tex",
                                 caption = "Univariate PGLS correlation structure model comparison across colour PCs.",
                                 label = "tab:corr_model_comparison_combined") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  stopifnot(all(c("Model","param","df","AICc","delta","weight") %in% names(pred)))
  
  tab <- pred %>%
    mutate(
      # nice names
      Predictor  = sub("_", "+", predictor),
      CorStruct  = dplyr::case_when(
        Model == "Pagel" ~ "Pagel ($\\lambda$)",
        Model == "OU"    ~ "OU ($\\alpha$)",
        Model == "BM"    ~ "BM",
        TRUE ~ Model
      ),
      Param      = dplyr::case_when(
        Model == "BM" ~ "---",
        is.na(param)      ~ "---",
        TRUE              ~ sprintf("%.3f", param)
      ),
      # numeric formatting
      AICc       = sprintf("%.3f", AICc),
      DeltaAICc  = sprintf("%.3f", delta),
      Weight     = sprintf("%.3f", weight)
    ) %>%
    arrange(PC, as.numeric(DeltaAICc)) %>%
    group_by(PC) %>%
    mutate(is_best = as.numeric(DeltaAICc) == min(as.numeric(DeltaAICc))) %>%
    ungroup() %>%
    mutate(
      Predictor  = ifelse(is_best, paste0("\\textbf{", Predictor, "}"), Predictor),
      CorStruct  = ifelse(is_best, paste0("\\textbf{", CorStruct, "}"), CorStruct),
      Param      = ifelse(is_best, paste0("\\textbf{", Param, "}"), Param),
      AICc       = ifelse(is_best, paste0("\\textbf{", AICc, "}"), AICc),
      DeltaAICc  = ifelse(is_best, paste0("\\textbf{", DeltaAICc, "}"), DeltaAICc),
      Weight     = ifelse(is_best, paste0("\\textbf{", Weight, "}"), Weight)
    ) %>%
    select(
      PC, Predictor, CorStruct, Param, df,
      AICc, DeltaAICc, Weight
    )
  
  # build latex
  latex_tab <-
    kable(tab,
          format = "latex",
          booktabs = TRUE,
          escape = FALSE,
          caption = caption,
          label = label,
          col.names = c("PC","Predictor","Correlation","Param","df","AICc","$\\Delta$AICc","Weight"),
          align = c("l","l","l","c","c","r","r","r")) |>
    kable_styling(latex_options = c("hold_position")) |>
    collapse_rows(columns = 1, latex_hline = "major")
  
  save_kable(latex_tab, file.path(outdir, out_tex))
  invisible(latex_tab)
}

# ===========================
# RUN: build all LaTeX tables for colour PGLS
# ===========================
make_coef_tables(col_pgls, outdir = "results/univariate", filename = "col_pgls_coefs.tex", caption = "")
make_coef_tables(morph_pgls, outdir = "results/univariate", filename = "morpho_pgls_coefs.tex", caption = "")

make_combined_pred_table(col_pgls$aic_pred,
                         outdir = "results/univariate",
                         out_tex = "color_pred_model_comparison_combined.tex",
                         caption = "Univariate PGLS predictor model comparison across colour PCs (best correlation structure per predictor). Param is $\\lambda$ for Pagel and $\\alpha$ for OU; BM has no parameter. Bold indicates the best (lowest AICc) model within each PC.",
                         label = "tab:color_pred_model_comparison_combined")
make_combined_pred_table(morph_pgls$aic_pred,
                         outdir = "results/univariate",
                         out_tex = "morpho_pred_model_comparison_combined.tex",
                         caption = "Univariate PGLS predictor model comparison across morphology PCs (best correlation structure per predictor). Param is $\\lambda$ for Pagel and $\\alpha$ for OU; BM has no parameter. Bold indicates the best (lowest AICc) model within each PC.",
                         label = "tab:color_pred_model_comparison_combined")

make_corr_cmp_tables(col_pgls$aic_corr,
                         outdir = "results/univariate",
                         out_tex = "color_corr_model_comparison_combined.tex",
                         caption = "Univariate PGLS correlation structure model comparison across colour PCs. Param is $\\lambda$ for Pagel and $\\alpha$ for OU; BM has no parameter. Bold indicates the best (lowest AICc) model within each PC.",
                         label = "tab:color_corr_model_comparison_combined")
make_corr_cmp_tables(morph_pgls$aic_corr,
                         outdir = "results/univariate",
                         out_tex = "morpho_corr_model_comparison_combined.tex",
                         caption = "Univariate PGLS correlation structure model comparison across morphology PCs. Param is $\\lambda$ for Pagel and $\\alpha$ for OU; BM has no parameter. Bold indicates the best (lowest AICc) model within each PC.",
                         label = "tab:color_corr_model_comparison_combined")

na.omit(col_pgls$coefs[col_pgls$coefs$p < 0.05,])
na.omit(col_pgls$coefs[col_pgls$coefs$p_FDR < 0.05,])

na.omit(morph_pgls$coefs[morph_pgls$coefs$p < 0.05,])
na.omit(morph_pgls$coefs[morph_pgls$coefs$p_FDR < 0.05,])

readr::write_csv(col_pgls$coefs, file.path(res_folder, "univariate/colorPC_pgls_coefs.csv"))
readr::write_csv(col_pgls$aic,   file.path(res_folder, "univariate/colorPC_pgls_aic.csv"))

readr::write_csv(morph_pgls$coefs, file.path(res_folder, "univariate/morphoPC_pgls_coefs.csv"))
readr::write_csv(morph_pgls$aic,   file.path(res_folder, "univariate/morphoPC_pgls_aic.csv"))

display_result <- function(df, pcol = "p_FDR", p.value = 0.05, title = "") {
  df %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(sig = !!as.name(pcol) < p.value) %>%
    dplyr::filter(sig) %>%
    dplyr::select(PC, term, estimate, SE, lower, upper, p, p_FDR, best_correlation) %>%
    knitr::kable(digits = 4, caption = title) %>%
    kableExtra::kable_styling(full_width = FALSE)
}

cat("\nSignificant PGLS (FDR<0.05) for colour PCs:\n")
print(display_result(col_pgls$coefs, title = "Significant PGLS (best correlation) — colour PCs"))

cat("\nSignificant PGLS (FDR<0.05) for morphology PCs:\n")
print(display_result(morph_pgls$coefs, title = "Significant PGLS (best correlation) — morphology PCs"))

#---------------------------#
# (i) CorStruct AICc tables #
#---------------------------#
make_corstruct_aic_table <- function(aic_corr_df, outfile_tex,
                                     title = "PGLS correlation-structure AICc (best predictor set per PC)") {
  # Keep only Model (BM/Pagel/OU), PC, AICc
  tab <- aic_corr_df %>%
    dplyr::select(PC, Model, AICc) %>%
    # optional rounding for prettiness
    dplyr::mutate(AICc = round(AICc, 2)) %>%
    tidyr::pivot_wider(names_from = PC, values_from = AICc) %>%
    dplyr::arrange(factor(Model, levels = c("BM","Pagel","OU")))
  
  # Pretty corStruct labels
  tab$Model <- dplyr::recode(tab$Model,
                             "BM" = "BM",
                             "Pagel" = "Pagel's $\\lambda$",
                             "OU" = "OU ($\\alpha$)")
  
  # Write LaTeX
  kable(tab, format = "latex", booktabs = TRUE, escape = FALSE,
        caption = title, label = NULL) %>%
    kable_styling(latex_options = c("hold_position","striped")) %>%
    writeLines(con = outfile_tex)
}

# Colour & Morphology AICc-by-PC tables
make_corstruct_aic_table(col_pgls$aic_corr,
                         outfile_tex = "pgls_color_corstruct_aic_byPC.tex",
                         title = "Colour PCs: PGLS correlation-structure AICc (best predictor set per PC)")

make_corstruct_aic_table(morph_pgls$aic_corr,
                         outfile_tex = "pgls_morpho_corstruct_aic_byPC.tex",
                         title = "Morphology PCs: PGLS correlation-structure AICc (best predictor set per PC)")

#---------------------------------------------------#
# (ii) Best-model summary per PC (with λ/α/— value) #
#---------------------------------------------------#
make_best_model_summary <- function(info_df, aic_corr_df, outfile_tex,
                                    title = "Best PGLS model per PC (predictor set, correlation, parameter)") {
  
  # param name & value by corStruct
  best <- info_df %>%
    mutate(param_name = dplyr::case_when(
      best_correlation == "Pagel" ~ "$\\lambda$",
      best_correlation == "OU"    ~ "$\\alpha$",
      TRUE                        ~ "---"
    ),
    param_value = dplyr::case_when(
      best_correlation == "Pagel" ~ lambda,
      best_correlation == "OU"    ~ alpha,
      TRUE                        ~ NA_real_
    )) %>%
    # Bring in the AICc of the selected corStruct (delta==0)
    left_join(
      aic_corr_df %>%
        group_by(PC) %>%
        slice_min(order_by = delta, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(PC, AICc),
      by = "PC"
    ) %>%
    transmute(
      PC,
      Predictor = predictor_model,
      CorStruct = dplyr::recode(best_correlation,
                                "BM"="BM", "Pagel"="Pagel's $\\lambda$", "OU"="OU"),
      Parameter = param_name,
      Value = ifelse(is.na(param_value), "---", sprintf("%.3f", param_value)),
      AICc = sprintf("%.2f", AICc)
    ) %>%
    arrange(PC)
  
  kable(best, format = "latex", booktabs = TRUE, escape = FALSE,
        caption = title, col.names = c("PC","Predictor set","Correlation","Parameter","Estimate","AICc")) %>%
    kable_styling(latex_options = c("hold_position","striped")) %>%
    writeLines(con = outfile_tex)
}

# ----------------------------
# Build a LaTeX coefficients table from run_pgls_models(...)$coefs
# ----------------------------
make_pgls_coef_table <- function(coef_df,
                                 outfile_tex,
                                 caption = "PGLS coefficients for best models (per PC)",
                                 digits = 3,
                                 show_only_significant = FALSE,     # set TRUE to filter by FDR<0.05
                                 sort_by = c("PC","predictor","term")) {
  
  stopifnot(all(c("PC","predictor","best_correlation","term",
                  "estimate","SE","t","p","p_FDR","lower","upper") %in% names(coef_df)))
  
  tab <- coef_df %>%
    # drop intercepts
    filter(term != "(Intercept)") %>%
    # optional significance filter
    { if (show_only_significant) filter(., !is.na(p_FDR) & p_FDR < 0.05) else . } %>%
    # pretty names / rounding
    dplyr::mutate(
      Correlation = recode(best_correlation,
                           "BM"="BM",
                           "Pagel"="Pagel's $\\lambda$",
                           "OU"="OU ($\\alpha$)"),
      Estimate = round(estimate, digits),
      SE       = round(SE, digits),
      t        = round(t, digits),
      p        = signif(p, 3),
      FDR      = ifelse(is.na(p_FDR), NA, signif(p_FDR, 3)),
      `95% CI` = ifelse(is.na(lower) | is.na(upper),
                        "---",
                        paste0("[", round(lower, digits), ", ", round(upper, digits), "]"))
    ) %>%
    dplyr::select(PC,
           predictor,
           Correlation,
           term,
           Estimate, SE, t, p, FDR, `95% CI`) 

  tab <- tab %>%
    arrange(dplyr::across(all_of(sort_by)))
  
  if (nrow(tab) == 0L) {
    warning("No rows to print (did you filter to only significant?). Writing an empty table.")
  }
  
  kable(tab, format = "latex", booktabs = TRUE, escape = FALSE,
        caption = caption,
        col.names = c("PC", "Predictor set", "Correlation", "Term",
                      "Estimate", "SE", "$t$", "$p$", "FDR", "95\\% CI")) %>%
    kable_styling(latex_options = c("hold_position","striped","repeat_header")) %>%
    collapse_rows(columns = 1:3, latex_hline = "major") %>%   # group by PC / predictor / correlation
    writeLines(con = outfile_tex)
}

# ----------------------------
# Produce tables for colour and morphology
# ----------------------------

# Colour coefficients (all terms)
make_pgls_coef_table(col_pgls$coefs,
                     outfile_tex = "pgls_color_coefficients.tex",
                     caption = "Colour PCs: coefficients of best PGLS models (per PC)",
                     digits = 3,
                     show_only_significant = FALSE)

# Colour coefficients (FDR-significant only) – optional extra
make_pgls_coef_table(col_pgls$coefs,
                     outfile_tex = "pgls_color_coefficients_sig.tex",
                     caption = "Colour PCs: FDR-significant coefficients (best PGLS per PC)",
                     digits = 3,
                     show_only_significant = TRUE)

# Morphology coefficients (all terms)
make_pgls_coef_table(morph_pgls$coefs,
                     outfile_tex = "pgls_morpho_coefficients.tex",
                     caption = "Morphology PCs: coefficients of best PGLS models (per PC)",
                     digits = 3,
                     show_only_significant = FALSE)

# Morphology coefficients (FDR-significant only) – optional extra
make_pgls_coef_table(morph_pgls$coefs,
                     outfile_tex = "pgls_morpho_coefficients_sig.tex",
                     caption = "Morphology PCs: FDR-significant coefficients (best PGLS per PC)",
                     digits = 3,
                     show_only_significant = TRUE)
