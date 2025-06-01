interactionHelper <- function(
    data,
    outcome_type = c("bt", "qt"),
    formula_str,
    prefix       = "results",
    use_firth    = FALSE,
    firth_engine = c("brglm2", "logistf"),
    nonlinear_E  = NULL
) {
  ## ---- package check ----------------------------------------------------
  pkgs <- c("ggplot2", "sandwich", "lmtest", "ggeffects", "rlang", "car")
  if (use_firth) pkgs <- c(pkgs, match.arg(firth_engine))
  miss <- setdiff(pkgs, rownames(installed.packages()))
  if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
  
  ## ---- parse formula & optionally add E^2 terms -------------------------
  fml_chr <- formula_str
  if (!is.null(nonlinear_E)) {
    for (v in nonlinear_E) {
      if (is.numeric(data[[v]]) && grepl(paste0("\\b", v, "\\b"), fml_chr)) {
        quad <- paste0("I(", v, "^2)")
        if (!grepl(quad, fml_chr, fixed = TRUE))
          fml_chr <- paste(fml_chr, "+", quad)
      }
    }
  }
  fml <- stats::as.formula(fml_chr)
  
  ## ---- fit model --------------------------------------------------------
  outcome_type <- match.arg(outcome_type)
  firth_engine <- match.arg(firth_engine)
  
  if (outcome_type == "bt" && use_firth) {
    if (firth_engine == "brglm2") {
      library(brglm2)
      fit <- stats::glm(fml, data = data,
                        family = binomial("logit"),
                        method  = "brglmFit",
                        type    = "AS_mean")    # Firth/adjusted-score
      vcov_use <- vcov(fit)                           # bias-reduced vcov
      rob_tab  <- lmtest::coeftest(fit, vcov. = vcov_use)
      lh_test  <- "LR"                                # profile-LR for contrasts
    } else {                                          # logistf
      fit <- logistf::logistf(fml, data = data)
      vcov_use <- fit$var
      est  <- coef(fit)
      se   <- sqrt(diag(vcov_use))
      rob_tab <- cbind(Estimate = est,
                       `Std. Error` = se,
                       `z value`    = est / se,
                       `Pr(>|z|)`   = 2*pnorm(-abs(est/se)))
      rownames(rob_tab) <- names(est)
      lh_test <- "LR"                                 # use LR for contrasts
    }
  } else {                                           # ordinary GLM / LM
    fit <- if (outcome_type == "bt")
      stats::glm(fml, data = data, family = binomial)
    else
      stats::lm (fml, data = data)
    vcov_use <- sandwich::vcovHC(fit, type = "HC3")
    rob_tab  <- lmtest::coeftest(fit, vcov. = vcov_use)
    lh_test  <- if (inherits(fit, "glm")) "Chisq" else "F"
  }
  
  ## ---- interaction coefficient names -----------------------------------
  int_rows <- grep(":", rownames(rob_tab))
  if (!length(int_rows)) stop("No interaction coefficients found.")
  int_coef_names <- rownames(rob_tab)[int_rows]
  
  ## ---- hypothesis tests -------------------------------------------------
  extract_lh <- function(R_mat, label) {
    q <- nrow(R_mat)
    beta <- coef(fit)
    W    <- R_mat %*% beta
    V_R  <- R_mat %*% vcov_use %*% t(R_mat)
    if (lh_test %in% c("Chisq", "LR")) {             # χ² or LR
      stat <- as.numeric(t(W) %*% solve(V_R, W))
      pval <- pchisq(stat, df = q, lower.tail = FALSE)
      data.frame(test_name = label, df1 = q, df2 = NA, statistic = stat,
                 p_value = pval, test_type = "Chisq",
                 contrast = label, row.names = NULL)
    } else {                                         # F
      df2  <- fit$df.residual
      stat <- (t(W) %*% solve(V_R, W)) / q
      pval <- pf(stat, q, df2, lower.tail = FALSE)
      data.frame(test_name = label, df1 = q, df2 = df2, statistic = stat,
                 p_value = pval, test_type = "F",
                 contrast = label, row.names = NULL)
    }
  }
  
  ## ---- global q-df test ---------------------------------------------------
  p <- length(coef(fit))
  q <- length(int_rows)
  R_global <- matrix(0, q, p, dimnames = list(NULL, names(coef(fit))))
  R_global[cbind(seq_len(q), int_rows)] <- 1
  lh_results <- extract_lh(R_global, "Global_interaction")
  
  ## pairwise contrasts ------------------------------------------------------
  q <- length(int_rows)
  if (q >= 2) {
    for (pair in combn(seq_len(q), 2, simplify = FALSE)) {   # <-- line changed
      R <- matrix(0, 1, length(coef(fit)),
                  dimnames = list(NULL, names(coef(fit))))
      R[1, int_rows[pair[1]]] <-  1                         # <-- map back
      R[1, int_rows[pair[2]]] <- -1
      lbl <- paste("Contrast",
                   int_coef_names[pair[1]], "vs",
                   int_coef_names[pair[2]], sep = "_")
      lh_results <- rbind(lh_results, extract_lh(R, lbl))
    }
  }
  
  ## ---- loop over interaction terms for plots ---------------------------
  term_labels <- grep(":", attr(stats::terms(fml), "term.labels"), value = TRUE)
  p_col <- grep("^Pr", colnames(rob_tab))
  plot_prob_list <- list(); plot_rd_list <- list()
  
  for (k in seq_along(term_labels)) {
    term_lbl <- term_labels[k]
    vars     <- strsplit(term_lbl, ":")[[1]]
    is_num   <- vapply(data[vars], is.numeric, logical(1))
    
    if (any(is_num) && any(!is_num)) {
      cont   <- vars[ is_num ][1]
      catg   <- vars[!is_num][1]
      gp_terms <- c(paste0(cont, " [all]"), catg)
    } else {
      gp_terms <- c(paste0(vars[1], " [all]"), vars[2])
      catg <- NULL
    }
    
    ## prediction grid
    grid <- ggeffects::ggpredict(fit, terms = gp_terms)
    
    ## first interaction-row p-value
    coef_row <- grep(paste0("^", vars[1], ".*:", vars[2], "|",
                            "^", vars[2], ".*:", vars[1]),
                     int_coef_names, value = TRUE)[1]
    int_p <- rob_tab[coef_row, p_col]
    
    ## probability / value plot
    pl_prob <- ggplot2::ggplot(grid,
                               ggplot2::aes(x = x, y = predicted, colour = group)) +
      ggplot2::geom_line(size = 0.9) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = conf.low, ymax = conf.high, fill = group),
        alpha = 0.15, colour = NA) +
      ggplot2::labs(
        x        = attr(grid, "x.title"),
        y        = if (outcome_type == "bt") "Predicted probability"
        else "Predicted value",
        colour   = attr(grid, "legend.title"),
        fill     = attr(grid, "legend.title"),
        subtitle = sprintf("%s  (p = %.3g)", coef_row, int_p)) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.key   = ggplot2::element_blank())
    
    ggplot2::ggsave(paste0(prefix, ".prob.", k, ".pdf"),
                    pl_prob, width = 7, height = 5)
    plot_prob_list[[k]] <- pl_prob
    
    ## risk difference (binary + categorical)
    pl_rd <- NULL
    if (outcome_type == "bt" && !is.null(catg)) {
      ref_level <- levels(data[[catg]])[1]
      ref_grid  <- subset(grid, group == ref_level)
      rd <- merge(grid, ref_grid, by = "x", suffixes = c("", ".ref"))
      
      se  <- (rd$conf.high - rd$conf.low) / (2*1.96)
      seR <- (rd$conf.high.ref - rd$conf.low.ref) / (2*1.96)
      
      rd$diff    <- rd$predicted - rd$predicted.ref
      rd$se_diff <- sqrt(se^2 + seR^2)
      rd$low     <- rd$diff - 1.96*rd$se_diff
      rd$high    <- rd$diff + 1.96*rd$se_diff
      
      pl_rd <- ggplot2::ggplot(rd,
                               ggplot2::aes(x = x, y = diff, colour = group, fill = group)) +
        ggplot2::geom_line(size = 0.9) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = low, ymax = high),
          alpha = 0.15, colour = NA) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::labs(
          x        = attr(grid, "x.title"),
          y        = "Risk difference vs reference",
          colour   = attr(grid, "legend.title"),
          fill     = attr(grid, "legend.title"),
          subtitle = sprintf("%s  (p = %.3g)", coef_row, int_p)) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       legend.key   = ggplot2::element_blank())
      
      ggplot2::ggsave(paste0(prefix, ".riskdiff.", k, ".pdf"),
                      pl_rd, width = 7, height = 5)
    }
    plot_rd_list[[k]] <- pl_rd
  }
  
  ## ---- write tables -----------------------------------------------------
  utils::write.table(rob_tab,
                     paste0(prefix, ".summary.tsv"),
                     sep = "\t", quote = FALSE, col.names = NA)
  
  utils::write.table(lh_results,
                     paste0(prefix, ".lh.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  
  invisible(list(model         = fit,
                 robust_table  = rob_tab,
                 lh_table      = lh_results,
                 plot_prob     = plot_prob_list,
                 plot_riskdiff = plot_rd_list))
}
