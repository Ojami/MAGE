# ============================================================================
# fitSurvCox.R  –  General‑purpose survival wrapper (rev‑6, Aug 2025)
# ----------------------------------------------------------------------------
#  * Supports alternate time‑scales (TTE / AOO) with optional left truncation.
#  * Handles one time‑varying covariate via counting‑process expansion.
#  * Loops over a predictor set (e.g. GWAS SNPs) and fits:
#      – univariable Cox  (HR)
#      – multivariable Cox (aHR)   [optionally with interaction term]
#      – sub‑distribution models (sHR) via crr / fastcmprsk
#  * Generates KM / log‑rank plots, CIF ribbons, PH diagnostics, overall KM & CIF.
#  * Computes median follow‑up by reverse Kaplan–Meier.
# ----------------------------------------------------------------------------

# ----- Imports ---------------------------------------------------------------
library(survival)
library(survminer)
library(ggsurvfit)
library(ggplot2)
library(dplyr)
library(data.table)
library(prodlim)
library(pbapply)
library(parallel)
library(riskRegression)
library(splines) 

# ----- Helper: pruning of model metric ----------------------------------------
.get_metrics <- function(fit, label, which_model) {
  g       <- broom::glance(fit)
  df_coef <- length(coef(fit))
  
  data.frame(
    Predictor   = label,
    Model       = which_model,
    N           = g$n,
    N_events    = g$nevent,
    Concordance = round(g$concordance, 3),
    LRT         = sprintf("%.3f (df = %d, p = %.3f)",
                          g$statistic.log, df_coef, g$p.value.log),
    check.names = FALSE
  )
}

# ----- Helper: pruning of model summaries -------------------------------------
.prune <- function(fit, tag) {
  s   <- summary(fit, digits = 3)
  out <- cbind(as.data.frame(s$coef), as.data.frame(s$conf.int))
  
  # rename columns
  nm <- names(out)
  nm[nm == "exp(coef)"]  <- tag          # HR / aHR
  nm[nm == "se(coef)"]   <- "SE"         # <‑‑ keep standard error
  nm[nm == "Pr(>|z|)"]   <- "P"
  nm[nm == "lower .95"]  <- "CIL"
  nm[nm == "upper .95"]  <- "CIH"
  names(out) <- nm
  
  out$CI   <- sprintf("[%.2f, %.2f]", out$CIL, out$CIH)
  out$Term <- rownames(out)
  
  out[, c("Term", tag, "SE", "P", "CI")]
  
}

# ----- Helper: merge a single TVC ---------------------------------------------
.merge_tvc <- function(df, tv_df, entry, exit, default = 0) {
  
  ## ----------------------------------------------------------
  ##  tv_df columns (required): eid | start | stop | value
  ##  Optional   tv_df$name   : variable label
  ##  default outside any interval : 0 for binary, NA for cont.
  ## ----------------------------------------------------------
  
  vname <- if ("name" %in% names(tv_df)) unique(tv_df$name)[1] else "TVC"

  ## 1. guarantee 'value' column (ever‑treated shortcut)
  if (!"value" %in% names(tv_df)) tv_df$value <- 1
  tv_df <- tv_df[, c("eid", "start", "stop", "value")]

  
  ## 2. explode intervals → change‑points ---------------------
  long_on  <- tv_df[, c("eid", "start", "value")]
  names(long_on) <- c("eid", "time", "value")
  
  long_off <- tv_df[, c("eid", "stop")]
  names(long_off) <- c("eid", "time")
  long_off$value  <- default                 # reset after interval
  
  changes <- rbind(long_on, long_off)
  changes <- changes[order(changes$eid, changes$time), ]
  changes <- changes[!duplicated(changes[, c("eid","time")]), ]
  
  
  ## 3. base counting‑process frame ---------------------------
  base <- df[, c("eid", entry, exit)]
  names(base) <- c("eid", "tstart", "tstop")
  idx_zero <- base$tstart >= base$tstop
  if (any(idx_zero)) base$tstop[idx_zero] <- base$tstop[idx_zero] + 1e-6
  
  ## 4. tmerge ------------------------------------------------
  tdat <- survival::tmerge(base, base,
                           id     = eid,
                           tstart = tstart,
                           tstop  = tstop)
  
  tdat <- survival::tmerge(tdat, changes,
                           id      = eid,
                           "..tmp" = tdc(time, value, init = default))
  
  tol <- sqrt(.Machine$double.eps)
  len <- tdat$tstop - tdat$tstart
  tdat <- tdat[len > tol, ]
  
  names(tdat)[names(tdat) == "..tmp"] <- vname
  dplyr::left_join(tdat, df, by = "eid")
}

# ----- Helper: delayed‑entry CIF with bootstrap --------------------------------
fit_cif <- function(dat, 
                    genoCol,
                    timeStop,
                    entryCol,
                    statusCol,
                    strataDF,
                    times,
                    causeLabel){
  
  form <- as.formula(
    sprintf("Hist(%s, %s, entry = %s) ~ %s",
            timeStop, statusCol, entryCol, genoCol))
  fit <- prodlim::prodlim(form, data = dat)
  pred <- predict(fit, newdata = strataDF, times = times,
                  cause = causeLabel, type = "risk")
  mat  <- do.call(cbind, lapply(pred, as.numeric))
  if (is.vector(mat)) mat <- matrix(mat, ncol = 1L)
  mat}


cif_boot_ci <- function(df, genoCol,
                        timeStop, entryCol, statusCol,
                        causeLabel = "event",
                        times   = NULL,
                        nBoot   = 500,
                        seed    = 1,
                        parallel = FALSE,
                        nCores   = NULL) {
  set.seed(seed)
  if (is.null(times)) {
    tmpT <- data.frame(
      time  = df[[timeStop]],
      event = df[[statusCol]],
      entry = df[[entryCol]])
    times <- sort(unique(prodlim::prodlim(
      Hist(time, event, entry = entry) ~ 1, data = tmpT)$time))
  }
  
  strataDF <- unique(df[, genoCol, drop = FALSE])
  pointMat <- fit_cif(df, genoCol,
                      timeStop, entryCol, statusCol,
                      strataDF, times, causeLabel)
  bootFun <- function(idx) {
    fit_cif(df[idx, , drop = FALSE], genoCol,
            timeStop, entryCol, statusCol,
            strataDF, times, causeLabel)
  }
  idxMat <- matrix(
    sample.int(nrow(df), nrow(df) * nBoot, replace = TRUE),
    nrow = nBoot)
  if (parallel) {
    if (is.null(nCores))
      nCores <- max(1L, parallel::detectCores() - 1L)
    if (.Platform$OS.type == "windows")
      nCores <- min(nCores, 60L)
    clType <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"
    cl <- parallel::makeCluster(nCores, type = clType)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE))
    ## ensure required packages are on each worker (PSOCK needs this)
    parallel::clusterEvalQ(cl, { library(prodlim); library(data.table) })
    parallel::clusterExport(
      cl,
      c("df", "genoCol", "timeStop", "entryCol", "statusCol",
        "strataDF", "times", "causeLabel", "fit_cif"),
      envir = environment())
    bootMats <- pbapply::pblapply(
      X   = seq_len(nBoot),
      FUN = function(b) bootFun(idxMat[b, ]),
      cl  = cl)
  } else {
    bootMats <- pbapply::pblapply(
      X   = seq_len(nBoot),
      FUN = function(b) bootFun(idxMat[b, ]))
  }
  arr <- array(unlist(bootMats),
               dim = c(length(times), ncol(pointMat), nBoot))
  lowMat  <- apply(arr, 1:2, quantile, .025, na.rm = TRUE)
  highMat <- apply(arr, 1:2, quantile, .975, na.rm = TRUE)
  kStrata    <- ncol(pointMat)
  labsStrata <- apply(strataDF, 1, paste, collapse = "_")
  data.table::data.table(
    time     = rep(times, kStrata),
    strata   = factor(rep(labsStrata, each = length(times)),
                      levels = labsStrata),
    cif      = as.numeric(pointMat),
    cif.low  = as.numeric(lowMat),
    cif.high = as.numeric(highMat)
  )
}

# ==============================================================================
fitSurvCox <- function(
    ## ----- CORE DATA --------------------------------------------------------
    df,                         # main data.frame
    outcomeVar,                 # column: 0 = censored, 1 = event
    competeVar   = NULL,        # multi‑level status; defaults to outcomeVar
    preds,                      # character vector of predictors (analysed 1‑by‑1)
    covars       = NULL,        # baseline covariates
    tv_df        = NULL,        # optional time‑varying‑covariate data.frame
    interactionVar = NULL,      # OPTIONAL: single column for interaction (e.g. "Sex")
    
    ## ----- TIME‑SCALE OPTIONS ----------------------------------------------
    timeOrigin = c("TTE", "AOO"),   # choose time‑scale logic
    LT         = FALSE,             # TRUE = left‑truncation / delayed entry
    prevalent  = FALSE,             # drop prevalent cases when FALSE
    landmark   = 0,                 # remove events ≤ landmark years
    
    ## ----- COMPETING‑RISK MODEL --------------------------------------------
    FGR   = c("none", "crr", "fastcrr"),  # sub‑distribution HR method
    nBoot = 5L,                         # bootstrap reps for CIF ribbons
    
    ## ----- PLOTTING TOGGLES & STYLE ----------------------------------------
    overallKM      = TRUE,          # KM curve for entire cohort
    plotCIF     = TRUE,          # overall/predictor specific CIF (Aalen–Johansen)
    logRank        = TRUE,          # per‑predictor KM / log‑rank plot
    logRankFun     = c("null","event","cumhaz","pct"), # KM transformation
    levelCutoff    = 5,             # max unique levels to stratify predictor
    kmEventYlimPad = 0.02,          # padding for 1‑S(t) y‑axis
    diag           = TRUE,          # proportional‑hazards diagnostics
    
    ## ----- I/O, REPRODUCIBILITY --------------------------------------------
    outDir = getwd(),               # folder for PDFs & results
    seed   = 1L,                    # RNG seed
    parallel = FALSE,
    nCores   = 5,
    
    ## arguments for nonlinearity and time‐dependency ---------------------
    nonlinVars   = NULL,        # e.g. c("Age") to model with splines
    splineDF     = 3,           # degrees of freedom for spline(s)
    timeDepVars  = NULL,        # e.g. c("Age") for time-dependent effects
    timeDepFun   = function(x, t, ...) x * log(t) # function for time transformation
){
  
  # Prep & checks --------------------------------------------------------------
  stopifnot(is.data.frame(df), is.character(preds), outcomeVar %in% names(df))
  if (is.null(competeVar)) competeVar <- outcomeVar
  stopifnot(competeVar %in% names(df))
  
  timeOrigin <- match.arg(timeOrigin)
  FGR        <- match.arg(FGR)
  logRankFun <- match.arg(logRankFun)
  
  if (LT) logRank <- FALSE # logRank is not supported with LT 
  
  if (!is.null(covars) && !is.character(covars))
    stop("`covars` must be character vector or NULL")
  if (!is.null(tv_df) && is.null(df$eid))
    stop("`df` must contain an `eid` column when a time‑varying covariate is supplied")
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  set.seed(seed)
  
  
  # Left‑truncation alignment ---------------------------------------------------
  if (LT) {                                             # delayed entry
    if (timeOrigin == "TTE") {
      entry_col <- "tt0"            # start  = age at diagnosis
      exit_col  <- "tt"             # stop   = time to event
    } else {                        # AOO + LT
      entry_col <- "tt"             # start  = time since birth to enrol
      exit_col  <- "tt0"            # stop   = age at event
      df[[entry_col]] <- df[[exit_col]] - df[[entry_col]]  # W = Y‑X
    }
  } else {                                              # no LT
    ## create a synthetic zero column so downstream code
    ## (padding, tmerge, CIF bootstrap) always has a start variable
    entry_col <- "entry0"
    if (!"entry0" %in% names(df))
      df$entry0 <- 0
    exit_col  <- "tt"
  }
  
  # Remove prevalent / invalid follow-up --------------------------------------
  if (!prevalent && "base" %in% names(df)) {
    df <- df[df$base <= 0, ]
    df$base <- NULL
  }
  # drop any negative or zero follow-up
  if (timeOrigin == "TTE")  idx_neg <- df[[exit_col]] <  0 else
    idx_neg <- df[[exit_col]] <= 0
  if (any(idx_neg)) df <- df[!idx_neg, ]
  
  # Landmark filter (now we know exactly which ‘clock’ to use) -----------------
  if (landmark > 0) {
    # drop anyone whose event occurred at or before the landmark
    # note: you could also drop those who are censored before L, if you like
    idx_landmark <- df[[exit_col]] <= landmark & df[[outcomeVar]] == 1
    df <- df[!idx_landmark, ]
  }

  
  # Merge TVC ------------------------------------------------------------------
  if (!is.null(tv_df)) {
    df <- .merge_tvc(df, tv_df, entry = entry_col, exit = exit_col)
    entry_col <- "tstart"; exit_col <- "tstop"
  } else {
    ## Pad zero‑length intervals  (exit == entry)
    idx_zero <- df[[entry_col]] >= df[[exit_col]]
    if (any(idx_zero))
      df[[exit_col]][idx_zero] <- df[[exit_col]][idx_zero] + 1e-6
  }
  
  # Reverse‑KM median follow‑up -------------------------------------------------
  if (LT) {
    surv_mf <- Surv(df[[entry_col]],
                    df[[exit_col]],
                    df[[competeVar]] == 0,
                    type = "counting")
  } else {
    surv_mf <- Surv(df[[exit_col]],
                    df[[competeVar]] == 0)   # two‑column form
  }
  
  mf <- survminer::surv_median(survival::survfit(surv_mf ~ 1))
  
  
  # overall KM plots -----------------------------------------------------------
  if (overallKM) {
    
    if (LT) {
      surv_overall <- survfit2(
        Surv(df[[entry_col]],              # start
             df[[exit_col]],               # stop
             df[[outcomeVar]] == 1,        # status
             type = "counting") ~ 1,
        data = df)
    } else {                               # no left‑truncation
      surv_overall <- survfit2(
        Surv(df[[exit_col]],               # time
             df[[outcomeVar]] == 1) ~ 1,
        data = df)
    }
    
    p <- ggsurvfit::ggsurvfit(surv_overall) +
      ggsurvfit::add_confidence_interval() +
      labs(x = "Years", y = "Overall survival probability") +
      ggsurvfit::add_risktable()
    
    fname <- sprintf("KM.%s.pdf", make.names(outcomeVar))
    grDevices::cairo_pdf(file.path(outDir, fname), onefile = FALSE)
    print(p)
    dev.off()
  }
  
  # overall CIF plot -----------------------------------------------------------
  if (plotCIF) {
    
    fac <- factor(df[[competeVar]],
                  levels = 0:3,
                  labels = c("censor", "event", "death", "competing"))
    df$..statusOverall <- droplevels(fac)
    
    if (LT) {
      ## ----- delayed entry – bootstrap helper ---------------------------
      keep <- complete.cases(df[, c(entry_col, exit_col, "..statusOverall")])
      tmp  <- as.data.frame(df[keep, , drop = FALSE])
      
      tmp$entry__  <- tmp[[entry_col]]
      tmp$exit__   <- tmp[[exit_col]]
      tmp$status__ <- tmp$..statusOverall
      tmp$all__    <- factor("overall")
      
      cif_dt <- cif_boot_ci(
        df         = tmp,
        genoCol    = "all__",
        timeStop   = "exit__",
        entryCol   = "entry__",
        statusCol  = "status__",
        causeLabel = "event",
        nBoot      = nBoot,
        seed       = seed,
        parallel   = parallel,
        nCores     = nCores)
      
      p <- ggplot(cif_dt,
                  aes(time, cif, colour = strata, fill = strata)) +
        geom_step(linewidth = .7) +
        geom_ribbon(aes(ymin = cif.low, ymax = cif.high),
                    alpha = .25, linewidth = 0) +
        theme_light() +
        labs(x = "Years", y = "Cumulative incidence") +
        guides(colour = "none", fill = "none")      # one stratum
      
    } else {
      ## ----- right‑censored – tidycmprsk + ggcuminc ---------------------
      keep <- complete.cases(df[, exit_col])
      cm   <- tidycmprsk::cuminc(
        Surv(df[[exit_col]][keep], df$..statusOverall[keep]) ~ 1,
        data = df[keep, , drop = FALSE])
      
      p <- ggcuminc(cm, outcome = setdiff(levels(df$..statusOverall), "censor")) +
        labs(x = "Years") +
        add_confidence_interval() +                 # <- ribbon
        add_risktable() +
        theme_light() +
        guides(colour = "none", fill = "none")
    }
    
    ## write with Cairo for crisp ribbon
    fname <- sprintf("CIF.%s.pdf", make.names(outcomeVar))
    grDevices::cairo_pdf(file.path(outDir, fname), onefile = FALSE)
    print(p)
    dev.off()
  }
  

  # ---- Run analysis per predictor -----------------------------------------
  res_list     <- vector("list", length(preds))
  metrics_list <- vector("list", length(preds))
  
  for (pred in preds) {
    
    ## ---------------------------------------------------------------
    ## 1 │ SURVIVAL OBJECTS (handles LT automatically)
    ## ---------------------------------------------------------------
    if (LT) {
      surv_cs <- Surv(df[[entry_col]], df[[exit_col]], df[[outcomeVar]] == 1, type = "counting")
      surv_cr <- Surv(df[[entry_col]], df[[exit_col]], df[[competeVar]])
    } else {
      surv_cs <- Surv(df[[exit_col]], df[[outcomeVar]] == 1)
      surv_cr <- Surv(df[[exit_col]], df[[competeVar]])
    }
    
    ## ------------------------------------------------------------------
    ## 2 │ CATEGORICAL KM / LOG‑RANK  (skipped when LT or too many levels)
    ## ------------------------------------------------------------------
    nlev <- data.table::uniqueN(df[[pred]], na.rm = TRUE)
    
    if (logRank && !LT && nlev <= levelCutoff) {
      
      # --- build a formula that references real columns -----------------
      fstr  <- sprintf("Surv(%s, `%s`==1) ~ `%s`",
                       exit_col, outcomeVar, pred)   # literal string
      
      form_km <- stats::as.formula(fstr)             # <- real formula object
      km_fit  <- survfit2(form_km, data = df)
      
      km_fit$call$formula <- form_km                 # <-- overwrite the call
      
      # median line
      sfit <- summary(km_fit)
      median_line <- ifelse(all(is.na(sfit$table[, "median"])), "none", "hv")
      
      # transformation
      fun_opt <- match.arg(logRankFun)
      fun     <- if (fun_opt == "null") NULL else fun_opt
      
      ylim <- NULL
      if (identical(fun, "event"))
        ylim <- c(0, 1 - min(km_fit$surv, na.rm = TRUE) + kmEventYlimPad)
      
      # x‑axis span
      tcols  <- names(df)[startsWith(names(df), "tt")]
      axlims <- range(unlist(df[, tcols]), na.rm = TRUE)
      if (timeOrigin == "AOO") axlims[2] <- axlims[2] + 2
      else                     axlims    <- c(0, axlims[2] + 1)
      
      g <- survminer::ggsurvplot(
        km_fit,
        fun             = fun,
        data            = df,
        conf.int        = TRUE,
        risk.table      = "absolute",
        xlab            = "Years",
        legend.title    = pred,
        ggtheme         = ggplot2::theme_bw(),
        break.time.by   = 10,
        ylim            = ylim,
        xlim            = axlims,
        pval            = TRUE,
        pval.coord      = c(1, mean(if (is.null(ylim)) c(0,1) else ylim)),
        surv.median.line = median_line,
        censor          = FALSE,
        size            = .5
      )

      fname <- sprintf("LogRankPlot.%s.%s.pdf",
                       make.names(outcomeVar), make.names(pred))
      grDevices::cairo_pdf(file.path(outDir, fname), onefile = FALSE)
      print(g); dev.off()

    }
    
    
    ## ---------------------------------------------------------------
    ## 3 │ BUILD MODEL FORMULAE (incorporating nonlinearity and time‐dependency)
    ## ---------------------------------------------------------------
    # if (is.null(interactionVar)) {
    #   expl <- c(pred, covars)
    # } else {
    #   pred_term <- paste0(pred, "*", interactionVar)
    #   expl <- c(pred_term, setdiff(covars, interactionVar))
    # }
    
    # base term list  (with or without interaction) ------------
    if (is.null(interactionVar)) {
      base_terms <- c(pred, covars)
    } else {
      base_terms <- c(pred,
                      interactionVar,
                      paste0(pred, ":", interactionVar),
                      setdiff(covars, interactionVar))
    }
    
    # baseline transformation (non-linear shape) ---------------
    baseline_term <- function(v) {
      if (v %in% nonlinVars)
        sprintf("ns(%s, df = %d)", v, splineDF)
      else
        v
    }
    expl_fixed <- vapply(base_terms, baseline_term, FUN.VALUE = "", USE.NAMES = FALSE)
    
    # add separate tt( ) terms  (RAW variable names only) ------
    tt_terms <- if (length(timeDepVars))
      sprintf("tt(%s)", timeDepVars) else character(0)
    
    # final explanatory set  -----------------------------------
    expl_all      <- unique(c(expl_fixed, tt_terms))
    formula_str   <- paste("surv_cs ~", paste(expl_all, collapse = " + "))
    model_formula <- as.formula(formula_str)
    
    # fit the Cox model  ---------------------------------------
    if (length(timeDepVars)) {
      ## you decide what g(x,t) is when calling the wrapper; e.g.
      ## timeDepFun = function(x, t, ...) x * log(t)
      multi <- coxph(model_formula, data = df, tt = timeDepFun)
    } else {
      multi <- coxph(model_formula, data = df)
    }
    
    # Univariable Cox remains as before:
    uni <- coxph(as.formula(paste("surv_cs ~", pred)), data = df)
    # 
    # # Univariable and multivariable Cox
    # uni   <- coxph(as.formula(paste("surv_cs ~", pred)), data = df)
    # multi <- coxph(reformulate(expl, response = "surv_cs"), data = df)
    
    tab0 <- .prune(uni,  "HR")
    tab1 <- .prune(multi, "aHR")
    out  <- merge(tab0, tab1, by = "Term", all = TRUE)
    
    # model‑level metrics 
    met_uni   <- .get_metrics(uni,   pred, "uni")
    met_multi <- .get_metrics(multi, pred, "multi")
    metrics_list[[pred]] <- rbind(met_uni, met_multi)
    
    ## ---------------------------------------------------------------
    ## 4 │ SUB‑DISTRIBUTION (Fine–Gray) if requested
    ## ---------------------------------------------------------------
    if (FGR != "none") {
      if (FGR == "crr") {
        formula_fg <- as.formula(paste("surv_cr ~", paste(expl, collapse = "+")))
        crr_fit    <- riskRegression::crrmulti(formula_fg, data = df)
      } else {
        formula_fg <- reformulate(expl,
                                  response = sprintf("Crisk(%s, %s)", exit_col, competeVar))
        crr_fit    <- fastcmprsk::fastCrr(formula_fg, data = df)
      }
      tab2 <- .prune(crr_fit, "sHR")
      out  <- merge(out, tab2, by = "Term", all = TRUE)
    }
    
    out$Predictor <- pred
    res_list[[pred]] <- out
    
    ## ---------------------------------------------------------------
    ## 5 │ CIF (bootstrap if LT, else tidycmprsk)
    ## ---------------------------------------------------------------
    if (plotCIF) {
      useBoot <- LT
      status_fac <- factor(df[[competeVar]], levels = 0:3,
                           labels = c("censor","event","death","competing"))
      
      if (useBoot) {
        df$..status <- droplevels(status_fac)
        df$..entry  <- df[[entry_col]]           # no shift
        df$..exit   <- df[[exit_col]]
        
        cif_dt <- cif_boot_ci(
          df         = df,
          genoCol    = pred,
          timeStop   = "..exit",
          entryCol   = "..entry",
          statusCol  = "..status",
          causeLabel = "event",
          nBoot      = nBoot,
          seed       = seed,
          parallel   = parallel,
          nCores     = nCores)
        
        p_cif <- ggplot(cif_dt,
                        aes(time, cif, colour = strata, fill = strata)) +
          geom_step(linewidth = .7) +
          geom_ribbon(aes(ymin = cif.low, ymax = cif.high),
                      alpha = .25, linewidth = 0) +
          theme_light() +
          labs(x = "Years",
               y = "Cumulative incidence",
               colour = pred)
        
      } else {
        keep <- complete.cases(df[, c(exit_col, pred)])
        tmp  <- data.frame(exit   = df[[exit_col]][keep],
                           status = status_fac[keep],
                           geno   = df[[pred]][keep])
        
        cm <- cmprsk::cuminc(tmp$exit, tmp$status, group = tmp$geno)  # cmprsk
        
        p_cif <- ggcuminc(cm, outcome = "event") +
          labs(x = "Years") +
          add_confidence_interval() +
          add_risktable() +
          ggplot2::theme_light() +
          ggplot2::labs(colour = pred)
      }
      
      axlims <- range(df[[exit_col]])
      if (timeOrigin == "AOO") axlims[2] <- axlims[2] + 2
      else                     axlims    <- c(0, axlims[2] + 1)
      
      p_cif <- p_cif + ggplot2::scale_x_continuous(limits = axlims, n.breaks = 6)
      
      fname <- sprintf("CIF.%s.%s.pdf",
                       make.names(outcomeVar), make.names(pred))
      
      grDevices::cairo_pdf(file.path(outDir, fname), onefile = FALSE)
      print(p_cif)
      dev.off()

    }
    
    ## ---------------------------------------------------------------
    ## 6 │ PH diagnostics  — Cairo PDF, one panel per term, no extras
    ## ---------------------------------------------------------------
    # if (diag) {
    #   cox.z  <- survival::cox.zph(multi)
    #   terms  <- setdiff(rownames(cox.z$table), "GLOBAL")   # exclude GLOBAL row
    #   
    #   ## 1 ─ build all term plots first
    #   plot_list <- lapply(terms, function(trm) survminer::ggcoxzph(cox.z, var = trm))
    #   
    #   ## 2 ─ open device after list is ready (avoids blank first page)
    #   fname <- sprintf("PHdiag.%s.%s.pdf",
    #                    make.names(outcomeVar), make.names(pred))
    #   grDevices::cairo_pdf(file.path(outDir, fname), onefile = TRUE)
    #   
    #   ## 3 ─ print each plot to a separate page
    #   for (p in plot_list) print(p)
    #   
    #   dev.off()
    # }
    
    # ---------------------------------------------------------------
    # 6 │ PH diagnostics — Cairo PDF, one panel per term, no extras
    # ---------------------------------------------------------------
    if (diag) {
      # Check if time-dependent variables are used, i.e. whether the model
      # contains any tt() terms. If so, we cannot apply cox.zph() on the original model.
      if (!is.null(timeDepVars) && length(timeDepVars) > 0) {
        
        ## helper: always get a ggplot from ggcoxzph()
        as_ggplot <- function(g) {
          if (inherits(g, "ggplot"))               return(g)
          if (is.list(g) && inherits(g[[1]], "ggplot")) return(g[[1]])
          stop("No ggplot object in ggcoxzph() output")
        }
        
        ## 1 ─ diagnostic model (tt() stripped)
        diag_formula <- gsub("tt\\(([^)]+)\\)", "\\1", formula_str)
        diag_model   <- coxph(as.formula(diag_formula), data = df)
        
        ## 2 ─ PH test on KM-rank scale
        zp     <- survival::cox.zph(diag_model, transform = "km")
        terms  <- setdiff(rownames(zp$table), "GLOBAL")
        x_grid <- zp$x                       # same grid ggcoxzph() uses
        
        ## 3 ─ open PDF
        if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
        if (length(dev.list()))  graphics.off()
        pdf(file.path(outDir,
                      sprintf("PHdiag.%s.%s.pdf",
                              make.names(outcomeVar), make.names(pred))),
            onefile = TRUE)
        
        ## 4 ─ loop over covariates
        for (trm in terms) {
          
          ## 4a – base panel (blue dots, LOESS, p-value)
          gp_base <- as_ggplot(survminer::ggcoxzph(zp, var = trm))
          
          ## 4b – coefficients
          beta_const <- coef(multi)[trm]                      # may be NA
          beta_diag  <- coef(diag_model)[trm]                 # always there
          if (is.na(beta_const)) beta_const <- beta_diag      # centre at 0
          
          idx      <- grep(sprintf("^tt\\(%s\\)", trm), names(coef(multi)))
          beta_tv  <- coef(multi)[idx]
          if (!length(beta_tv) || all(is.na(beta_tv))) beta_tv <- 0
          
          ## 4c – g(t) matrix (n × m)
          g_t <- timeDepFun(1, x_grid)        # vector or matrix
          g_t <- as.matrix(g_t)               # force matrix
          m    <- ncol(g_t)                   # number of spline columns
          
          ## 4d – beta_tv, padded to match m
          idx      <- grep(sprintf("^tt\\(%s\\)", trm), names(coef(multi)))
          beta_tv  <- as.numeric(coef(multi)[idx])      # numeric, may be length 0/1/m
          
          if (length(beta_tv) == 0L)            beta_tv <- rep(0, m)
          if (length(beta_tv) == 1L && m > 1L)  beta_tv <- rep(beta_tv, m)
          
          ## final dimension check (defensive)
          stopifnot(length(beta_tv) == m)
          
          ## 4e – model curve
          model_curve <- as.vector(g_t %*% beta_tv) + (beta_const - beta_diag)
          
          ## 4f – plot overlay if informative
          if (any(is.finite(model_curve) & model_curve != 0, na.rm = TRUE)) {
            gp_base <- gp_base +
              ggplot2::geom_line(
                data = data.frame(time = x_grid, mc = model_curve),
                ggplot2::aes(time, mc),
                colour = "red", linetype = 2, linewidth = .8
              )
          } else {
            gp_base <- gp_base +
              ggplot2::geom_hline(yintercept = 0,
                                  colour = "red", linetype = 2, linewidth = .6)
          }
          
          print(gp_base)
        }
        
        dev.off()

      } else {
        # Standard PH diagnostics when no tt() terms are used
        cox.z <- survival::cox.zph(multi)
        terms <- setdiff(rownames(cox.z$table), "GLOBAL")   # exclude GLOBAL row
        plot_list <- lapply(terms, function(trm) survminer::ggcoxzph(cox.z, var = trm))
        
        fname <- sprintf("PHdiag.%s.%s.pdf",
                         make.names(outcomeVar), make.names(pred))
        grDevices::cairo_pdf(file.path(outDir, fname), onefile = TRUE)
        for (p in plot_list) print(p)
        dev.off()
        
        
      }
    }
    
    
  }
  
  res      <- dplyr::bind_rows(res_list)
  metrics  <- dplyr::bind_rows(metrics_list)
  
  invisible(list(results = res,
                 followUp = mf,
                 metrics  = metrics))
  
}

# END ==========================================================================
