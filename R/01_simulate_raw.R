# ---------------------------------------------------------------------------
# 01_simulate_raw.R - generates SYNTHETIC SDTM-like raw data for a Phase I
# first-in-human dose-escalation study in advanced solid tumours.
#
# Built-in "real world" features the downstream pipeline must handle:
#   * BLQ concentrations stored as character ("<50") and non-linear (TMDD-like) PK
#   * Two sites reporting haemoglobin in different units (g/L vs g/dL)
#   * Partial AE start dates (YYYY-MM)
#   * Not-done samples / missing results
#   * Semi-structured (nested JSON) biomarker data from two assays
#   * A mid-study protocol amendment adding an exploratory ctDNA endpoint:
#     patients already on study have no baseline ctDNA sample
#   * Staggered cohort enrolment and a data cut-off with ongoing subjects
# ---------------------------------------------------------------------------

simulate_raw <- function(cfg, out_dir = "data/raw") {
  set.seed(cfg$seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cut <- as.Date(cfg$data_cut_date); amend <- as.Date(cfg$amendment1_date)
  cyc <- cfg$cycle_length_days; doses <- cfg$doses_mg
  n <- sum(cfg$cohort_sizes)
  fmt_d  <- function(d) format(d, "%Y-%m-%d")
  fmt_dt <- function(t) format(t, "%Y-%m-%dT%H:%M")

  # ---- subject level ------------------------------------------------------
  s <- data.frame(COHORT = rep(seq_along(doses), cfg$cohort_sizes))
  s$DOSE    <- doses[s$COHORT]
  s$STUDYID <- cfg$study_id
  s$SUBJID  <- sprintf("%03d", seq_len(n))
  s$USUBJID <- paste(s$STUDYID, s$SUBJID, sep = "-")
  s$SITEID  <- sample(c("01", "02"), n, TRUE, prob = c(.6, .4))
  s$AGE     <- as.integer(pmin(pmax(round(rnorm(n, 60, 11)), 28), 82))
  s$SEX     <- sample(c("M", "F"), n, TRUE)
  s$RACE    <- sample(c("WHITE", "ASIAN", "BLACK OR AFRICAN AMERICAN", "OTHER"), n, TRUE,
                      prob = c(.7, .15, .1, .05))
  s$ECOG    <- sample(0:1, n, TRUE, prob = c(.45, .55))
  s$TUMTYPE <- sample(c("NSCLC", "COLORECTAL", "RENAL CELL", "MELANOMA"), n, TRUE)
  idx <- ave(seq_len(n), s$COHORT, FUN = seq_along)
  s$FDOSE   <- as.Date(cfg$first_dose_date) + cfg$cohort_stagger_days * (s$COHORT - 1) + 2 * (idx - 1)
  s$ARM     <- sprintf("DRUGX %d mg Q3W", s$DOSE)
  s$ARMCD   <- paste0("C", s$COHORT)

  # ---- latent individual parameters ---------------------------------------
  B   <- pmin(pmax(rnorm(n, 68, 22), 20), 160)                        # baseline SLD (mm)
  kd  <- 0.012 * s$DOSE / (s$DOSE + 200) * rlnorm(n, 0, 0.6)           # drug-induced shrinkage
  kg  <- 0.0022 * rlnorm(n, 0, 0.7)                                     # tumour growth
  sld_true <- function(i, t) B[i] * (exp(-kd[i] * t) + exp(kg[i] * t) - 1)
  V   <- 3.5 * rlnorm(n, 0, 0.2)                                       # volume (L)
  CL  <- 0.2 * rlnorm(n, 0, 0.3) * (1 + 6 / (1 + s$DOSE / 30))         # non-linear CL (L/day)
  kel <- CL / V

  ad <- cfg$tumor$assessment_days
  tr_rows <- list(); ex_rows <- list(); ds_rows <- list()
  eot_day <- numeric(n); doses_off <- vector("list", n); vis_keep <- vector("list", n)
  max_fu <- as.integer(cut - s$FDOSE)

  for (i in seq_len(n)) {
    base <- round(B[i] * exp(rnorm(1, 0, 0.03)), 1)
    nadir <- base; newles <- FALSE; pd_day <- Inf; vis <- list()
    for (j in seq_along(ad)) {
      v <- round(sld_true(i, ad[j]) * exp(rnorm(1, 0, 0.04)), 1)
      if (!newles && v > 1.1 * nadir && runif(1) < 0.3) newles <- TRUE
      vis[[j]] <- list(j = j, day = ad[j] + sample(-3:3, 1), val = v, nl = newles)
      pd <- newles || (v >= (1 + cfg$tumor$pd_threshold_pct / 100) * nadir &&
                         (v - nadir) >= cfg$tumor$pd_min_abs_mm)
      nadir <- min(nadir, v)
      if (pd) { pd_day <- vis[[j]]$day; break }
    }
    p_ae <- 0.05 + 0.15 * s$DOSE[i] / 900
    ae_day <- if (runif(1) < p_ae) round(runif(1, 22, 150)) else Inf
    wd_day <- if (runif(1) < 0.10) round(runif(1, 22, 200)) else Inf
    ev <- c(PD = pd_day, AE = ae_day, WD = wd_day)
    last_ad <- max(ad) + 7
    ev_min <- min(ev)
    if (is.finite(ev_min) && ev_min <= max_fu[i]) {
      eot <- ev_min
      rsn <- c(PD = "PROGRESSIVE DISEASE", AE = "ADVERSE EVENT",
               WD = "WITHDRAWAL BY SUBJECT")[[names(which.min(ev))]]
    } else if (!is.finite(ev_min) && last_ad <= max_fu[i]) {
      eot <- last_ad; rsn <- "STUDY TREATMENT COMPLETED"
    } else {
      eot <- max_fu[i]; rsn <- "ONGOING"
    }
    eot_day[i] <- eot

    # dosing: Q3W, occasional 0-2 day delays, only before end of treatment / cut-off
    off <- cyc * (0:19) + c(0, sample(0:2, 19, TRUE, prob = c(.7, .2, .1)))
    off <- off[off < eot & off <= max_fu[i]]
    doses_off[[i]] <- off
    ex_rows[[i]] <- data.frame(
      STUDYID = s$STUDYID[i], USUBJID = s$USUBJID[i], EXSEQ = seq_along(off),
      EXTRT = "DRUGX", EXDOSE = s$DOSE[i], EXDOSU = "mg", EXROUTE = "INTRAVENOUS",
      EXSTDTC = fmt_dt(as.POSIXct(paste(s$FDOSE[i] + off, "09:00"), tz = "UTC")),
      EXENDTC = fmt_dt(as.POSIXct(paste(s$FDOSE[i] + off, "10:00"), tz = "UTC")),
      stringsAsFactors = FALSE)
    ds_rows[[i]] <- data.frame(STUDYID = s$STUDYID[i], USUBJID = s$USUBJID[i], DSSEQ = 1,
                               DSDECOD = rsn, DSSTDTC = fmt_d(s$FDOSE[i] + eot),
                               stringsAsFactors = FALSE)
    keep <- Filter(function(z) z$day <= eot && z$day <= max_fu[i], vis)
    vis_keep[[i]] <- keep
    rows <- data.frame(VISIT = "BASELINE", TRDTC = fmt_d(s$FDOSE[i] - sample(5:14, 1)),
                       TRORRES = base, NEWLES = "N", stringsAsFactors = FALSE)
    for (z in keep)
      rows <- rbind(rows, data.frame(VISIT = paste0("WK", 6 * z$j), TRDTC = fmt_d(s$FDOSE[i] + z$day),
                                     TRORRES = z$val, NEWLES = ifelse(z$nl, "Y", "N"),
                                     stringsAsFactors = FALSE))
    tr_rows[[i]] <- cbind(STUDYID = s$STUDYID[i], USUBJID = s$USUBJID[i], TRSEQ = seq_len(nrow(rows)),
                          TRTESTCD = "SUMDIAM", TRTEST = "Sum of Diameters of Target Lesions",
                          rows, TRORRESU = "mm", stringsAsFactors = FALSE)
  }
  ex <- do.call(rbind, ex_rows); ds <- do.call(rbind, ds_rows); tr <- do.call(rbind, tr_rows)

  # ---- DM -----------------------------------------------------------------
  dm <- data.frame(STUDYID = s$STUDYID, USUBJID = s$USUBJID, SUBJID = s$SUBJID, SITEID = s$SITEID,
                   AGE = s$AGE, AGEU = "YEARS", SEX = s$SEX, RACE = s$RACE, ECOG = s$ECOG,
                   TUMTYPE = s$TUMTYPE, ARMCD = s$ARMCD, ARM = s$ARM, RFSTDTC = fmt_d(s$FDOSE),
                   stringsAsFactors = FALSE)

  # ---- PC (cycle 1 PK profile; 504 h sample = pre-dose cycle 2) ------------
  hrs <- cfg$pk$hours; lloq <- cfg$pk$lloq_ng_ml
  pc_rows <- list(); k <- 0
  for (i in seq_len(n)) {
    t0 <- as.POSIXct(paste(s$FDOSE[i], "09:00"), tz = "UTC")
    has2 <- length(doses_off[[i]]) >= 2
    for (h in hrs) {
      if (h == 504 && !has2) next
      mins <- if (h == 0) -30 + round(rnorm(1, 0, 3)) else if (h <= 4) h * 60 + round(rnorm(1, 0, 2))
              else h * 60 + round(runif(1, -60, 60))
      smp <- t0 + mins * 60
      if (as.Date(smp) > cut) next
      conc <- if (h == 0) 0 else (s$DOSE[i] * 1000 / V[i]) * exp(-kel[i] / 24 * max(h - 1, 0)) *
                exp(rnorm(1, 0, 0.15))
      stat <- ""
      res <- if (h == 0 || conc < lloq) paste0("<", lloq) else format(signif(conc, 4), scientific = FALSE, trim = TRUE)
      if (h > 0 && runif(1) < 0.02) { res <- ""; stat <- "NOT DONE" }
      k <- k + 1
      pc_rows[[k]] <- data.frame(STUDYID = s$STUDYID[i], USUBJID = s$USUBJID[i], PCSEQ = 0,
        PCTESTCD = "DRUGX", PCTEST = "DRUGX serum concentration", PCORRES = res, PCORRESU = "ng/mL",
        PCLLOQ = lloq, PCSTAT = stat, PCTPT = paste0(h, "H"), PCDTC = fmt_dt(smp),
        stringsAsFactors = FALSE)
    }
  }
  pc <- do.call(rbind, pc_rows)
  pc$PCSEQ <- ave(seq_len(nrow(pc)), pc$USUBJID, FUN = seq_along)

  # ---- AE -------------------------------------------------------------------
  terms <- data.frame(
    PT = c("FATIGUE", "NAUSEA", "DIARRHOEA", "RASH MACULO-PAPULAR", "HEADACHE", "PYREXIA",
           "ALANINE AMINOTRANSFERASE INCREASED", "ANAEMIA", "NEUTROPHIL COUNT DECREASED",
           "CYTOKINE RELEASE SYNDROME", "INFUSION RELATED REACTION", "DECREASED APPETITE",
           "VOMITING", "ARTHRALGIA"),
    verb = c("Fatigue", "Nausea", "Diarrhea", "Rash", "Headache", "Fever", "ALT increased", "Anemia",
             "Neutropenia", "Cytokine release syndrome", "Infusion reaction", "Poor appetite",
             "Vomiting", "Joint pain"),
    w = c(10, 8, 6, 5, 5, 5, 4, 4, 3, 2, 3, 3, 4, 2),
    dosew = c(0, 0, 0, 1, 0, 1, 2, 1, 2, 6, 4, 0, 0, 0), stringsAsFactors = FALSE)
  ae_rows <- list(); k <- 0
  for (i in seq_len(n)) {
    lam <- 1 + 0.15 * length(doses_off[[i]]) + 1.5 * s$DOSE[i] / 900
    na <- rpois(1, lam); dosef <- s$DOSE[i] / 900
    for (a in seq_len(na)) {
      tm <- terms[sample(nrow(terms), 1, prob = terms$w + terms$dosew * s$DOSE[i] / 300), ]
      infus <- tm$PT %in% c("CYTOKINE RELEASE SYNDROME", "INFUSION RELATED REACTION")
      st_off <- if (infus && runif(1) < 0.8) sample(0:1, 1) else floor(runif(1) * (eot_day[i] + 14))
      sdate <- s$FDOSE[i] + st_off
      if (sdate > cut) next
      gr <- sample(1:4, 1, prob = c(.5, .3, .15 * dosef + .03, .02 * dosef + .005))
      edate <- sdate + sample(1:14, 1)
      partial <- runif(1) < 0.06
      k <- k + 1
      ae_rows[[k]] <- data.frame(STUDYID = s$STUDYID[i], USUBJID = s$USUBJID[i], AESEQ = 0,
        AETERM = ifelse(runif(1) < .5, tm$verb, tolower(tm$verb)), AEDECOD = tm$PT,
        AESTDTC = if (partial) format(sdate, "%Y-%m") else fmt_d(sdate),
        AEENDTC = if (edate > cut) "" else fmt_d(edate), AETOXGR = gr,
        AESER = ifelse(gr >= 3 && runif(1) < .4, "Y", "N"),
        AEREL = ifelse(infus || runif(1) < .5, "RELATED", "NOT RELATED"),
        AEOUT = ifelse(edate > cut, "NOT RECOVERED/NOT RESOLVED", "RECOVERED/RESOLVED"),
        sdate = sdate, stringsAsFactors = FALSE)
    }
  }
  ae <- do.call(rbind, ae_rows)
  ae <- ae[order(ae$USUBJID, ae$sdate), ]
  ae$AESEQ <- ave(seq_len(nrow(ae)), ae$USUBJID, FUN = seq_along)
  ae$sdate <- NULL

  # ---- LB (two sites; site 02 reports haemoglobin in g/dL) -------------------
  labs <- data.frame(cd = c("ALT", "AST", "HGB", "ANC", "CREAT"),
    nm = c("Alanine Aminotransferase", "Aspartate Aminotransferase", "Hemoglobin",
           "Absolute Neutrophil Count", "Creatinine"),
    u = c("U/L", "U/L", "g/L", "10^9/L", "umol/L"),
    lo = c(7, 10, 120, 1.5, 45), hi = c(40, 40, 160, 7.5, 110), stringsAsFactors = FALSE)
  lb_rows <- list(); k <- 0
  for (i in seq_len(n)) {
    nd <- length(doses_off[[i]]); dosef <- s$DOSE[i] / 900
    offs <- sort(unique(c(-7, 7, 14, cyc * seq_len(max(nd - 1, 0)))))
    offs <- offs[offs <= eot_day[i] & offs <= max_fu[i]]
    vl <- ifelse(offs == -7, "SCREENING", ifelse(offs == 7, "C1D8", ifelse(offs == 14, "C1D15",
          paste0("C", offs / cyc + 1, "D1"))))
    vn <- seq_along(offs) - 1
    rsn <- ds$DSDECOD[i]
    if (rsn != "ONGOING" && !(eot_day[i] %in% offs)) {
      offs <- c(offs, eot_day[i]); vl <- c(vl, "EOT"); vn <- c(vn, 99)
    }
    b <- c(ALT = 25 * rlnorm(1, 0, .3), AST = 24 * rlnorm(1, 0, .3), HGB = rnorm(1, 135, 12),
           ANC = 4.5 * rlnorm(1, 0, .3), CREAT = rnorm(1, 78, 12))
    for (v in seq_along(offs)) {
      post <- offs[v] > 0
      val <- c(
        ALT = b[["ALT"]] * rlnorm(1, 0, .15) * (1 + 2.5 * dosef * post * runif(1, .3, 1)),
        AST = b[["AST"]] * rlnorm(1, 0, .15) * (1 + 2.0 * dosef * post * runif(1, .3, 1)),
        HGB = b[["HGB"]] - 0.25 * max(offs[v], 0) / 7 - 6 * dosef * post + rnorm(1, 0, 4),
        ANC = b[["ANC"]] * (1 - 0.35 * dosef * post) * rlnorm(1, 0, .25),
        CREAT = b[["CREAT"]] + rnorm(1, 0, 5))
      for (l in seq_len(nrow(labs))) {
        g_dl <- labs$cd[l] == "HGB" && s$SITEID[i] == "02"
        x <- val[[labs$cd[l]]]; lo <- labs$lo[l]; hi <- labs$hi[l]; u <- labs$u[l]
        if (g_dl) { x <- x / 10; lo <- lo / 10; hi <- hi / 10; u <- "g/dL" }
        dec <- c(ALT = 0, AST = 0, HGB = ifelse(g_dl, 1, 0), ANC = 2, CREAT = 0)[[labs$cd[l]]]
        nd_flag <- runif(1) < 0.01
        k <- k + 1
        lb_rows[[k]] <- data.frame(STUDYID = s$STUDYID[i], USUBJID = s$USUBJID[i], LBSEQ = 0,
          LBTESTCD = labs$cd[l], LBTEST = labs$nm[l], LBORRES = ifelse(nd_flag, NA, round(x, dec)),
          LBORRESU = u, LBSTNRLO = lo, LBSTNRHI = hi, LBSTAT = ifelse(nd_flag, "NOT DONE", ""),
          VISIT = vl[v], VISITNUM = vn[v], LBDTC = fmt_d(s$FDOSE[i] + offs[v]), stringsAsFactors = FALSE)
      }
    }
  }
  lb <- do.call(rbind, lb_rows)
  lb$LBSEQ <- ave(seq_len(nrow(lb)), lb$USUBJID, FUN = seq_along)

  # ---- BM: semi-structured JSON lines (two assays, nested records) ------------
  bm_lines <- character(); sid <- 0
  cxh <- c(0, 24, 72, 168, 504); cxlab <- c("C1D1 PRE", "C1D2 24H", "C1D4 72H", "C1D8 168H", "C2D1 PRE")
  toj <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null",
                                                  digits = NA))
  for (i in seq_len(n)) {
    t0 <- as.POSIXct(paste(s$FDOSE[i], "09:00"), tz = "UTC")
    has2 <- length(doses_off[[i]]) >= 2
    bcx <- 90 * rlnorm(1, 0, .35)
    for (q in seq_along(cxh)) {
      if (cxh[q] == 504 && !has2) next
      smp <- if (cxh[q] == 0) t0 - (30 + round(rnorm(1, 0, 3))) * 60 else
        t0 + (cxh[q] * 60 + round(runif(1, -60, 60))) * 60
      if (as.Date(smp) > cut) next
      fold <- 1 + 5 * s$DOSE[i] / (s$DOSE[i] + 150) * (1 - exp(-cxh[q] / 24)) * exp(-cxh[q] / 200)
      val <- bcx * fold * rlnorm(1, 0, .2)
      sid <- sid + 1
      blq <- val < cfg$biomarker$cxcl10_lloq_pg_ml
      rec <- list(usubjid = s$USUBJID[i],
        sample = list(id = sprintf("BM%05d", sid), matrix = "serum", collected = fmt_dt(smp),
                      nominal_visit = cxlab[q], nominal_hr = cxh[q]),
        assay = list(name = "CXCL10 Luminex", analyte = "CXCL10", unit = "pg/mL",
                     lloq = cfg$biomarker$cxcl10_lloq_pg_ml, run_id = sprintf("R%03d", sample(1:12, 1))),
        result = if (blq) list(value = NULL, flag = "BLQ") else list(value = round(val, 1), flag = NULL))
      bm_lines <- c(bm_lines, toj(rec))
    }
    # exploratory ctDNA (added by amendment 1): only samples collected on/after amendment date
    b_vaf <- if (runif(1) < 0.25) 0 else 1.2 * rlnorm(1, 0, .8)
    ct_days <- c(0, vapply(vis_keep[[i]], function(z) as.numeric(z$day), 0))
    ct_lab  <- c("C1D1 PRE", vapply(vis_keep[[i]], function(z) paste0("WK", 6 * z$j), ""))
    for (q in seq_along(ct_days)) {
      d <- s$FDOSE[i] + ct_days[q]
      if (d < amend || d > cut) next
      smp <- as.POSIXct(paste(d, ifelse(ct_days[q] == 0, "08:30", "10:15")), tz = "UTC")
      vaf <- b_vaf * (sld_true(i, ct_days[q]) / B[i])^1.5 * rlnorm(1, 0, .25)
      sid <- sid + 1
      blq <- vaf < cfg$biomarker$ctdna_lloq_pct
      vars <- list()
      if (!blq) {
        ng <- sample(1:3, 1); genes <- sample(c("KRAS", "TP53", "EGFR", "BRAF", "PIK3CA", "APC"), ng)
        vv <- c(vaf, vaf * runif(ng - 1, .2, .9))[sample(ng)]
        vars <- lapply(seq_len(ng), function(g) list(gene = genes[g], vaf = round(vv[g], 3)))
      }
      rec <- list(usubjid = s$USUBJID[i],
        sample = list(id = sprintf("BM%05d", sid), matrix = "plasma", collected = fmt_dt(smp),
                      nominal_visit = ct_lab[q], nominal_hr = NULL),
        assay = list(name = "ctDNA NGS panel v1", analyte = "ctDNA", unit = "%",
                     lloq = cfg$biomarker$ctdna_lloq_pct, run_id = sprintf("N%03d", sample(1:6, 1))),
        result = list(flag = if (blq) "BLQ" else NULL, variants = vars))
      bm_lines <- c(bm_lines, toj(rec))
    }
  }

  # ---- write -------------------------------------------------------------------
  wr <- function(x, f) write.csv(x, file.path(out_dir, f), row.names = FALSE, na = "")
  wr(dm, "dm.csv"); wr(ex, "ex.csv"); wr(ds, "ds.csv"); wr(pc, "pc.csv"); wr(ae, "ae.csv")
  wr(lb, "lb.csv"); wr(tr, "tr.csv"); writeLines(bm_lines, file.path(out_dir, "bm.jsonl"))
  log_msg(sprintf("Synthetic raw data written to %s (%d subjects)", out_dir, n))
  invisible(TRUE)
}
