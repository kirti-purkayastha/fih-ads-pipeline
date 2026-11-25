# ---------------------------------------------------------------------------
# 03_adlb_adpc_adpp.R - laboratory, PK concentration and PK parameter datasets
# ---------------------------------------------------------------------------

# ---- ADLB -------------------------------------------------------------------
build_adlb <- function(raw, adsl, cfg) {
  raw$lb %>%
    left_join(adsl %>% select(USUBJID, COHORT, TRT01PN, TRTSDT), by = "USUBJID") %>%
    mutate(
      ADT = to_date(LBDTC), ADY = study_day(ADT, TRTSDT),
      # unit harmonisation (SAP 7.1): haemoglobin g/dL -> g/L
      fac   = ifelse(LBTESTCD == "HGB" & LBORRESU == "g/dL", 10, 1),
      AVAL  = LBORRES * fac, ANRLO = LBSTNRLO * fac, ANRHI = LBSTNRHI * fac,
      AVALU = ifelse(LBTESTCD == "HGB", "g/L", LBORRESU),
      PARAMCD = LBTESTCD, PARAM = paste0(LBTEST, " (", AVALU, ")"),
      AVISIT = VISIT, AVISITN = as.integer(VISITNUM),
      ANRIND = case_when(is.na(AVAL) ~ NA_character_, AVAL < ANRLO ~ "LOW",
                         AVAL > ANRHI ~ "HIGH", TRUE ~ "NORMAL")) %>%
    group_by(USUBJID, PARAMCD) %>% arrange(ADT, AVISITN, .by_group = TRUE) %>%
    mutate(
      .bl   = flag_last(ADT, ADT <= TRTSDT & !is.na(AVAL)),
      ABLFL = ifelse(.bl, "Y", NA_character_),
      BASE  = if (any(.bl)) AVAL[.bl][1] else NA_real_,
      BNRIND = if (any(.bl)) ANRIND[.bl][1] else NA_character_,
      post  = is.na(ABLFL) & ADT > TRTSDT,
      CHG   = ifelse(post & !is.na(BASE), AVAL - BASE, NA_real_),
      PCHG  = ifelse(post & !is.na(BASE) & BASE != 0, 100 * CHG / BASE, NA_real_),
      SHIFT1 = ifelse(post & !is.na(BNRIND) & !is.na(ANRIND), paste(BNRIND, "to", ANRIND), NA_character_)) %>%
    ungroup() %>% arrange(USUBJID, LBSEQ) %>%
    select(STUDYID, USUBJID, COHORT, TRT01PN, LBSEQ, PARAMCD, PARAM, AVISIT, AVISITN, ADT, ADY,
           AVAL, AVALU, ANRLO, ANRHI, ANRIND, BASE, BNRIND, CHG, PCHG, SHIFT1, ABLFL,
           LBORRES, LBORRESU)
}

# ---- ADPC -------------------------------------------------------------------
build_adpc <- function(raw, adsl, cfg) {
  raw$pc %>%
    left_join(adsl %>% select(USUBJID, COHORT, TRT01PN, TRTSDTM), by = "USUBJID") %>%
    mutate(
      ADTM   = to_dtm(PCDTC),
      ARELTM = round(as.numeric(difftime(ADTM, TRTSDTM, units = "hours")), 2),
      ATPTN  = as.numeric(sub("H$", "", PCTPT)),
      ATPT   = ifelse(ATPTN == 0, "PRE-DOSE", ifelse(ATPTN == 1, "1H (END OF INFUSION)",
                 ifelse(ATPTN == 504, "504H (C2D1 PRE-DOSE)", PCTPT))),
      notdone = PCSTAT == "NOT DONE",
      BLQFL  = ifelse(notdone, NA_character_, ifelse(startsWith(PCORRES, "<"), "Y", "N")),
      # BLQ -> 0 (SAP 5.1); not-done samples stay missing
      AVAL   = ifelse(notdone, NA_real_,
                      ifelse(BLQFL == "Y", 0, suppressWarnings(as.numeric(PCORRES)))),
      PARAMCD = "DRUGX", PARAM = "DRUGX serum concentration (ng/mL)", AVALU = PCORRESU,
      AVISIT = "CYCLE 1", LLOQ = PCLLOQ,
      ANL01FL = ifelse(!is.na(AVAL), "Y", NA_character_)) %>%
    select(STUDYID, USUBJID, COHORT, TRT01PN, PCSEQ, PARAMCD, PARAM, AVISIT, ATPT, ATPTN,
           ADTM, ARELTM, AVAL, AVALU, LLOQ, BLQFL, PCORRES, ANL01FL)
}

# ---- NCA --------------------------------------------------------------------
# Non-compartmental analysis of one profile (SAP 5.2):
#  * AUC: linear-up / log-down trapezoidal from time 0 to Tlast (last conc > 0)
#  * lambda_z: log-linear regression through the last `min_pts` quantifiable
#    post-Tmax points; accepted if slope < 0 and adjusted R2 >= `min_r2`
#  * CL (L/h) = Dose(ng) / AUCinf (ng*h/mL) / 1000
nca_profile <- function(time, conc, dose_mg, min_pts = 3, min_r2 = 0.8) {
  ok <- !is.na(time) & !is.na(conc)
  time <- pmax(time[ok], 0); conc <- conc[ok]
  o <- order(time); time <- time[o]; conc <- conc[o]
  out <- c(CMAX = NA, TMAX = NA, AUCLST = NA, AUCIFO = NA, LAMZ = NA, LAMZNPT = NA,
           R2ADJ = NA, THALF = NA, CL = NA, VZ = NA, CMAXDN = NA)
  if (length(conc) < 2 || !any(conc > 0)) return(out)
  out["CMAX"] <- max(conc); out["TMAX"] <- time[which.max(conc)]
  out["CMAXDN"] <- out[["CMAX"]] / dose_mg
  ilast <- max(which(conc > 0))
  auc <- 0
  for (k in seq_len(ilast - 1)) {
    dt <- time[k + 1] - time[k]; c1 <- conc[k]; c2 <- conc[k + 1]
    auc <- auc + if (c2 < c1 && c2 > 0) (c1 - c2) / log(c1 / c2) * dt else (c1 + c2) / 2 * dt
  }
  out["AUCLST"] <- auc
  idx <- which(time > out[["TMAX"]] & conc > 0 & seq_along(conc) <= ilast)
  if (length(idx) >= min_pts) {
    idx <- tail(idx, min_pts)
    fit <- lm(log(conc[idx]) ~ time[idx]); s <- summary(fit)
    r2adj <- 1 - (1 - s$r.squared) * (min_pts - 1) / (min_pts - 2)
    out["R2ADJ"] <- r2adj; out["LAMZNPT"] <- min_pts
    lamz <- -unname(coef(fit)[2])
    if (is.finite(lamz) && lamz > 0 && r2adj >= min_r2) {
      out["LAMZ"] <- lamz; out["THALF"] <- log(2) / lamz
      out["AUCIFO"] <- auc + conc[ilast] / lamz
      out["CL"] <- dose_mg * 1e6 / out[["AUCIFO"]] / 1000
      out["VZ"] <- out[["CL"]] / lamz
    }
  }
  out
}

pp_meta <- tibble::tibble(
  PARAMN  = 1:11,
  PARAMCD = c("CMAX", "TMAX", "AUCLST", "AUCIFO", "LAMZ", "LAMZNPT", "R2ADJ", "THALF", "CL", "VZ", "CMAXDN"),
  PARAM   = c("Maximum concentration", "Time of maximum concentration",
              "AUC from time 0 to last quantifiable concentration",
              "AUC from time 0 extrapolated to infinity", "Terminal elimination rate constant",
              "Number of points used for terminal slope", "Adjusted R-squared of terminal fit",
              "Terminal half-life", "Total clearance", "Terminal volume of distribution",
              "Dose-normalised Cmax"),
  AVALU   = c("ng/mL", "h", "h*ng/mL", "h*ng/mL", "1/h", NA, NA, "h", "L/h", "L", "ng/mL/mg"))

build_adpp <- function(adpc, adsl, cfg) {
  pk_subj <- adsl %>% filter(PKFL == "Y") %>% select(USUBJID, COHORT, TRT01PN)
  res <- lapply(seq_len(nrow(pk_subj)), function(i) {
    u <- pk_subj$USUBJID[i]
    d <- adpc[adpc$USUBJID == u & !is.na(adpc$AVAL), ]
    v <- nca_profile(d$ARELTM, d$AVAL, pk_subj$TRT01PN[i], cfg$pk$lamz_min_points, cfg$pk$lamz_min_r2adj)
    tibble::tibble(USUBJID = u, PARAMCD = names(v), AVAL = unname(v))
  })
  bind_rows(res) %>%
    left_join(pk_subj, by = "USUBJID") %>% left_join(pp_meta, by = "PARAMCD") %>%
    mutate(STUDYID = cfg$study_id, AVISIT = "CYCLE 1") %>%
    arrange(USUBJID, PARAMN) %>%
    select(STUDYID, USUBJID, COHORT, TRT01PN, PARAMN, PARAMCD, PARAM, AVAL, AVALU, AVISIT)
}
