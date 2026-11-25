# ---------------------------------------------------------------------------
# 04_adbm_adtr_adrs.R - exploratory biomarkers (semi-structured JSON input),
# tumour measurements and best overall response
# ---------------------------------------------------------------------------

# Flatten nested JSON-lines biomarker data into one row per sample/analyte.
read_bm <- function(path) {
  recs <- lapply(readLines(path), jsonlite::fromJSON, simplifyVector = FALSE)
  rows <- lapply(recs, function(r) {
    vars <- r$result$variants
    tibble::tibble(
      USUBJID = r$usubjid, SAMPID = r$sample$id, MATRIX = r$sample$matrix,
      BMDTC = r$sample$collected, BMTPT = r$sample$nominal_visit,
      ANAME = r$assay$name, ANALYTE = r$assay$analyte, BMORRESU = r$assay$unit,
      BMLLOQ = as.numeric(r$assay$lloq), RUNID = r$assay$run_id,
      VALUE = if (is.null(r$result$value)) NA_real_ else as.numeric(r$result$value),
      FLAG  = if (is.null(r$result$flag)) NA_character_ else r$result$flag,
      NVAR  = length(vars),
      MAXVAF = if (length(vars)) max(vapply(vars, function(v) as.numeric(v$vaf), 0)) else NA_real_)
  })
  bind_rows(rows)
}

build_adbm <- function(raw, adsl, cfg) {
  read_bm(raw$bm_path) %>%
    left_join(adsl %>% select(USUBJID, COHORT, TRT01PN, TRTSDTM), by = "USUBJID") %>%
    mutate(
      ADTM = to_dtm(BMDTC),
      ARELTM = round(as.numeric(difftime(ADTM, TRTSDTM, units = "hours")), 2),
      PARAMCD = ifelse(ANALYTE == "CXCL10", "CXCL10", "CTDNAVAF"),
      PARAM = ifelse(ANALYTE == "CXCL10", "CXCL10 (pg/mL)", "ctDNA max variant allele fraction (%)"),
      AVALU = BMORRESU, LLOQ = BMLLOQ, AVISIT = BMTPT,
      BLQFL = ifelse(!is.na(FLAG) & FLAG == "BLQ", "Y", "N"),
      # BLQ -> LLOQ/2 (SAP 8.2); ctDNA value = max VAF over nested variant calls
      AVAL  = ifelse(BLQFL == "Y", LLOQ / 2, coalesce(VALUE, MAXVAF)),
      DTYPE = ifelse(BLQFL == "Y", "LLOQ/2", NA_character_)) %>%
    group_by(USUBJID, PARAMCD) %>% arrange(ADTM, .by_group = TRUE) %>%
    mutate(
      .bl   = flag_last(ADTM, ADTM <= TRTSDTM & !is.na(AVAL)),
      ABLFL = ifelse(.bl, "Y", NA_character_),
      BASE  = if (any(.bl)) AVAL[.bl][1] else NA_real_,
      .bblq = if (any(.bl)) BLQFL[.bl][1] == "Y" else NA,
      BASEAVL = ifelse(!is.na(BASE), "Y", "N"),
      post  = is.na(ABLFL) & ADTM > TRTSDTM,
      ok    = post & !is.na(BASE) & !coalesce(.bblq, TRUE),   # no ratios off a BLQ baseline (SAP 8.3)
      CHG   = ifelse(ok, AVAL - BASE, NA_real_),
      PCHG  = ifelse(ok, 100 * CHG / BASE, NA_real_),
      R2BASE = ifelse(ok, AVAL / BASE, NA_real_)) %>%
    ungroup() %>% arrange(USUBJID, PARAMCD, ADTM) %>%
    mutate(STUDYID = cfg$study_id) %>%
    select(STUDYID, USUBJID, COHORT, TRT01PN, SAMPID, PARAMCD, PARAM, AVISIT, ADTM, ARELTM,
           AVAL, AVALU, LLOQ, BLQFL, DTYPE, BASE, CHG, PCHG, R2BASE, ABLFL, BASEAVL,
           RUNID, ANAME, MATRIX, NVAR)
}

build_adtr <- function(raw, adsl, cfg) {
  tp <- cfg$tumor
  raw$tr %>%
    left_join(adsl %>% select(USUBJID, COHORT, TRT01PN, TRTSDT), by = "USUBJID") %>%
    mutate(
      ADT = to_date(TRDTC), ADY = study_day(ADT, TRTSDT),
      PARAMCD = "SUMDIAM", PARAM = "Sum of diameters of target lesions (mm)",
      AVISIT = VISIT, AVISITN = ifelse(VISIT == "BASELINE", 0L, suppressWarnings(as.integer(sub("WK", "", VISIT)))),
      AVAL = TRORRES, AVALU = TRORRESU, NEWLESFL = NEWLES) %>%
    group_by(USUBJID) %>% arrange(ADT, .by_group = TRUE) %>%
    mutate(
      .bl   = flag_last(ADT, ADT <= TRTSDT & !is.na(AVAL)),
      ABLFL = ifelse(.bl, "Y", NA_character_),
      BASE  = if (any(.bl)) AVAL[.bl][1] else NA_real_,
      NADIR = dplyr::lag(cummin(AVAL)),            # smallest earlier value incl. baseline
      post  = is.na(ABLFL) & ADT > TRTSDT,
      CHG   = ifelse(post, AVAL - BASE, NA_real_),
      PCHG  = ifelse(post, 100 * CHG / BASE, NA_real_),
      ABSINC = ifelse(post, AVAL - NADIR, NA_real_),
      PCHGNAD = ifelse(post, 100 * (AVAL - NADIR) / NADIR, NA_real_),
      # simplified RECIST 1.1 target-lesion response (SAP 9.1); no CR (non-target not collected)
      AVALC = case_when(
        !post ~ NA_character_,
        NEWLESFL == "Y" | (PCHGNAD >= tp$pd_threshold_pct - 1e-9 & ABSINC >= tp$pd_min_abs_mm - 1e-9) ~ "PD",
        PCHG <= tp$pr_threshold_pct + 1e-9 ~ "PR",
        TRUE ~ "SD")) %>%
    ungroup() %>% arrange(USUBJID, TRSEQ) %>%
    select(STUDYID, USUBJID, COHORT, TRT01PN, TRSEQ, PARAMCD, PARAM, AVISIT, AVISITN, ADT, ADY,
           AVAL, AVALU, BASE, CHG, PCHG, NADIR, PCHGNAD, ABSINC, NEWLESFL, AVALC, ABLFL)
}

build_adrs <- function(adtr, adsl) {
  rk <- c(PD = 1, SD = 2, PR = 3)
  bor <- adtr %>% filter(!is.na(AVALC)) %>% mutate(rank = rk[AVALC]) %>%
    group_by(USUBJID) %>% arrange(ADT, .by_group = TRUE) %>%
    summarise(NASSESS = n(), BESTPCHG = min(PCHG), AVALC = c("PD", "SD", "PR")[max(rank)],
              BORDT = ADT[which(rank == max(rank))[1]], .groups = "drop")
  adsl %>% select(STUDYID, USUBJID, COHORT, TRT01PN, TRTSDT) %>%
    left_join(bor, by = "USUBJID") %>%
    mutate(PARAMCD = "BOR", PARAM = "Best overall response (unconfirmed)",
           AVALC = coalesce(AVALC, "NE"), NASSESS = coalesce(NASSESS, 0L),
           BORADY = study_day(BORDT, TRTSDT),
           ORRFL = ifelse(AVALC == "PR", "Y", "N"),
           DCRFL = ifelse(AVALC %in% c("PR", "SD"), "Y", "N")) %>%
    select(STUDYID, USUBJID, COHORT, TRT01PN, PARAMCD, PARAM, AVALC, NASSESS, BESTPCHG,
           BORDT, BORADY, ORRFL, DCRFL)
}
