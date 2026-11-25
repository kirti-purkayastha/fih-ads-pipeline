# ---------------------------------------------------------------------------
# 02_adsl_adae.R - ADSL (subject level) and ADAE (adverse events)
# Derivations follow docs/mock_SAP.md; variable-level specs in spec/adam_spec.psv
# ADSL is built in two stages because DLT flags depend on ADAE, which in turn
# needs ADSL treatment dates (common ADaM dependency pattern).
# ---------------------------------------------------------------------------

build_adsl_base <- function(raw, cfg) {
  ex <- raw$ex %>% mutate(EXSTDT = to_date(EXSTDTC), EXSTDTM = to_dtm(EXSTDTC))
  trt <- ex %>% group_by(USUBJID) %>% arrange(EXSTDTM, .by_group = TRUE) %>%
    summarise(TRTSDT = min(EXSTDT), TRTSDTM = min(EXSTDTM), TRTEDT = max(EXSTDT),
              NDOSES = n(), TRT01AN = first(EXDOSE), CUMDOSE = sum(EXDOSE), .groups = "drop")

  # PK-evaluable: >=1 quantifiable post-dose concentration (SAP 3.2)
  pk <- raw$pc %>% left_join(trt %>% select(USUBJID, TRTSDTM), by = "USUBJID") %>%
    mutate(post = to_dtm(PCDTC) > TRTSDTM,
           quant = nzchar(PCORRES) & !startsWith(PCORRES, "<")) %>%
    group_by(USUBJID) %>%
    summarise(PKFL = ifelse(any(post & quant, na.rm = TRUE), "Y", "N"), .groups = "drop")

  dsp <- raw$ds %>% transmute(USUBJID, EOSDT = to_date(DSSTDTC), DCSREAS = DSDECOD,
    EOSSTT = case_when(DSDECOD == "ONGOING" ~ "ONGOING",
                       DSDECOD == "STUDY TREATMENT COMPLETED" ~ "COMPLETED",
                       TRUE ~ "DISCONTINUED"))

  raw$dm %>%
    left_join(trt, by = "USUBJID") %>% left_join(pk, by = "USUBJID") %>% left_join(dsp, by = "USUBJID") %>%
    mutate(
      COHORT  = as.integer(sub("^C", "", ARMCD)),
      TRT01P  = ARM,
      TRT01PN = as.integer(sub("^DRUGX ([0-9]+) mg.*$", "\\1", ARM)),
      TRT01A  = ifelse(is.na(TRT01AN), NA_character_, sprintf("DRUGX %d mg Q3W", as.integer(TRT01AN))),
      AGEGR1  = ifelse(AGE >= 65, ">=65", "<65"),
      TRTDURD = as.integer(TRTEDT - TRTSDT) + 1L,
      SAFFL   = ifelse(!is.na(TRTSDT), "Y", "N"),
      PKFL    = coalesce(PKFL, "N"))
}

build_adae <- function(raw, adsl, cfg) {
  ae <- raw$ae %>%
    left_join(adsl %>% select(USUBJID, COHORT, TRT01PN, TRTSDT, TRTEDT), by = "USUBJID")
  st <- impute_start(ae$AESTDTC, ae$TRTSDT)
  ae %>% mutate(
    ASTDT  = st$date, ASTDTF = st$flag,
    AENDT  = to_date(AEENDTC),
    ASTDY  = study_day(ASTDT, TRTSDT), AENDY = study_day(AENDT, TRTSDT),
    ATOXGRN = as.integer(AETOXGR), ATOXGR = paste("GRADE", AETOXGR),
    AESERFL = AESER, RELFL = ifelse(AEREL == "RELATED", "Y", "N"),
    TRTEMFL = ifelse(!is.na(ASTDT) & ASTDT >= TRTSDT & ASTDT <= TRTEDT + cfg$safety_window_days, "Y", "N"),
    AESIFL  = ifelse(AEDECOD %in% c("CYTOKINE RELEASE SYNDROME", "INFUSION RELATED REACTION"), "Y", "N"),
    DLTFL   = ifelse(TRTEMFL == "Y" & ATOXGRN >= 3 & RELFL == "Y" & ASTDY <= cfg$dlt_window_days, "Y", "N")) %>%
    group_by(USUBJID) %>% arrange(ASTDT, AESEQ, .by_group = TRUE) %>%
    mutate(AOCCFL = ifelse(TRTEMFL == "Y" & cumsum(TRTEMFL == "Y") == 1, "Y", NA_character_)) %>%
    group_by(USUBJID, AEDECOD) %>% arrange(ASTDT, AESEQ, .by_group = TRUE) %>%
    mutate(AOCCPFL = ifelse(TRTEMFL == "Y" & cumsum(TRTEMFL == "Y") == 1, "Y", NA_character_)) %>%
    ungroup() %>% arrange(USUBJID, AESEQ) %>%
    select(STUDYID, USUBJID, AESEQ, COHORT, TRT01PN, AETERM, AEDECOD, AESTDTC, ASTDT, ASTDTF, AENDT,
           ASTDY, AENDY, ATOXGR, ATOXGRN, AESERFL, RELFL, AEOUT, TRTEMFL, AESIFL, DLTFL, AOCCFL, AOCCPFL)
}

adsl_add_ae_flags <- function(adsl, adae, cfg) {
  dlt <- adae %>% group_by(USUBJID) %>%
    summarise(DLTFL = ifelse(any(DLTFL == "Y"), "Y", "N"), .groups = "drop")
  adsl %>% left_join(dlt, by = "USUBJID") %>%
    mutate(DLTFL = coalesce(DLTFL, "N"),
           DLTEVLFL = ifelse(SAFFL == "Y" &
             (DLTFL == "Y" | (as.integer(EOSDT - TRTSDT) + 1L >= cfg$dlt_window_days)), "Y", "N")) %>%
    select(STUDYID, USUBJID, SUBJID, SITEID, COHORT, TRT01P, TRT01PN, TRT01A, TRT01AN, AGE, AGEGR1,
           SEX, RACE, ECOG, TUMTYPE, TRTSDT, TRTSDTM, TRTEDT, NDOSES, CUMDOSE, TRTDURD, EOSDT, EOSSTT,
           DCSREAS, SAFFL, PKFL, DLTFL, DLTEVLFL)
}
