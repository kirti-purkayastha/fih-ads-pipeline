# ---------------------------------------------------------------------------
# 06_checks.R - data integrity & consistency checks on the final ADaM datasets
# Each check returns the NUMBER OF FAILING RECORDS; 0 = PASS.
# Used by the pipeline (hard stop on failure), by testthat, and by CI.
# ---------------------------------------------------------------------------
run_checks <- function(adam, raw, spec, cfg) {
  out <- list()
  chk <- function(id, ds, desc, fails) {
    out[[length(out) + 1]] <<- data.frame(check_id = id, dataset = ds, description = desc,
      n_fail = as.integer(sum(fails, na.rm = TRUE)),
      status = ifelse(sum(fails, na.rm = TRUE) == 0, "PASS", "FAIL"), stringsAsFactors = FALSE)
  }
  cut <- as.Date(cfg$data_cut_date); amend <- as.Date(cfg$amendment1_date)
  adsl <- adam$ADSL

  # 1. Specification conformance (names, order, label, class)
  cls <- c(text = "character", integer = "integer", float = "numeric", date = "Date", datetime = "POSIXct")
  for (ds in names(adam)) {
    sp <- spec[spec$dataset == ds, ]; d <- adam[[ds]]
    bad_names <- !identical(names(d), sp$variable)
    bad_lab <- sum(vapply(seq_len(nrow(sp)), function(i) !identical(attr(d[[sp$variable[i]]], "label"), sp$label[i]), TRUE))
    bad_cls <- sum(vapply(seq_len(nrow(sp)), function(i) !inherits(d[[sp$variable[i]]], cls[[sp$type[i]]]), TRUE))
    chk("SPEC01", ds, "Variables, order, labels and classes match specification", bad_names + bad_lab + bad_cls)
  }

  # 2. Key uniqueness
  keys <- list(ADSL = "USUBJID", ADAE = c("USUBJID", "AESEQ"), ADLB = c("USUBJID", "LBSEQ"),
               ADPC = c("USUBJID", "PCSEQ"), ADPP = c("USUBJID", "PARAMCD"),
               ADBM = c("SAMPID", "PARAMCD"), ADTR = c("USUBJID", "TRSEQ"), ADRS = "USUBJID")
  for (ds in names(keys)) chk("KEY01", ds, paste("No duplicate records by", paste(keys[[ds]], collapse = "+")),
                              sum(duplicated(adam[[ds]][keys[[ds]]])))

  # 3. Referential integrity
  for (ds in setdiff(names(adam), "ADSL"))
    chk("REF01", ds, "All USUBJID exist in ADSL", sum(!adam[[ds]]$USUBJID %in% adsl$USUBJID))

  # 4. Record reconciliation raw -> ADaM
  chk("REC01", "ADSL", "Records = DM records", abs(nrow(adsl) - nrow(raw$dm)))
  chk("REC02", "ADAE", "Records = AE records", abs(nrow(adam$ADAE) - nrow(raw$ae)))
  chk("REC03", "ADLB", "Records = LB records", abs(nrow(adam$ADLB) - nrow(raw$lb)))
  chk("REC04", "ADPC", "Records = PC records", abs(nrow(adam$ADPC) - nrow(raw$pc)))
  chk("REC05", "ADTR", "Records = TR records", abs(nrow(adam$ADTR) - nrow(raw$tr)))
  chk("REC06", "ADBM", "Records = JSON biomarker lines", abs(nrow(adam$ADBM) - length(readLines(raw$bm_path))))
  chk("REC07", "ADSL", "Sum of NDOSES = EX records", abs(sum(adsl$NDOSES) - nrow(raw$ex)))

  # 5. Date logic
  chk("DT01", "ADSL", "TRTSDT <= TRTEDT <= EOSDT <= data cut-off", sum(adsl$TRTEDT < adsl$TRTSDT | adsl$EOSDT < adsl$TRTEDT | adsl$EOSDT > cut))
  chk("DT02", "ADAE", "AENDT >= ASTDT when both present", sum(adam$ADAE$AENDT < adam$ADAE$ASTDT, na.rm = TRUE))
  chk("DT03", "ADAE", "No analysis start date after data cut-off", sum(adam$ADAE$ASTDT > cut, na.rm = TRUE))
  chk("DT04", "ADLB", "No lab date after data cut-off", sum(adam$ADLB$ADT > cut, na.rm = TRUE))
  chk("DT05", "ADBM", "No ctDNA sample before protocol amendment 1", sum(as.Date(adam$ADBM$ADTM) < amend & adam$ADBM$PARAMCD == "CTDNAVAF"))
  chk("DT06", "ADAE", "Imputed start dates (ASTDTF) never precede first dose month", sum(!is.na(adam$ADAE$ASTDTF) & adam$ADAE$ASTDT < adsl$TRTSDT[match(adam$ADAE$USUBJID, adsl$USUBJID)]))

  # 6. Derivation consistency
  for (ds in c("ADLB", "ADBM", "ADTR")) {
    d <- adam[[ds]]; ok <- !is.na(d$CHG)
    chk("DER01", ds, "CHG = AVAL - BASE", sum(abs(d$CHG[ok] - (d$AVAL[ok] - d$BASE[ok])) > 1e-8))
  }
  for (ds in c("ADLB", "ADTR")) {
    d <- adam[[ds]]; g <- if (ds == "ADLB") paste(d$USUBJID, d$PARAMCD) else d$USUBJID
    n_bl <- tapply(d$ABLFL %in% "Y", g, sum)
    chk("DER02", ds, "At most one baseline record per subject/parameter", sum(n_bl > 1))
    chk("DER03", ds, "BASE equals AVAL of the baseline-flagged record", sum(abs(d$BASE[d$ABLFL %in% "Y"] - d$AVAL[d$ABLFL %in% "Y"]) > 1e-8))
  }
  chk("DER04", "ADLB", "Haemoglobin harmonised to g/L", sum(adam$ADLB$PARAMCD == "HGB" & adam$ADLB$AVALU != "g/L"))
  chk("DER05", "ADLB", "Haemoglobin values within physiological range (60-200 g/L)",
      sum(adam$ADLB$PARAMCD == "HGB" & !is.na(adam$ADLB$AVAL) & (adam$ADLB$AVAL < 60 | adam$ADLB$AVAL > 200)))
  chk("DER06", "ADSL", "TRTDURD = TRTEDT - TRTSDT + 1", sum(adsl$TRTDURD != as.integer(adsl$TRTEDT - adsl$TRTSDT) + 1L))
  chk("DER07", "ADPC", "AVAL missing only for not-done samples", sum(is.na(adam$ADPC$AVAL) != (adam$ADPC$PCORRES == "")))
  chk("DER08", "ADPC", "BLQ records have AVAL = 0", sum(adam$ADPC$BLQFL %in% "Y" & adam$ADPC$AVAL != 0, na.rm = TRUE))
  pp <- tidyr::pivot_wider(adam$ADPP[, c("USUBJID", "PARAMCD", "AVAL")], names_from = PARAMCD, values_from = AVAL)
  chk("DER09", "ADPP", "AUCIFO >= AUCLST", sum(pp$AUCIFO < pp$AUCLST, na.rm = TRUE))
  mx <- tapply(adam$ADPC$AVAL, adam$ADPC$USUBJID, max, na.rm = TRUE)
  chk("DER10", "ADPP", "CMAX equals highest observed concentration", sum(abs(pp$CMAX - mx[pp$USUBJID]) > 1e-8, na.rm = TRUE))

  # 7. Controlled terminology / flags
  chk("CT01", "ADLB", "ANRIND in LOW/NORMAL/HIGH", sum(!adam$ADLB$ANRIND %in% c("LOW", "NORMAL", "HIGH", NA)))
  chk("CT02", "ADAE", "ATOXGRN between 1 and 5", sum(!adam$ADAE$ATOXGRN %in% 1:5))
  chk("CT03", "ADRS", "AVALC in PR/SD/PD/NE", sum(!adam$ADRS$AVALC %in% c("PR", "SD", "PD", "NE")))
  chk("CT04", "ADSL", "Flags are Y/N", sum(!unlist(adsl[c("SAFFL", "PKFL", "DLTFL", "DLTEVLFL")]) %in% c("Y", "N")))
  chk("CT05", "ADSL", "SEX in M/F", sum(!adsl$SEX %in% c("M", "F")))
  chk("CT06", "ADTR", "Visit response in PR/SD/PD (post-baseline)", sum(!adam$ADTR$AVALC %in% c("PR", "SD", "PD", NA)))

  # 8. Cross-dataset consistency
  chk("XD01", "ADSL", "SAFFL = Y if and only if NDOSES >= 1", sum((adsl$SAFFL == "Y") != (adsl$NDOSES >= 1)))
  dlt_ae <- unique(adam$ADAE$USUBJID[adam$ADAE$DLTFL == "Y"])
  chk("XD02", "ADSL", "DLTFL consistent with ADAE", sum((adsl$DLTFL == "Y") != (adsl$USUBJID %in% dlt_ae)))
  last_resp <- tapply(seq_len(nrow(adam$ADTR)), adam$ADTR$USUBJID, function(i) { d <- adam$ADTR[i, ]; d$AVALC[which.max(d$ADT)] })
  pd_dc <- adsl$USUBJID[adsl$DCSREAS == "PROGRESSIVE DISEASE"]
  chk("XD03", "ADTR", "Subjects discontinued for PD have PD at last assessment", sum(last_resp[pd_dc] != "PD", na.rm = TRUE))
  chk("XD04", "ADRS", "ORRFL = Y only if a PR visit exists in ADTR",
      sum(adam$ADRS$ORRFL == "Y" & !adam$ADRS$USUBJID %in% adam$ADTR$USUBJID[adam$ADTR$AVALC %in% "PR"]))
  chk("XD05", "ADPP", "PK-evaluable subjects (PKFL = Y) all have ADPP records", sum(!adsl$USUBJID[adsl$PKFL == "Y"] %in% adam$ADPP$USUBJID))
  chk("XD06", "ADAE", "All AEs on or after first dose are flagged TEAE (within safety window)",
      sum(adam$ADAE$ASTDT >= adsl$TRTSDT[match(adam$ADAE$USUBJID, adsl$USUBJID)] &
          adam$ADAE$ASTDT <= adsl$TRTEDT[match(adam$ADAE$USUBJID, adsl$USUBJID)] + cfg$safety_window_days &
          adam$ADAE$TRTEMFL != "Y", na.rm = TRUE))
  res <- do.call(rbind, out)
  for (i in which(res$status == "FAIL")) log_msg("CHECK FAILED: ", res$check_id[i], " ", res$dataset[i], " - ", res$description[i], " (", res$n_fail[i], ")")
  log_msg(sprintf("Data integrity checks: %d run, %d passed, %d failed", nrow(res), sum(res$status == "PASS"), sum(res$status == "FAIL")))
  res
}
