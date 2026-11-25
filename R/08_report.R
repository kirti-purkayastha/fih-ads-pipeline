# ---------------------------------------------------------------------------
# 08_report.R - markdown study summary + run manifest (audit trail)
# ---------------------------------------------------------------------------
md_table <- function(df) {
  df[] <- lapply(df, as.character); df[is.na(df)] <- ""
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

write_report <- function(adam, tabs, qc, checks, cfg, path = "output/study_summary.md") {
  adsl <- adam$ADSL
  slopes <- read.csv("output/tables/dose_proportionality.csv")
  L <- c(
    sprintf("# %s - Exploratory Early-Development Summary (SYNTHETIC DATA)", cfg$study_id), "",
    sprintf("Data cut-off %s. Generated %s by `run_pipeline.R`; every number below is computed from the ADaM datasets in `output/adam/`.",
            cfg$data_cut_date, format(Sys.Date())), "",
    "## Populations", "",
    sprintf("- Enrolled/treated (SAFFL): **%d**; PK-evaluable (PKFL): **%d**; DLT-evaluable: **%d**", sum(adsl$SAFFL == "Y"), sum(adsl$PKFL == "Y"), sum(adsl$DLTEVLFL == "Y")),
    sprintf("- Ongoing at cut-off: **%d**; discontinued: **%d**; subjects with a DLT: **%d**", sum(adsl$EOSSTT == "ONGOING"), sum(adsl$EOSSTT == "DISCONTINUED"), sum(adsl$DLTFL == "Y")), "",
    "## Demographics", "", md_table(tabs$demographics), "",
    "## Safety - treatment-emergent adverse events (subjects, n)", "", md_table(tabs$teae_overview), "",
    "## Pharmacokinetics - geometric mean (geometric CV%)", "",
    "Units: Cmax ng/mL; AUCinf h*ng/mL; t1/2 h; CL L/h. NE = not estimable.", "", md_table(tabs$pk_summary), "",
    sprintf("Dose proportionality (power model, 90%% CI): Cmax slope %.2f (%.2f-%.2f); AUCinf slope %.2f (%.2f-%.2f).",
            slopes$beta[slopes$PARAMCD == "CMAX"], slopes$lo[slopes$PARAMCD == "CMAX"], slopes$hi[slopes$PARAMCD == "CMAX"],
            slopes$beta[slopes$PARAMCD == "AUCIFO"], slopes$lo[slopes$PARAMCD == "AUCIFO"], slopes$hi[slopes$PARAMCD == "AUCIFO"]), "",
    "![PK profiles](figures/01_pk_profiles.png)", "", "![Dose proportionality](figures/02_dose_proportionality.png)", "",
    "## Best overall response (unconfirmed, target lesions only)", "", md_table(tabs$response), "",
    "![Waterfall](figures/04_waterfall.png)", "", "![Swimmer](figures/03_swimmer.png)", "",
    "## Pharmacodynamic and exploratory biomarkers", "", md_table(tabs$cxcl10_peak), "",
    "![CXCL10](figures/05_cxcl10_fold_change.png)", "",
    sprintf("ctDNA was added by protocol amendment 1 (%s): **%d of %d** subjects have a ctDNA baseline; the others contribute post-amendment samples without a baseline and are excluded from change-from-baseline summaries.",
            cfg$amendment1_date, n_distinct(adam$ADBM$USUBJID[adam$ADBM$PARAMCD == "CTDNAVAF" & adam$ADBM$BASEAVL == "Y"]), nrow(adsl)), "",
    "## Quality control", "", "Independent double programming (production vs QC):", "", md_table(qc[, c("dataset", "n_production", "n_qc", "variables_compared", "value_differences", "status")]), "",
    sprintf("Data integrity checks: **%d** run, **%d** passed, **%d** failed (see `output/qc/integrity_checks.csv`).", nrow(checks), sum(checks$status == "PASS"), sum(checks$status == "FAIL")), "",
    "> All data are simulated. No real patient information is used and no clinical conclusions should be drawn.")
  writeLines(L, path)
  log_msg("Report written: ", path)
}

write_manifest <- function(dirs = c("config", "spec", "data/raw", "output/adam", "output/qc", "output/tables", "output/figures"),
                           path = "output/run_manifest.csv") {
  f <- unlist(lapply(dirs, list.files, full.names = TRUE, recursive = TRUE))
  m <- data.frame(file = f, bytes = file.size(f), md5 = unname(tools::md5sum(f)), stringsAsFactors = FALSE)
  write.csv(m, path, row.names = FALSE)
  writeLines(c(capture.output(sessionInfo())), "output/session_info.txt")
  writeLines(.pipeline_env$log, "output/run_log.txt")
  log_msg(sprintf("Manifest written: %d files with MD5 checksums", nrow(m)))
}
