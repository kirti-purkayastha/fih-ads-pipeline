#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# run_pipeline.R - end-to-end, one-command reproducible run
#   Rscript run_pipeline.R            # use existing data/raw (or simulate if empty)
#   Rscript run_pipeline.R --regen    # regenerate synthetic raw data first
# Exit status is non-zero if QC compare or integrity checks fail (CI gate).
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
for (f in c(sort(list.files("R", pattern = "^0[0-9].*\\.R$", full.names = TRUE)),
            list.files("R/qc", pattern = "\\.R$", full.names = TRUE))) source(f)

cfg  <- load_config(); spec <- read_spec()
log_msg("Pipeline start | study ", cfg$study_id, " | seed ", cfg$seed)

if ("--regen" %in% args || !file.exists("data/raw/dm.csv")) simulate_raw(cfg)
raw <- read_raw()

# ---- ADaM derivations ---------------------------------------------------------
adsl0 <- build_adsl_base(raw, cfg)
adae  <- build_adae(raw, adsl0, cfg)
adsl  <- adsl_add_ae_flags(adsl0, adae, cfg)
adpc  <- build_adpc(raw, adsl, cfg)
adam <- list(ADSL = adsl, ADAE = adae, ADLB = build_adlb(raw, adsl, cfg), ADPC = adpc,
             ADPP = build_adpp(adpc, adsl, cfg), ADBM = build_adbm(raw, adsl, cfg),
             ADTR = build_adtr(raw, adsl, cfg))
adam$ADRS <- build_adrs(adam$ADTR, adsl)
adam <- Map(function(d, n) write_adam(d, n, spec), adam, names(adam))

# ---- QC ------------------------------------------------------------------------
qc <- run_qc(raw, cfg, adam)
checks <- run_checks(adam, raw, spec, cfg)
dir.create("output/qc", recursive = TRUE, showWarnings = FALSE)
write.csv(checks, "output/qc/integrity_checks.csv", row.names = FALSE)

# ---- outputs ----------------------------------------------------------------------
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
tabs <- make_tables(adam, cfg)
make_figures(adam, cfg)
write_report(adam, tabs, qc, checks, cfg)
write_manifest()

if (any(qc$status != "MATCH") || any(checks$status != "PASS")) {
  log_msg("PIPELINE FAILED QUALITY GATES"); quit(status = 1)
}
log_msg("Pipeline finished: all quality gates passed")
