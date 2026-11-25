if (basename(getwd()) == "testthat") setwd("../..")   # tests run from the project root

# Integration tests: build ADaM from raw in a temp copy, then prove the QC machinery
# actually catches errors (a QC process that never fails is not evidence of anything).

build_all <- function(raw, cfg) {
  adsl0 <- build_adsl_base(raw, cfg); adae <- build_adae(raw, adsl0, cfg)
  adsl <- adsl_add_ae_flags(adsl0, adae, cfg); adpc <- build_adpc(raw, adsl, cfg)
  adam <- list(ADSL = adsl, ADAE = adae, ADLB = build_adlb(raw, adsl, cfg), ADPC = adpc,
               ADPP = build_adpp(adpc, adsl, cfg), ADBM = build_adbm(raw, adsl, cfg),
               ADTR = build_adtr(raw, adsl, cfg))
  adam$ADRS <- build_adrs(adam$ADTR, adsl)
  spec <- read_spec()
  Map(function(d, n) apply_spec(d, n, spec), adam, names(adam))
}

cfg <- load_config(); spec <- read_spec()
tmp <- tempfile(); dir.create(tmp); suppressMessages(simulate_raw(cfg, tmp))
raw <- read_raw(tmp); adam <- build_all(raw, cfg)

test_that("all integrity checks pass on a clean run", {
  ck <- suppressMessages(run_checks(adam, raw, spec, cfg))
  expect_gte(nrow(ck), 50)
  expect_equal(sum(ck$status == "FAIL"), 0)
})

test_that("production and independent QC datasets match", {
  qc <- suppressMessages(run_qc(raw, cfg, adam, out_dir = file.path(tmp, "qc")))
  expect_true(all(qc$status == "MATCH"))
})

test_that("pipeline is deterministic for a fixed seed", {
  tmp2 <- tempfile(); dir.create(tmp2); suppressMessages(simulate_raw(cfg, tmp2))
  f <- list.files(tmp, full.names = FALSE, pattern = "\\.(csv|jsonl)$")
  expect_equal(unname(tools::md5sum(file.path(tmp, f))), unname(tools::md5sum(file.path(tmp2, f))))
})

test_that("the specification is enforced (extra or missing variables stop the build)", {
  bad <- adam$ADSL; bad$EXTRA <- 1
  expect_error(apply_spec(bad, "ADSL", spec), "not in spec")
  expect_error(apply_spec(adam$ADSL[, -1], "ADSL", spec), "not derived")
})

test_that("seeded errors are detected by the integrity checks", {
  a <- adam; a$ADLB$AVAL[a$ADLB$PARAMCD == "HGB"][1] <- 13.5          # unit error (g/dL value left in g/L field)
  ck <- suppressMessages(run_checks(a, raw, spec, cfg)); expect_true("DER05" %in% ck$check_id[ck$status == "FAIL"])

  a <- adam; a$ADSL$TRTEDT[1] <- a$ADSL$TRTSDT[1] - 5                    # date logic
  ck <- suppressMessages(run_checks(a, raw, spec, cfg)); expect_true("DT01" %in% ck$check_id[ck$status == "FAIL"])

  a <- adam; a$ADAE <- rbind(a$ADAE, a$ADAE[1, ])                        # duplicate record
  ck <- suppressMessages(run_checks(a, raw, spec, cfg)); expect_true("KEY01" %in% ck$check_id[ck$status == "FAIL"])

})

test_that("seeded errors are detected by independent double programming", {
  a <- adam; a$ADPP$AVAL[a$ADPP$PARAMCD == "AUCIFO"][1] <- a$ADPP$AVAL[a$ADPP$PARAMCD == "AUCIFO"][1] * 1.01
  qc <- suppressMessages(run_qc(raw, cfg, a, out_dir = file.path(tmp, "qc2")))
  expect_equal(qc$status[qc$dataset == "ADPP (NCA)"], "MISMATCH")
  a <- adam; a$ADSL$PKFL[1] <- ifelse(a$ADSL$PKFL[1] == "Y", "N", "Y")
  qc <- suppressMessages(run_qc(raw, cfg, a, out_dir = file.path(tmp, "qc3")))
  expect_equal(qc$status[qc$dataset == "ADSL"], "MISMATCH")
  # a self-consistent derivation error (integrity checks cannot see it) is caught only by double programming
  a <- adam; a$ADRS$AVALC[a$ADRS$AVALC == "PR"][1] <- "SD"
  qc <- suppressMessages(run_qc(raw, cfg, a, out_dir = file.path(tmp, "qc4")))
  expect_equal(qc$status[qc$dataset == "ADRS"], "MISMATCH")
})

test_that("amendment-driven missing baselines are handled, not hidden", {
  ct <- adam$ADBM[adam$ADBM$PARAMCD == "CTDNAVAF", ]
  expect_true(any(ct$BASEAVL == "N"))                                    # pre-amendment subjects exist
  expect_true(all(is.na(ct$CHG[ct$BASEAVL == "N"])))                     # and get no change-from-baseline
  expect_true(all(as.Date(ct$ADTM) >= as.Date(cfg$amendment1_date)))
})

test_that("unit harmonisation removes the site-specific haemoglobin unit difference", {
  raw_hgb <- raw$lb[raw$lb$LBTESTCD == "HGB", ]
  expect_setequal(unique(raw_hgb$LBORRESU), c("g/L", "g/dL"))            # raw data really are mixed
  expect_equal(unique(adam$ADLB$AVALU[adam$ADLB$PARAMCD == "HGB"]), "g/L")
})
