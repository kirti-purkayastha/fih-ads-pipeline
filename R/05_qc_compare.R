# ---------------------------------------------------------------------------
# 05_qc_compare.R - PROC COMPARE-style comparison of production vs QC datasets
# ---------------------------------------------------------------------------
qc_compare <- function(prod, qc, keys, name, vars = NULL, tol = 1e-6) {
  prod <- as.data.frame(prod); qc <- as.data.frame(qc)
  if (is.null(vars)) vars <- setdiff(intersect(names(prod), names(qc)), keys)
  kp <- do.call(paste, c(prod[keys], sep = "|")); kq <- do.call(paste, c(qc[keys], sep = "|"))
  only_p <- setdiff(kp, kq); only_q <- setdiff(kq, kp); common <- intersect(kp, kq)
  pm <- prod[match(common, kp), , drop = FALSE]; qm <- qc[match(common, kq), , drop = FALSE]
  diffs <- list()
  for (v in vars) {
    a <- pm[[v]]; b <- qm[[v]]
    if (is.numeric(a) && is.numeric(b)) {
      same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol * pmax(1, abs(b)))
    } else {
      same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & as.character(a) == as.character(b))
    }
    bad <- !same | is.na(same)
    if (any(bad)) diffs[[v]] <- data.frame(dataset = name, variable = v, key = common[bad],
                                           production = as.character(a[bad]), qc = as.character(b[bad]))
  }
  nval <- if (length(diffs)) sum(vapply(diffs, nrow, 1L)) else 0L
  list(summary = data.frame(dataset = name, n_production = nrow(prod), n_qc = nrow(qc),
         keys_only_production = length(only_p), keys_only_qc = length(only_q),
         variables_compared = length(vars), variables_with_differences = length(diffs),
         value_differences = nval,
         status = ifelse(!length(only_p) && !length(only_q) && !length(diffs), "MATCH", "MISMATCH")),
       diffs = if (length(diffs)) do.call(rbind, diffs) else NULL)
}

run_qc <- function(raw, cfg, adam, out_dir = "output/qc") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  q_adsl <- qc_adsl(raw, cfg); q_adpp <- qc_adpp(raw, cfg); q_adrs <- qc_adrs(raw, cfg)
  pp_vars <- c("CMAX", "TMAX", "AUCLST", "AUCIFO", "LAMZ", "THALF", "CL")
  res <- list(
    qc_compare(adam$ADSL, q_adsl, "USUBJID", "ADSL"),
    qc_compare(adam$ADPP[adam$ADPP$PARAMCD %in% pp_vars & !is.na(adam$ADPP$AVAL), ],
               q_adpp[!is.na(q_adpp$AVAL), ], c("USUBJID", "PARAMCD"), "ADPP (NCA)"),
    qc_compare(adam$ADRS, q_adrs, "USUBJID", "ADRS"))
  summ <- do.call(rbind, lapply(res, `[[`, "summary"))
  diffs <- do.call(rbind, lapply(res, `[[`, "diffs"))
  write.csv(summ, file.path(out_dir, "qc_compare_summary.csv"), row.names = FALSE)
  write.csv(if (is.null(diffs)) data.frame() else diffs, file.path(out_dir, "qc_compare_differences.csv"), row.names = FALSE)
  for (i in seq_len(nrow(summ))) log_msg(sprintf("QC compare %-12s %s (%d obs, %d vars)", summ$dataset[i],
                                                  summ$status[i], summ$n_production[i], summ$variables_compared[i]))
  summ
}
