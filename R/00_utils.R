# ---------------------------------------------------------------------------
# 00_utils.R - shared helpers: config, logging, date handling, spec handling
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(yaml); library(jsonlite); library(ggplot2)
})

.pipeline_env <- new.env()
.pipeline_env$log <- character()

load_config <- function(path = "config/study_config.yml") yaml::read_yaml(path)

log_msg <- function(...) {
  m <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  message(m)
  .pipeline_env$log <- c(.pipeline_env$log, m)
  invisible(m)
}

# ---- date helpers ---------------------------------------------------------
to_date <- function(x) {
  x <- substr(x, 1, 10)
  x[is.na(x) | nchar(x) < 10] <- NA_character_
  as.Date(x)
}
to_dtm <- function(x) {
  x[is.na(x) | !nzchar(x)] <- NA_character_
  as.POSIXct(x, format = "%Y-%m-%dT%H:%M", tz = "UTC")
}
# CDISC study-day convention: there is no day 0
study_day <- function(date, ref) {
  d <- as.integer(date - ref)
  ifelse(is.na(d), NA_integer_, ifelse(d >= 0L, d + 1L, d))
}
# TRUE for the single latest candidate record (by date) within a group
flag_last <- function(date, cand) {
  out <- rep(FALSE, length(date))
  cand[is.na(cand)] <- FALSE
  if (any(cand)) {
    m <- max(date[cand]); i <- which(cand & date == m)
    out[i[length(i)]] <- TRUE
  }
  out
}
# Partial-date imputation for AE start (SAP 6.2):
#   YYYY-MM -> first of month, but never earlier than first dose when same month
impute_start <- function(dtc, ref_start) {
  n <- nchar(dtc); full <- n >= 10; ym <- n == 7
  d <- rep(as.Date(NA), length(dtc))
  d[full] <- as.Date(substr(dtc[full], 1, 10))
  d[ym]   <- as.Date(paste0(dtc[ym], "-01"))
  idx <- ym & !is.na(ref_start) & format(d, "%Y-%m") == format(ref_start, "%Y-%m")
  d[idx] <- ref_start[idx]
  list(date = d, flag = ifelse(full, NA_character_, ifelse(ym, "D", "M")))
}

# ---- specification handling ----------------------------------------------
read_spec <- function(path = "spec/adam_spec.psv") {
  read.delim(path, sep = "|", stringsAsFactors = FALSE, quote = "", comment.char = "",
             strip.white = TRUE)
}

coerce_type <- function(x, type) {
  switch(type,
    text = as.character(x), integer = as.integer(x), float = as.numeric(x),
    date = as.Date(x), datetime = x, x)
}

# Enforce the specification: same variables, same order, types and labels.
apply_spec <- function(df, ds, spec) {
  sp <- spec[spec$dataset == ds, , drop = FALSE]
  if (!nrow(sp)) stop(ds, ": no specification found")
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  miss  <- setdiff(sp$variable, names(df))
  extra <- setdiff(names(df), sp$variable)
  if (length(miss))  stop(ds, ": in spec but not derived: ", paste(miss, collapse = ", "))
  if (length(extra)) stop(ds, ": derived but not in spec: ", paste(extra, collapse = ", "))
  df <- df[, sp$variable, drop = FALSE]
  for (i in seq_len(nrow(sp))) {
    v <- sp$variable[i]
    df[[v]] <- coerce_type(df[[v]], sp$type[i])
    attr(df[[v]], "label") <- sp$label[i]
  }
  df
}

write_adam <- function(df, name, spec, dir = "output/adam") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  df <- apply_spec(df, name, spec)
  saveRDS(df, file.path(dir, paste0(tolower(name), ".rds")))
  write.csv(df, file.path(dir, paste0(tolower(name), ".csv")), row.names = FALSE, na = "")
  log_msg(sprintf("%s written: %d rows x %d cols", name, nrow(df), ncol(df)))
  invisible(df)
}

read_adam <- function(name, dir = "output/adam") readRDS(file.path(dir, paste0(tolower(name), ".rds")))

# ---- raw data reader ------------------------------------------------------
read_raw <- function(dir = "data/raw") {
  num_cols <- c("AGE","ECOG","EXDOSE","EXSEQ","PCSEQ","PCLLOQ","LBSEQ","LBORRES","LBSTNRLO",
                "LBSTNRHI","VISITNUM","TRSEQ","TRORRES","AESEQ","AETOXGR")
  rd <- function(f) {
    x <- read.csv(file.path(dir, f), colClasses = "character", stringsAsFactors = FALSE,
                  na.strings = character(0))
    for (v in intersect(num_cols, names(x))) x[[v]] <- suppressWarnings(as.numeric(x[[v]]))
    x
  }
  list(dm = rd("dm.csv"), ex = rd("ex.csv"), ds = rd("ds.csv"), pc = rd("pc.csv"),
       ae = rd("ae.csv"), lb = rd("lb.csv"), tr = rd("tr.csv"),
       bm_path = file.path(dir, "bm.jsonl"))
}
