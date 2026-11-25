# ---------------------------------------------------------------------------
# qc/qc_programs.R - INDEPENDENT QC programs (independent double programming)
#
# Written in base R only, deliberately using different logic and different
# source paths from the production code in R/02-R/04 (production = dplyr,
# vectorised; QC = base R, explicit loops, straight from RAW data).
# Results are compared with PROC COMPARE-style logic in R/05_qc_compare.R.
# ---------------------------------------------------------------------------

qc_adsl <- function(raw, cfg) {
  dm <- raw$dm; ex <- raw$ex; ds <- raw$ds; pc <- raw$pc
  do.call(rbind, lapply(dm$USUBJID, function(u) {
    e <- ex[ex$USUBJID == u, , drop = FALSE]; e <- e[order(e$EXSTDTC), , drop = FALSE]
    first <- as.Date(substr(e$EXSTDTC[1], 1, 10)); last <- as.Date(substr(e$EXSTDTC[nrow(e)], 1, 10))
    t0 <- as.POSIXct(e$EXSTDTC[1], format = "%Y-%m-%dT%H:%M", tz = "UTC")
    p  <- pc[pc$USUBJID == u, , drop = FALSE]
    tt <- as.POSIXct(p$PCDTC, format = "%Y-%m-%dT%H:%M", tz = "UTC")
    num <- suppressWarnings(as.numeric(p$PCORRES))           # "<50" and "" become NA
    arm <- dm$ARM[dm$USUBJID == u]; armcd <- dm$ARMCD[dm$USUBJID == u]
    d <- ds[ds$USUBJID == u, , drop = FALSE]
    data.frame(USUBJID = u, TRTSDT = first, TRTEDT = last, NDOSES = nrow(e),
               TRTDURD = as.integer(last - first) + 1L, EOSDT = as.Date(d$DSSTDTC),
               EOSSTT = if (d$DSDECOD == "ONGOING") "ONGOING" else if (d$DSDECOD == "STUDY TREATMENT COMPLETED") "COMPLETED" else "DISCONTINUED",
               SAFFL = if (nrow(e) > 0) "Y" else "N",
               PKFL = if (any(tt > t0 & !is.na(num), na.rm = TRUE)) "Y" else "N",
               TRT01PN = as.integer(strsplit(arm, " ")[[1]][2]),
               COHORT = as.integer(substring(armcd, 2)),
               AGEGR1 = if (dm$AGE[dm$USUBJID == u] >= 65) ">=65" else "<65",
               stringsAsFactors = FALSE)
  }))
}

# independent NCA: vectorised segment AUCs and closed-form least squares
qc_nca <- function(t, c, dose, npts, r2min) {
  t <- pmax(t, 0); o <- order(t); t <- t[o]; c <- c[o]
  tl <- max(which(c > 0)); tt <- t[seq_len(tl)]; cc <- c[seq_len(tl)]
  c1 <- cc[-tl]; c2 <- cc[-1]; dt <- diff(tt)
  seg <- ifelse(c2 < c1 & c2 > 0, (c1 - c2) * dt / (log(c1) - log(c2)), (c1 + c2) * dt / 2)
  res <- c(CMAX = max(c), TMAX = t[which(c == max(c))[1]], AUCLST = sum(seg),
           LAMZ = NA, THALF = NA, AUCIFO = NA, CL = NA)
  pts <- which(t > res[["TMAX"]] & c > 0 & seq_along(c) <= tl)
  if (length(pts) >= npts) {
    pts <- pts[(length(pts) - npts + 1):length(pts)]
    x <- t[pts]; y <- log(c[pts])
    slope <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
    r2 <- cor(x, y)^2; r2adj <- 1 - (1 - r2) * (npts - 1) / (npts - 2)
    if (slope < 0 && r2adj >= r2min) {
      res["LAMZ"] <- -slope; res["THALF"] <- log(2) / -slope
      res["AUCIFO"] <- res[["AUCLST"]] + c[tl] / -slope
      res["CL"] <- dose * 1e6 / res[["AUCIFO"]] / 1000
    }
  }
  res
}

qc_adpp <- function(raw, cfg) {
  dm <- raw$dm; ex <- raw$ex; pc <- raw$pc
  out <- list()
  for (u in dm$USUBJID) {
    p <- pc[pc$USUBJID == u & nzchar(pc$PCORRES), , drop = FALSE]
    e <- ex[ex$USUBJID == u, , drop = FALSE]; e <- e[order(e$EXSTDTC), , drop = FALSE]
    t0 <- as.POSIXct(e$EXSTDTC[1], format = "%Y-%m-%dT%H:%M", tz = "UTC")
    tt <- as.numeric(difftime(as.POSIXct(p$PCDTC, format = "%Y-%m-%dT%H:%M", tz = "UTC"), t0, units = "hours"))
    tt <- round(tt, 2)                                        # spec: ARELTM rounded to 2 dp
    cc <- ifelse(grepl("^<", p$PCORRES), 0, suppressWarnings(as.numeric(p$PCORRES)))
    if (!any(cc > 0, na.rm = TRUE)) next
    if (!any(tt > 0 & cc > 0)) next                       # PK population rule
    r <- qc_nca(tt, cc, e$EXDOSE[1], cfg$pk$lamz_min_points, cfg$pk$lamz_min_r2adj)
    out[[u]] <- data.frame(USUBJID = u, PARAMCD = names(r), AVAL = unname(r), stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

qc_adrs <- function(raw, cfg) {
  tr <- raw$tr; tp <- cfg$tumor
  do.call(rbind, lapply(unique(raw$dm$USUBJID), function(u) {
    x <- tr[tr$USUBJID == u, , drop = FALSE]
    base <- x$TRORRES[x$VISIT == "BASELINE"]
    post <- x[x$VISIT != "BASELINE", , drop = FALSE]; post <- post[order(post$TRDTC), , drop = FALSE]
    resp <- character(0); pch <- numeric(0); prev <- base
    for (k in seq_len(nrow(post))) {
      v <- post$TRORRES[k]; nad <- min(prev)
      pd <- post$NEWLES[k] == "Y" || (100 * (v - nad) / nad >= tp$pd_threshold_pct - 1e-9 && (v - nad) >= tp$pd_min_abs_mm - 1e-9)
      pc <- 100 * (v - base) / base
      resp <- c(resp, if (pd) "PD" else if (pc <= tp$pr_threshold_pct + 1e-9) "PR" else "SD")
      pch <- c(pch, pc); prev <- c(prev, v)
    }
    bor <- if (!length(resp)) "NE" else if ("PR" %in% resp) "PR" else if ("SD" %in% resp) "SD" else "PD"
    data.frame(USUBJID = u, AVALC = bor, NASSESS = length(resp),
               BESTPCHG = if (length(pch)) min(pch) else NA_real_,
               ORRFL = if (bor == "PR") "Y" else "N", DCRFL = if (bor %in% c("PR", "SD")) "Y" else "N",
               stringsAsFactors = FALSE)
  }))
}
