# ---------------------------------------------------------------------------
# 07_outputs.R - exploratory figures and summary tables (all driven by ADaM)
# ---------------------------------------------------------------------------
geo_mean <- function(x) { x <- x[!is.na(x) & x > 0]; if (!length(x)) NA_real_ else exp(mean(log(x))) }
geo_cv   <- function(x) { x <- x[!is.na(x) & x > 0]; if (length(x) < 2) NA_real_ else 100 * sqrt(exp(var(log(x))) - 1) }
dose_fac <- function(x) factor(x, levels = sort(unique(x)), labels = paste(sort(unique(x)), "mg"))
theme_study <- function() theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
save_fig <- function(p, name, w = 8, h = 5, dir = "output/figures") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, name), p, width = w, height = h, dpi = 150, bg = "white")
}

make_figures <- function(adam, cfg) {
  adsl <- adam$ADSL
  # 1. Geometric-mean PK profiles ------------------------------------------------
  pk <- adam$ADPC %>% filter(!is.na(AVAL), ATPTN > 0) %>%
    group_by(TRT01PN, ATPTN) %>% summarise(gm = geo_mean(AVAL), n = sum(AVAL > 0), .groups = "drop") %>%
    filter(!is.na(gm)) %>% mutate(Dose = dose_fac(TRT01PN))
  p1 <- ggplot(pk, aes(ATPTN, gm, colour = Dose, group = Dose)) +
    geom_hline(yintercept = cfg$pk$lloq_ng_ml, linetype = "dashed", colour = "grey40") +
    geom_line() + geom_point(size = 1.8) + scale_y_log10(labels = scales::label_comma()) +
    annotate("text", x = max(pk$ATPTN), y = cfg$pk$lloq_ng_ml, label = "LLOQ", vjust = -0.6, hjust = 1, size = 3) +
    labs(title = "DRUGX serum concentration - geometric mean by dose (Cycle 1)",
         subtitle = "Synthetic data. BLQ values (0) excluded from geometric means.",
         x = "Nominal time since start of infusion (h)", y = "Concentration (ng/mL, log scale)") + theme_study()
  save_fig(p1, "01_pk_profiles.png")

  # 2. Dose proportionality (power model) ----------------------------------------
  dp <- adam$ADPP %>% filter(PARAMCD %in% c("CMAX", "AUCIFO"), !is.na(AVAL), AVAL > 0)
  slopes <- dp %>% group_by(PARAMCD) %>% group_modify(~ {
    f <- lm(log(AVAL) ~ log(TRT01PN), data = .x); ci <- confint(f, "log(TRT01PN)", level = 0.90)
    tibble(beta = coef(f)[2], lo = ci[1], hi = ci[2], lab = sprintf("slope = %.2f (90%% CI %.2f-%.2f)", coef(f)[2], ci[1], ci[2]))
  })
  write.csv(slopes %>% select(-lab), "output/tables/dose_proportionality.csv", row.names = FALSE)
  dp <- dp %>% mutate(Parameter = recode(PARAMCD, CMAX = "Cmax (ng/mL)", AUCIFO = "AUCinf (h*ng/mL)"))
  sl <- slopes %>% mutate(Parameter = recode(PARAMCD, CMAX = "Cmax (ng/mL)", AUCIFO = "AUCinf (h*ng/mL)"))
  p2 <- ggplot(dp, aes(TRT01PN, AVAL)) + geom_point(aes(colour = dose_fac(TRT01PN)), alpha = .7) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", linewidth = .6) +
    geom_text(data = sl, aes(x = min(dp$TRT01PN), y = Inf, label = lab), hjust = 0, vjust = 1.5, size = 3, inherit.aes = FALSE) +
    scale_x_log10() + scale_y_log10(labels = scales::label_comma()) + facet_wrap(~Parameter, scales = "free_y") +
    labs(title = "Dose proportionality (power model)", subtitle = "AUC slope > 1 indicates saturable (target-mediated) clearance; Cmax is ~proportional (Cmax ~ dose / volume). Synthetic data.",
         x = "Dose (mg, log scale)", y = "Exposure (log scale)", colour = "Dose") + theme_study()
  save_fig(p2, "02_dose_proportionality.png", w = 9, h = 4.5)

  # 3. Swimmer plot ---------------------------------------------------------------
  sw <- adsl %>% mutate(dur = as.integer(EOSDT - TRTSDT) + 1L, Dose = dose_fac(TRT01PN)) %>%
    arrange(COHORT, dur) %>% mutate(row = row_number())
  ev_resp <- adam$ADTR %>% filter(!is.na(AVALC)) %>% inner_join(sw %>% select(USUBJID, row), by = "USUBJID")
  ev_ae <- adam$ADAE %>% filter(TRTEMFL == "Y", ATOXGRN >= 3) %>% inner_join(sw %>% select(USUBJID, row), by = "USUBJID")
  p3 <- ggplot(sw) +
    geom_segment(aes(x = 0, xend = dur, y = row, yend = row, colour = Dose), linewidth = 2.6, alpha = .75) +
    geom_point(data = ev_resp, aes(ADY, row, shape = AVALC, fill = AVALC), size = 2.3, colour = "black") +
    geom_point(data = ev_ae, aes(ASTDY, row + 0.38), shape = 25, size = 1.8, fill = "red", colour = "red") +
    geom_text(data = sw %>% filter(EOSSTT == "ONGOING"), aes(dur + 4, row, label = ">"), size = 3.5) +
    scale_shape_manual(values = c(PR = 21, SD = 22, PD = 24), name = "Visit response") +
    scale_fill_manual(values = c(PR = "#2ca02c", SD = "#bdbdbd", PD = "#d62728"), name = "Visit response") +
    scale_y_continuous(breaks = sw$row, labels = sub(".*-", "", sw$USUBJID), expand = expansion(add = 1)) +
    labs(title = "Swimmer plot - treatment duration, visit responses and grade >= 3 TEAEs",
         subtitle = "Red triangles = grade >= 3 TEAE; '>' = ongoing at data cut-off. Synthetic data.",
         x = "Study day", y = "Subject") + theme_study() + theme(axis.text.y = element_text(size = 6))
  save_fig(p3, "03_swimmer.png", w = 9, h = 7)

  # 4. Waterfall ---------------------------------------------------------------------
  wf <- adam$ADRS %>% filter(!is.na(BESTPCHG)) %>% arrange(desc(BESTPCHG)) %>%
    mutate(ord = row_number(), Dose = dose_fac(TRT01PN))
  p4 <- ggplot(wf, aes(ord, BESTPCHG, fill = Dose)) + geom_col() +
    geom_hline(yintercept = c(-30, 20), linetype = "dashed") +
    geom_text(aes(label = AVALC, y = ifelse(BESTPCHG >= 0, BESTPCHG + 3, BESTPCHG - 3)), size = 2.6) +
    labs(title = "Best percentage change in sum of target-lesion diameters",
         subtitle = "Dashed lines: -30% (PR) and +20% (PD). Responses unconfirmed. Synthetic data.",
         x = "Subjects (ordered)", y = "Best change from baseline (%)") + theme_study() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  save_fig(p4, "04_waterfall.png")

  # 5. CXCL10 PD biomarker --------------------------------------------------------------
  cx <- adam$ADBM %>% filter(PARAMCD == "CXCL10", !is.na(R2BASE)) %>%
    group_by(TRT01PN, AVISIT) %>% summarise(t = median(ARELTM), fold = geo_mean(R2BASE), .groups = "drop") %>%
    mutate(Dose = dose_fac(TRT01PN))
  p5 <- ggplot(cx, aes(t, fold, colour = Dose, group = Dose)) + geom_hline(yintercept = 1, linetype = "dashed") +
    geom_line() + geom_point() + scale_y_log10() +
    labs(title = "CXCL10 fold change from baseline (pharmacodynamic response)",
         subtitle = "Geometric mean ratio to pre-dose; synthetic data", x = "Hours since first dose",
         y = "Fold change (log scale)") + theme_study()
  save_fig(p5, "05_cxcl10_fold_change.png")

  # 6. ctDNA vs tumour change (evolving endpoint; only subjects with a ctDNA baseline) ---------
  ct <- adam$ADBM %>% filter(PARAMCD == "CTDNAVAF", !is.na(PCHG)) %>% select(USUBJID, AVISIT, TRT01PN, ct_pchg = PCHG)
  tm <- adam$ADTR %>% filter(!is.na(PCHG)) %>% select(USUBJID, AVISIT, tum_pchg = PCHG)
  j <- inner_join(ct, tm, by = c("USUBJID", "AVISIT")) %>% mutate(Dose = dose_fac(TRT01PN))
  if (nrow(j) > 2) {
    p6 <- ggplot(j, aes(tum_pchg, ct_pchg)) + geom_point(aes(colour = Dose)) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = .5) +
      labs(title = "ctDNA change vs tumour change (matched visits)",
           subtitle = sprintf("n = %d matched visits from %d subjects with an available ctDNA baseline (post-amendment enrolment)", nrow(j), n_distinct(j$USUBJID)),
           x = "Tumour sum of diameters, % change", y = "ctDNA max VAF, % change") + theme_study()
    save_fig(p6, "06_ctdna_vs_tumour.png")
  }
  invisible(TRUE)
}

make_tables <- function(adam, cfg, dir = "output/tables") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  adsl <- adam$ADSL; adae <- adam$ADAE
  # demographics by cohort
  demog <- adsl %>% group_by(Dose = dose_fac(TRT01PN)) %>% summarise(
    N = n(), `Age, median (range)` = sprintf("%.0f (%d-%d)", median(AGE), min(AGE), max(AGE)),
    `Male, n` = sum(SEX == "M"), `ECOG 1, n` = sum(ECOG == 1), .groups = "drop")
  tot <- adsl %>% summarise(Dose = "Total", N = n(), `Age, median (range)` = sprintf("%.0f (%d-%d)", median(AGE), min(AGE), max(AGE)),
                            `Male, n` = sum(SEX == "M"), `ECOG 1, n` = sum(ECOG == 1))
  demog <- bind_rows(mutate(demog, Dose = as.character(Dose)), tot)
  # TEAE overview by cohort
  te <- adae %>% filter(TRTEMFL == "Y")
  by_subj <- function(f) vapply(split(adsl$USUBJID, dose_fac(adsl$TRT01PN)), function(u) length(unique(te$USUBJID[f(te) & te$USUBJID %in% u])), 1L)
  teae <- tibble(Dose = levels(dose_fac(adsl$TRT01PN)), N = as.integer(table(dose_fac(adsl$TRT01PN))),
    `Any TEAE` = by_subj(function(d) rep(TRUE, nrow(d))), `Grade >= 3` = by_subj(function(d) d$ATOXGRN >= 3),
    `Related` = by_subj(function(d) d$RELFL == "Y"), `Serious` = by_subj(function(d) d$AESERFL == "Y"),
    `CRS / IRR (AESI)` = by_subj(function(d) d$AESIFL == "Y"), `DLT` = by_subj(function(d) d$DLTFL == "Y"))
  # PK summary
  pk <- adam$ADPP %>% filter(PARAMCD %in% c("CMAX", "AUCIFO", "THALF", "CL")) %>%
    group_by(Dose = dose_fac(TRT01PN), PARAMCD) %>%
    summarise(n = sum(!is.na(AVAL)), gm = geo_mean(AVAL), cv = geo_cv(AVAL), .groups = "drop") %>%
    mutate(txt = ifelse(n == 0, "NE", sprintf("%s (%.0f%%)", formatC(gm, digits = 3, format = "fg", big.mark = ","), cv))) %>%
    select(Dose, PARAMCD, txt) %>% tidyr::pivot_wider(names_from = PARAMCD, values_from = txt) %>%
    mutate(Dose = as.character(Dose)) %>% select(Dose, CMAX, AUCIFO, THALF, CL)
  # response
  resp <- adam$ADRS %>% group_by(Dose = dose_fac(TRT01PN)) %>% summarise(
    N = n(), PR = sum(AVALC == "PR"), SD = sum(AVALC == "SD"), PD = sum(AVALC == "PD"), NE = sum(AVALC == "NE"),
    `ORR, %` = round(100 * mean(ORRFL == "Y"), 0), `DCR, %` = round(100 * mean(DCRFL == "Y"), 0), .groups = "drop") %>%
    mutate(Dose = as.character(Dose))
  # CXCL10 peak fold change
  pd <- adam$ADBM %>% filter(PARAMCD == "CXCL10", !is.na(R2BASE)) %>% group_by(USUBJID, TRT01PN) %>%
    summarise(peak = max(R2BASE), .groups = "drop") %>% group_by(Dose = dose_fac(TRT01PN)) %>%
    summarise(n = n(), `Peak CXCL10 fold change (geo mean)` = round(geo_mean(peak), 2), .groups = "drop") %>%
    mutate(Dose = as.character(Dose))
  tabs <- list(demographics = demog, teae_overview = teae, pk_summary = pk, response = resp, cxcl10_peak = pd)
  for (n in names(tabs)) write.csv(tabs[[n]], file.path(dir, paste0(n, ".csv")), row.names = FALSE)
  tabs
}
