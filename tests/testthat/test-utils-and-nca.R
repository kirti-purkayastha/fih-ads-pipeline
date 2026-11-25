# Unit tests for helper functions and the NCA routine (known analytic answers)

test_that("study_day follows the CDISC convention (no day 0)", {
  ref <- as.Date("2025-03-03")
  expect_equal(study_day(as.Date("2025-03-03"), ref), 1L)
  expect_equal(study_day(as.Date("2025-03-04"), ref), 2L)
  expect_equal(study_day(as.Date("2025-03-02"), ref), -1L)
  expect_true(is.na(study_day(as.Date(NA), ref)))
})

test_that("partial AE start dates are imputed per SAP 6.2", {
  ref <- as.Date("2025-03-10")
  r <- impute_start(c("2025-05-17", "2025-05", "2025-03"), rep(ref, 3))
  expect_equal(r$date, as.Date(c("2025-05-17", "2025-05-01", "2025-03-10")))
  expect_equal(r$flag, c(NA, "D", "D"))
})

test_that("flag_last picks exactly one latest candidate, ignoring non-candidates", {
  d <- as.Date("2025-01-01") + c(0, 5, 10, 15)
  expect_equal(flag_last(d, c(TRUE, TRUE, TRUE, FALSE)), c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(sum(flag_last(d, rep(FALSE, 4))), 0)
})

test_that("NCA recovers exact parameters for a mono-exponential IV profile", {
  k <- 0.05; dose <- 100; c0 <- 1000                       # ng/mL
  t <- c(0, 1, 2, 4, 8, 12, 24, 48)
  conc <- c0 * exp(-k * t)
  r <- suppressWarnings(nca_profile(t, conc, dose))   # perfect fit warning is expected here
  expect_equal(unname(r["LAMZ"]), k, tolerance = 1e-8)
  expect_equal(unname(r["THALF"]), log(2) / k, tolerance = 1e-8)
  expect_equal(unname(r["CMAX"]), c0)
  # log-down trapezoid is exact for mono-exponential decay
  expect_equal(unname(r["AUCLST"]), (c0 - conc[length(conc)]) / k, tolerance = 1e-8)
  expect_equal(unname(r["AUCIFO"]), c0 / k, tolerance = 1e-8)
  expect_equal(unname(r["CL"]), dose * 1e6 / (c0 / k) / 1000, tolerance = 1e-8)
})

test_that("NCA reports not-estimable rather than guessing when the terminal phase is unusable", {
  r <- nca_profile(c(0, 1, 4), c(0, 500, 300), 10)         # too few post-Tmax points
  expect_true(is.na(r["LAMZ"])); expect_true(is.na(r["AUCIFO"]))
  r2 <- nca_profile(c(0, 1, 4), c(0, 0, 0), 10)            # all BLQ
  expect_true(is.na(r2["CMAX"]))
})

test_that("independent QC NCA agrees with production NCA on a noisy profile", {
  set.seed(1)
  t <- c(0, 1, 4, 24, 72, 168, 336); conc <- 800 * exp(-0.02 * t) * exp(rnorm(length(t), 0, 0.05))
  a <- nca_profile(t, conc, 50); b <- qc_nca(t, conc, 50, 3, 0.8)
  for (p in c("CMAX", "TMAX", "AUCLST", "AUCIFO", "THALF", "CL"))
    expect_equal(unname(a[p]), unname(b[p]), tolerance = 1e-8, label = p)
})
