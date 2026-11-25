# Run from the project root:  Rscript tests/testthat.R
library(testthat)
for (f in c(sort(list.files("R", pattern = "^0[0-8].*\\.R$", full.names = TRUE)),
            list.files("R/qc", pattern = "\\.R$", full.names = TRUE))) source(f)
res <- test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\n%d tests, %d failed, %d errors\n", sum(df$nb), sum(df$failed), sum(df$error)))
if (sum(df$failed) + sum(df$error) > 0) quit(status = 1)
