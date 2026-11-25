# Reproducible environment. Pin the R version; packages come from the Posit Public Package Manager
# snapshot date baked into the base image (rocker/r-ver), so the build is repeatable.
FROM rocker/r-ver:4.3.3
RUN install2.r --error dplyr tidyr purrr tibble ggplot2 plotly shiny testthat yaml jsonlite scales
WORKDIR /project
COPY . /project
CMD ["Rscript", "run_pipeline.R", "--regen"]
