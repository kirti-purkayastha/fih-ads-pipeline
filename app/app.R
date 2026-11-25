# ---------------------------------------------------------------------------
# Exploration layer: interactive Shiny app built ONLY on the ADaM outputs.
# Run from the project root after the pipeline:  shiny::runApp("app")
# Every displayed value can be traced back through the ADaM spec to raw data
# (see the "Lineage" tab).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(shiny); library(dplyr); library(ggplot2); library(plotly) })

root <- if (file.exists("output/adam/adsl.rds")) "." else ".."
ad <- function(n) readRDS(file.path(root, "output/adam", paste0(tolower(n), ".rds")))
adsl <- ad("ADSL"); adae <- ad("ADAE"); adpc <- ad("ADPC"); adpp <- ad("ADPP")
adbm <- ad("ADBM"); adtr <- ad("ADTR"); adrs <- ad("ADRS"); adlb <- ad("ADLB")
spec <- read.delim(file.path(root, "spec/adam_spec.psv"), sep = "|", quote = "", comment.char = "", stringsAsFactors = FALSE)
doses <- sort(unique(adsl$TRT01PN))
dfac <- function(x) factor(x, levels = doses, labels = paste(doses, "mg"))

ui <- fluidPage(
  titlePanel("ROC-FIH-SYN-001 - exploratory early-development viewer (synthetic data)"),
  sidebarLayout(
    sidebarPanel(width = 3,
      checkboxGroupInput("dose", "Dose cohorts", choices = doses, selected = doses),
      checkboxInput("logy", "Log scale (PK / biomarker)", TRUE),
      hr(), helpText("All data are simulated. Built exclusively on ADaM datasets in output/adam/.")),
    mainPanel(width = 9, tabsetPanel(
      tabPanel("PK", plotlyOutput("pk", height = 420), br(), tableOutput("ppTab")),
      tabPanel("Biomarkers",
        selectInput("bm", "Analyte", c("CXCL10 (PD)" = "CXCL10", "ctDNA max VAF (exploratory)" = "CTDNAVAF")),
        plotlyOutput("bmPlot", height = 420), uiOutput("bmNote")),
      tabPanel("Response", plotlyOutput("wf", height = 420), plotlyOutput("spider", height = 380)),
      tabPanel("Safety", plotlyOutput("ae", height = 420), tableOutput("aeTab")),
      tabPanel("Labs", selectInput("lab", "Parameter", sort(unique(adlb$PARAMCD)), "ALT"), plotlyOutput("lbPlot", height = 420)),
      tabPanel("Lineage",
        selectInput("lds", "Dataset", sort(unique(spec$dataset))),
        p("Variable-level specification: every derived variable with its source and derivation rule."),
        tableOutput("lin"))))))

server <- function(input, output, session) {
  sel <- reactive(adsl$USUBJID[adsl$TRT01PN %in% as.integer(input$dose)])
  scl <- function(p, lg) if (lg) p + scale_y_log10() else p

  output$pk <- renderPlotly({
    d <- adpc %>% filter(USUBJID %in% sel(), !is.na(AVAL), AVAL > 0) %>% mutate(Dose = dfac(TRT01PN))
    p <- ggplot(d, aes(ARELTM, AVAL, group = USUBJID, colour = Dose, text = USUBJID)) + geom_line(alpha = .6) + geom_point(size = 1) +
      labs(x = "Hours since first dose", y = "DRUGX (ng/mL)") + theme_bw()
    ggplotly(scl(p, input$logy), tooltip = c("text", "x", "y"))
  })
  output$ppTab <- renderTable({
    adpp %>% filter(USUBJID %in% sel(), PARAMCD %in% c("CMAX", "AUCIFO", "THALF", "CL"), !is.na(AVAL)) %>%
      group_by(Dose = paste(TRT01PN, "mg"), Parameter = PARAMCD) %>%
      summarise(n = n(), `Geo mean` = signif(exp(mean(log(AVAL))), 3), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Parameter, values_from = `Geo mean`)
  })

  output$bmPlot <- renderPlotly({
    d <- adbm %>% filter(USUBJID %in% sel(), PARAMCD == input$bm) %>% mutate(Dose = dfac(TRT01PN))
    p <- ggplot(d, aes(ARELTM / 24, AVAL, group = USUBJID, colour = Dose, text = paste(USUBJID, AVISIT))) +
      geom_line(alpha = .5) + geom_point(aes(shape = BLQFL == "Y")) + labs(x = "Days since first dose", y = unique(d$PARAM)[1], shape = "BLQ (LLOQ/2)") + theme_bw()
    ggplotly(scl(p, input$logy), tooltip = c("text", "x", "y"))
  })
  output$bmNote <- renderUI({
    if (input$bm != "CTDNAVAF") return(NULL)
    n0 <- length(unique(adbm$USUBJID[adbm$PARAMCD == "CTDNAVAF" & adbm$BASEAVL == "N"]))
    helpText(sprintf("ctDNA sampling was added by Protocol Amendment 1. %d subject(s) have no ctDNA baseline because they were dosed before the amendment; they are shown but excluded from change-from-baseline derivations.", n0))
  })

  output$wf <- renderPlotly({
    d <- adrs %>% filter(USUBJID %in% sel(), !is.na(BESTPCHG)) %>% arrange(desc(BESTPCHG)) %>%
      mutate(ord = factor(sub(".*-", "", USUBJID), levels = sub(".*-", "", USUBJID)), Dose = dfac(TRT01PN))
    p <- ggplot(d, aes(ord, BESTPCHG, fill = Dose, text = paste(USUBJID, AVALC))) + geom_col() +
      geom_hline(yintercept = c(-30, 20), linetype = "dashed") + labs(x = "Subject", y = "Best % change in target lesions") + theme_bw()
    ggplotly(p, tooltip = c("text", "y"))
  })
  output$spider <- renderPlotly({
    d <- adtr %>% filter(USUBJID %in% sel(), !is.na(PCHG) | ABLFL %in% "Y") %>% mutate(PCHG = ifelse(is.na(PCHG), 0, PCHG), Dose = dfac(TRT01PN))
    p <- ggplot(d, aes(ADY, PCHG, group = USUBJID, colour = Dose, text = USUBJID)) + geom_line(alpha = .6) + geom_point(size = 1) +
      geom_hline(yintercept = c(-30, 20), linetype = "dashed") + labs(x = "Study day", y = "% change from baseline") + theme_bw()
    ggplotly(p, tooltip = c("text", "x", "y"))
  })

  output$ae <- renderPlotly({
    d <- adae %>% filter(USUBJID %in% sel(), TRTEMFL == "Y") %>% count(AEDECOD, Grade = factor(ATOXGRN)) %>%
      group_by(AEDECOD) %>% mutate(tot = sum(n)) %>% ungroup() %>% mutate(AEDECOD = reorder(AEDECOD, tot))
    ggplotly(ggplot(d, aes(AEDECOD, n, fill = Grade)) + geom_col() + coord_flip() + labs(x = NULL, y = "TEAE events") + theme_bw())
  })
  output$aeTab <- renderTable({
    adae %>% filter(USUBJID %in% sel(), TRTEMFL == "Y", DLTFL == "Y" | AESIFL == "Y") %>%
      transmute(Subject = USUBJID, Term = AEDECOD, Grade = ATOXGRN, `Study day` = ASTDY, DLT = DLTFL, AESI = AESIFL)
  })

  output$lbPlot <- renderPlotly({
    d <- adlb %>% filter(USUBJID %in% sel(), PARAMCD == input$lab, !is.na(AVAL)) %>% mutate(Dose = dfac(TRT01PN))
    p <- ggplot(d, aes(ADY, AVAL, group = USUBJID, colour = Dose, text = USUBJID)) + geom_line(alpha = .4) + geom_point(size = .8) +
      geom_hline(aes(yintercept = ANRHI), linetype = "dashed", colour = "grey40") + labs(x = "Study day", y = unique(d$PARAM)[1]) + theme_bw()
    ggplotly(p, tooltip = c("text", "x", "y"))
  })

  output$lin <- renderTable(spec[spec$dataset == input$lds, c("variable", "label", "source", "derivation")])
}

shinyApp(ui, server)
