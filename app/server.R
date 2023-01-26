library(shiny)
library(shinyFeedback)
library(glue)
library(tidyverse)
library(markdown)
library(phacking)

# ------------------------------------------------------------------------------
# helper functions for formatting
# ------------------------------------------------------------------------------

.str <- function(s) {
  paste(strwrap(glue(s, .envir = parent.frame())), collapse = " ")
}

ci_text <- function(estimate, ci_lower, ci_upper, sig = 2) {
  .str("{signif(estimate, sig)} (95% CI [{signif(ci_lower, sig)},
        {signif(ci_upper, sig)}])")
}
estimate_text <- function(model_label, model_result, sig = 2) {
  if (is.null(model_result)) ci <- ""
  else ci <- ci_text(model_result$estimate, model_result$ci_lower,
                     model_result$ci_upper, sig = sig)
  p(strong(glue("{str_to_sentence(model_label)} estimate:")), br(), ci)
}

danger <- function(inputId, show, text) {
  feedbackDanger(inputId, show, text, color = "#e74c3c", icon = NULL)
}

warn <- function(inputId, show, text) {
  feedbackWarning(inputId, show, text, color = "#f39c12", icon = NULL)
}

shinyServer(function(input, output) {
  
  # ----------------------------------------------------------------------------
  # overall input elements
  # ----------------------------------------------------------------------------
  
  output$y_cols <- renderUI({
    # req(input$meta_data)
    req(meta_data_raw())
    selectInput("y_col", "Column of point estimates",
                choices = c("Select a column" = "", names(meta_data_raw())))
  })
  
  output$v_cols <- renderUI({
    # req(input$meta_data)
    req(meta_data_raw())
    selectInput("v_col", "Column of estimated variances",
                choices = c("Select a column" = "", names(meta_data_raw())))
  })
  
  output$directions <- renderUI({
    # req(input$meta_data, input$v_col)
    req(meta_data_raw())
    selectInput("direction", "Direction",
                choices = c("favor positive", "favor negative"))
  })
  
  # ----------------------------------------------------------------------------
  # reactive values based on overall inputs
  # ----------------------------------------------------------------------------
  
  meta_data_raw <- reactive({
    # meta_file <- input$meta_data
    # read_csv(meta_file$datapath, show_col_types = FALSE)
    phacking::money_priming_meta
  })
  
  meta_data <- reactive({
    # req(input$meta_data, input$y_col, input$v_col)
    req(meta_data_raw(), input$y_col, input$v_col)
    # phacking::money_priming_meta
    meta_data_raw() |>
      filter(!is.na(.data[[input$y_col]]), !is.na(.data[[input$v_col]]))
  })
  
  positive <- reactive({
    req(input$direction)
    str_detect(input$direction, "positive")
  })
  
  y_vals <- reactive({
    req(meta_data(), input$y_col)
    meta_data()[[input$y_col]]
  })
  
  v_vals <- reactive({
    req(meta_data(), input$v_col)
    meta_data()[[input$v_col]]
  })
  
  # ----------------------------------------------------------------------------
  # input validation
  # ----------------------------------------------------------------------------
  
  valid_y <- reactive({
    req(y_vals())
    y_valid <- is.numeric(y_vals())
    danger("y_col", !y_valid, "values must be numeric")
    req(y_valid)
  })
  
  valid_v <- reactive({
    req(v_vals())
    v_valid <- is.numeric(v_vals()) & all(v_vals() > 0)
    danger("v_col", !v_valid, "values must be numeric & positive")
    req(v_valid)
  })
  
  valid_affirm <- reactive({
    req(y_vals(), v_vals(), input$direction)
    
    if (positive()) yi = y_vals() else yi = -y_vals()
    # TODO: could this not duplicate affirm calculation?
    pvals <- 2 * (1 - pnorm(abs(yi) / sqrt(v_vals())))
    alpha <- formals(phacking::phacking_meta)$alpha_select
    affirm <- (pvals < alpha) & (yi > 0)
    no_aff <- sum(affirm) == 0
    no_nonaff <- sum(!affirm) == 0
    no_either <- no_aff | no_nonaff
    no_dir <- if (no_aff) "affirmative" else if (no_nonaff) "nonaffirmative"
    error <- .str("There are zero {no_dir} studies – double check your columns
                  and direction.")
    danger("error", no_either, error)
    req(!no_either)
  })
  
  # ----------------------------------------------------------------------------
  # phacking_meta
  # ----------------------------------------------------------------------------
  
  # TODO: include? correct?
  uncorrected_model <- reactive({
    req(valid_y(), valid_v(), valid_affirm())
    robu_formula <- as.formula(glue("{input$y_col} ~ 1"))
    meta_model <- robumeta::robu(robu_formula,
                                 # studynum = cluster_col(),
                                 studynum = 1:nrow(meta_data()),
                                 data = meta_data(),
                                 var.eff.size = v_vals(),
                                 small = TRUE)
    meta_result <- metabias::robu_ci(meta_model)
    
    opposite_dir <- meta_result$estimate < 0 & positive() |
      meta_result$estimate > 0 & !positive()
    warn("error", opposite_dir,
         "Warning: favored direction is opposite of the pooled estimate.")
    meta_result
  })
  
  output$uncorrected <- renderUI({
    req(uncorrected_model())
    estimate_text("uncorrected", uncorrected_model())
  })
  
  corrected_model <- reactive({
    req(valid_y(), valid_v(), valid_affirm())
    meta <- phacking_meta(yi = meta_data()[[input$y_col]],
                          vi = meta_data()[[input$v_col]],
                          favor_positive = positive())
    meta$stats |> rename(estimate = mode)
  })
  
  output$corrected_mu <- renderUI({
    req(corrected_model())
    estimate_text("corrected mean (μ)", corrected_model() |> filter(param == "mu"))
  })

  output$corrected_tau <- renderUI({
    req(corrected_model())
    estimate_text("corrected heterogeneity (τ)", corrected_model() |> filter(param == "tau"))
  })
  
  corrected_summary <- reactive({
    req(corrected_model())
    # cm <- corrected_model()
    mu <- corrected_model() |> filter(param == "mu")
    tau <- corrected_model() |> filter(param == "tau")
    .str("Accounting for potential <em>p</em>-hacking and publication bias that
          favor affirmative results, the estimated meta-analytic mean (μ) is
          {ci_text(mu$estimate, mu$ci_lower, mu$ci_upper)} and the estimated standard
          deviation of the effects, i.e., heterogeneity (τ), is
          {ci_text(tau$estimate, tau$ci_lower, tau$ci_upper)}.")
  })
  
  output$corrected_summary <- renderUI({
    req(corrected_summary())
    p(em(HTML(corrected_summary())))
  })
  
  output$clip_corrected <- renderUI({
    req(corrected_summary())
    rclipButton(
      inputId = "clipbtn_corrected",
      label = "Copy summary",
      clipText = corrected_summary(), 
      icon = icon("clipboard")
    )
  })
  
  # ----------------------------------------------------------------------------
  # qqplot
  # ----------------------------------------------------------------------------
  
  # funnel_plot <- function() {
  #   significance_funnel(yi = y_vals(), vi = v_vals(),
  #                       favor_positive = positive(),
  #                       est_all = uncorrected_model()$estimate,
  #                       est_worst = worst_model()$estimate) +
  #     theme_classic(base_family = "Lato") +
  #     theme(legend.position = "top",
  #           legend.title = element_blank())
  # }
  # 
  # fp_res <- 300
  # fp_width <- 1200
  # fp_height <- 1100
  # 
  # output$funnel <- renderPlot({
  #   req(uncorrected_model(), worst_model())
  #   funnel_plot()
  # }, res = fp_res, height = fp_height, width = fp_width)
  # 
  # output$download_funnel <- downloadHandler(
  #   filename = function() {
  #     paste0(tools::file_path_sans_ext(input$meta_data$name), "_funnel", ".png")
  #   },
  #   content = function(file) {
  #     ggsave(file, plot = funnel_plot(), device = "png", dpi = fp_res,
  #            height = fp_height, width = fp_width, units = "px")
  #   }
  # )
  # 
  # output$download_funnel_button <- renderUI({
  #   req(uncorrected_model(), worst_model())
  #   downloadButton("download_funnel")
  # })
  
})
