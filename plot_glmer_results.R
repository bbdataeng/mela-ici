plot_glmer_fitted_PCA <- function(
    model, # binomial model fitted with lme4::glmer()
    model_formula, # formula used to fit the model
    results_table, # results_table, output of get_results_table()
    data, # data used to fit the model
    bootstrap_object, # bootstrap object returned by boot.glmm.pred()
    pca_results, # PCA results returned by prcomp()
    link = "logit", # link function
    fitted_resolution = 100, # resolution of calculated fitted values
    PC_to_plot = "PC1", # which PC should be plotted
    show_loadings = TRUE, # add a panel with arrows showing feature loadings for that PC
    file_path = NULL, # optional: path to save the plot
    res_ppi = 300, # resolution (pixels per inch)
    width_in = 8, height_in = 8 # width and height of the plot in inches
    ) {
  # check that the directory exists
  if (!is.null(file_path)) {
    if (!dir.exists(dirname(file_path))) {
      stop( # check that directory exists
        "The specificied directory does not exist."
      )
    }
    # check that file format is supported
    supported_formats <- c("pdf", "png", "jpeg", "jpg", "tiff", "bmp")
    file_format <- tolower(tools::file_ext(file_path))
    if (!file_format %in% supported_formats) {
      stop(
        "File format not supported or missing. Supported formats are: ",
        paste(supported_formats, collapse = ", ")
      )
    }
  }
  # TODO: other if-statements checking correct usage

  # extract model terms
  response_term <- as.character(model_formula)[2]
  predictors <- as.character(model_formula)[3] |>
    strsplit(" \\+ ") |>
    unlist()
  fe_predictors <- grepv("^\\(", predictors, invert = TRUE)
  re_predictors <- grepv("^\\(", predictors, invert = FALSE)

  # obtain fitted values and their confidence intervals
  # over the range of values for that PC
  # with all other model predictors being centered
  fit <- fit.fun(
    m = model, data = data, pred.data = NULL,
    boots = bootstrap_object$all.boots, link = link,
    use = PC_to_plot, # get fitted values over the range of the PC to be plotted
    center = setdiff(fe_predictors, PC_to_plot), # center all other predictors
    resol = fitted_resolution, level = 0.95, offset2add = NULL, keep.boots = FALSE
  )


  # get maximum absolute rotation (for scaling the plot axis)
  xx <- grepv("^PC", fe_predictors) # PCs included in the model formula
  max_rotation <- max(abs(pca_results[["rotation"]][, xx])) * 1.1

  # extract rotation values
  rot <- pca_results[["rotation"]][, PC_to_plot]
  rot <- rot[order(abs(rot), decreasing = TRUE)] # order by decreasing magnitude

  ### prepare stats ###
  r <- results_table[PC_to_plot, , drop = FALSE]
  beta <- as.numeric(r[, "Estimate"])
  lcl <- as.numeric(r[, "lower_CI"])
  ucl <- as.numeric(r[, "upper_CI"])
  OR <- as.numeric(r[, "OR"])
  OR_ci <- exp(c(lcl, ucl))
  chisq <- r[, "Chisq"]
  p_val <- as.character(r[, "p.val"]) # may include stars, e.g. "0.003**"
  # outcome counts
  tab <- table(data[, response_term])
  n_tot <- sum(tab)
  lvl_names <- names(tab)
  # build lines
  lines_txt <- c(
    sprintf("\u03B2 = %.3f  [%.3f, %.3f]", beta, lcl, ucl),
    sprintf("OR = %.2f  [%.2f, %.2f]", OR, OR_ci[1], OR_ci[2]),
    sprintf(
      "LRT: \u03C7\u00B2 = %s, p = %s",
      ifelse(is.na(chisq), "NA", sprintf("%.3f", as.numeric(chisq))), p_val
    ),
    sprintf(
      "N = %d  (%s = %d, %s = %d)",
      n_tot, lvl_names[1], tab[1], lvl_names[2], tab[2]
    )
  )

  # set seed
  set.seed(123)

  # open appropriate graphic device
  if (!is.null(file_path)) {
    if (file_format == "pdf") {
      pdf(file = file_path, width = width_in, height = height_in)
    } else if (file_format == "png") {
      png(
        filename = file_path, width = res_ppi * width_in,
        height = res_ppi * height_in, res = res_ppi
      )
    } else if (file_format == "jpeg" || file_format == "jpg") {
      jpeg(
        filename = file_path, width = res_ppi * width_in,
        height = res_ppi * height_in, res = res_ppi
      )
    } else if (file_format == "tiff") {
      tiff(
        filename = file_path, width = res_ppi * width_in,
        height = res_ppi * height_in, res = res_ppi
      )
    } else if (file_format == "bmp") {
      bmp(
        filename = file_path, width = res_ppi * width_in,
        height = res_ppi * height_in, res = res_ppi
      )
    } else {
      stop("Unknown error when opening the graphic device")
    }
  }

  # graphical parameters
  par(
    las = 1, mfrow = c(2, 1), mar = c(3.5, 4, 0.5, 4),
    mgp = c(2.5, 0.8, 0), tcl = -0.3, xpd = FALSE
  )

  # plot
  plot(
    x = data[, PC_to_plot],
    y = jitter(as.numeric(data[, response_term]) - 1, amount = 0.06),
    ylim = c(-0.2, 1.2), yaxs = "i",
    xlab = "", ylab = "Response",
    yaxt = "n", pch = 21,
    bg = grey(0.2, 0.25), col = adjustcolor("black", 0.4), cex = 1
  )
  mtext(1, text = PC_to_plot, line = 2)
  axis(2, at = 0:1, labels = levels(data[, response_term]))
  axis(4, at = pretty(0:1))
  mtext(4, text = "p(R)", las = 0, line = par("mgp")[1])

  # CIs of the fitted values
  polygon(
    x = c(fit[, PC_to_plot], rev(fit[, PC_to_plot])),
    y = c(fit$lwr, rev(fit$upr)),
    border = NA, col = adjustcolor("#023E8A", 0.2)
  )
  lines(x = fit[, PC_to_plot], y = fit$fit, lwd = 2, col = "#023E8A")



  # bottom plot: PCA loadings
  par(mar = c(2.5, 2, 0.5, 12), mgp = c(1.5, 0.6, 0))
  plot(NULL,
    xlab = paste("PCA loadings for", PC_to_plot), ylab = "",
    xlim = c(-max_rotation, max_rotation),
    ylim = c(length(rot), 0) + 0.5, yaxt = "n", bty = "n", yaxs = "i"
  )
  mtext(side = 2, text = "Feature", las = 0)
  cols <- ifelse(rot >= 0, "#2A9D8F", "#E76F51")
  abline(v = 0, lty = 1)
  arrows(
    x0 = 0, x1 = rot,
    y0 = seq_along(rot), y1 = seq_along(rot),
    length = 0.05, lwd = 1.6, col = cols
  )
  text(
    x = rep(0, length(rot)), y = seq_along(rot),
    pos = ifelse(rot >= 0, 2, 4),
    labels = names(rot),
    cex = scales::rescale(abs(rot), to = c(0.4, 0.8)),
    xpd = NA
  )
  legend("bottomright",
    bty = "n", inset = 0.01, lwd = 2,
    col = c("#2A9D8F", "#E76F51"), legend = c("+", "−")
  )


  # add statistics
  par(xpd = TRUE)
  legend(
    x = mean(c(par("usr")[2], grconvertX(1, from = "ndc", to = "user"))),
    y = mean(par("usr")[3:4]),
    legend = lines_txt, text.col = "black", cex = 0.9, bty = "n",
    xjust = 0.5, yjust = 0.5,
    title = paste("Model results for", PC_to_plot)
  )


  # close graphic device
  if (!is.null(file_path)) dev.off()
}
