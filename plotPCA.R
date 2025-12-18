##### function to plot PCA data, with colors coding a factor #####
plot_pca_factor <- function(
    pca_res, # output of prcomp()
    PCs_to_plot = c("PC1", "PC2"), # PCs to plot on X and Y axis respectively
    metadata, # dataframe of metadata (rows matching PCA data)
    vars_to_use, # which variables (columns of metadata) should be used?
    colors_vars_list, # list of colors for the variables
    par_mar = c(2.5, 2.5, 2.5, 5), # graphical parameters par("mar")
    dark_theme = FALSE, # set to false to have white theme
    pt.cex = 1, # point size
    equal_axis_scale = FALSE, # set to TRUE to force same scaling of the X and Y axes
    legend_cex = 1, # relative dimension of the legend
    legend_ncols = 1, # columns to split the legend in
    output_folder, # folder to save the plots
    file_format = "png", # file format
    res_ppi = 300, # resolution (pixels per inch)
    width_in = 8, height_in = 8 # width and height of the plot in inches
    ) {
  ##### quality check #####

  # check that metadata is a dataframe
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a dataframe.")
  }
  # check that metadata has the right number of rows
  if (nrow(metadata) != nrow(pca_res$x)) {
    stop("Different number of rows in 'metadata' and 'pca_res$x'")
  }
  # check that vars_to_use is a character
  if (!is.character(vars_to_use)) {
    stop("'vars_to_use' must be a single character or vector of characters.")
  }
  # check that vars_to_use are all present in the names of metadata
  if (!all(vars_to_use %in% names(metadata))) {
    xx <- vars_to_use[which(!vars_to_use %in% names(metadata))]
    stop("Missing columns in 'metadata': ", paste(xx, collapse = ", "))
  }
  # check that PCs_to_plot are valid
  if (!all(PCs_to_plot %in% colnames(pca_res$x))) {
    stop("Unvalid PCs specified in 'PCs_to_plot'.")
  }
  # check that the directory exists
  if (!dir.exists(output_folder)) {
    stop("The specified directory does not exist.")
  }
  # check that file format is supported
  supported_formats <- c("pdf", "png", "jpeg", "jpg", "tiff", "bmp")
  if (!file_format %in% supported_formats) {
    stop(
      "File format not supported. supported formats are: ",
      paste(supported_formats, collapse = ", ")
    )
  }

  # set columns of metadata as factors if not already
  for (i in ncol(metadata)) {
    if (!is.factor(metadata[, i])) {
      metadata[, i] <- as.factor(metadata[, i])
      message("Transformed ", names(metadata)[i], " into a factor with levels in alphabetical order")
    }
  }

  # check that colors_vars_list is a named list
  if (!is.list(colors_vars_list) | is.null(names(colors_vars_list))) {
    stop("'colors_vars_list' must be a named list.")
  }
  # check that the names of colors_vars_list correspond to vars_to_use
  if (!all(vars_to_use %in% names(colors_vars_list))) {
    xx <- vars_to_use[which(!vars_to_use %in% names(colors_vars_list))]
    stop("Missing variables in 'colors_vars_list': ", paste(xx, collapse = ", "))
  }
  # check that colors are ok
  for (xvar in vars_to_use) {
    # check length
    if (length(colors_vars_list[[xvar]]) != nlevels(metadata[, xvar])) {
      stop(
        "The length of 'colors_vars_list$", xvar,
        "' must equal the number of levels of 'metadata$", xvar, "'."
      )
    }
    # check if all entries are valid colors
    valid <- tryCatch(
      {
        col2rgb(colors_vars_list[[xvar]])
        TRUE
      },
      error = function(e) FALSE
    )
    if (!valid) {
      stop(
        "'colors_vars_list$", xvar,
        "' must contain valid R color names or hex codes."
      )
    }
    # check that all names of the colors correspond to levels of the variable
    if (!all(names(colors_vars_list[[xvar]]) %in% levels(metadata[, xvar]))) {
      stop(
        "The names of 'colors_vars_list$", xvar,
        "' don't match the levels of 'metadata$", xvar, "'."
      )
    }
  }

  ##### data preparation #####

  # extract PCs
  xpcs <- pca_res$x

  # calculate % of variation in each PC
  percent_var <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)
  names(percent_var) <- colnames(xpcs)

  # reorder colors in the same order as the levels in metadata
  for (xvar in vars_to_use) {
    colors_vars_list[[xvar]] <- colors_vars_list[[xvar]][levels(metadata[, xvar])]
  }

  ##### create plot for each variable to use #####
  for (xvar in vars_to_use) {
    # open appropriate graphic device
    filename <- file.path(output_folder, paste0(
      "pca_", xvar, "_", paste0(PCs_to_plot, collapse = "-")
    ))
    if (file_format == "pdf") {
      pdf(
        file = paste(filename, file_format, sep = "."),
        width = width_in, height = height_in
      )
    } else if (file_format == "png") {
      png(
        filename = paste(filename, file_format, sep = "."),
        width = res_ppi * width_in, height = res_ppi * height_in, res = res_ppi
      )
    } else if (file_format == "jpeg" | file_format == "jpg") {
      jpeg(
        filename = paste(filename, file_format, sep = "."),
        width = res_ppi * width_in, height = res_ppi * height_in, res = res_ppi
      )
    } else if (file_format == "tiff") {
      tiff(
        filename = paste(filename, file_format, sep = "."),
        width = res_ppi * width_in, height = res_ppi * height_in, res = res_ppi
      )
    } else if (file_format == "bmp") {
      bmp(
        filename = paste(filename, file_format, sep = "."),
        width = res_ppi * width_in, height = res_ppi * height_in, res = res_ppi
      )
    } else {
      stop("Unknown error when opening the graphic device")
    }

    # set dark theme if appropriate
    if (dark_theme) {
      par(
        bg = "black", fg = "white", col = "white", col.axis = "white",
        col.lab = "white", col.main = "white", col.sub = "white"
      )
    }

    # set graphical parameters
    par(mar = par_mar, mgp = c(1.5, 0.4, 0), tcl = -0.2, las = 1, xpd = FALSE)

    # create plot
    plot(
      x = xpcs[, PCs_to_plot[1]],
      y = xpcs[, PCs_to_plot[2]],
      pch = 21, cex = pt.cex, bty = "o", axes = TRUE,
      asp = ifelse(equal_axis_scale, 1, NA),
      ylab = paste0(
        PCs_to_plot[2], ": ", percent_var[PCs_to_plot[2]],
        "% of total variation"
      ),
      xlab = paste0(
        PCs_to_plot[1], ": ", percent_var[PCs_to_plot[1]],
        "% of total variation"
      ),
      main = paste0("PCA by ", xvar),
      bg = colors_vars_list[[xvar]][as.numeric(metadata[, xvar])]
    )

    # add legend
    par(xpd = TRUE)
    legend(
      x = mean(c(grconvertX(1, from = "ndc", to = "user"), par("usr")[2])),
      y = mean(par("usr")[3:4]),
      xjust = 0.5, yjust = 0.5,
      legend = levels(metadata[, xvar]),
      pt.bg = colors_vars_list[[xvar]], pch = 21, pt.cex = pt.cex,
      border = NA, cex = legend_cex, ncol = legend_ncols,
      title = xvar, bty = "n"
    )

    # close device
    dev.off()
  }
}
