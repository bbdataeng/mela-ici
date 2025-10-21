# Function to make stacked barplots ---------------------------------------

make_stacked_barplot <- function(
    x_var, # x-axis variable
    x_var_name = "", # name of the x-axis variable
    col_var, # color-axis variable
    col_var_name = "", # name of the color-axis variable
    type = c("frequency", "proportion"),
    color_palette, # vector of colors matching levels of col_var
    rotate_x_var_labels = FALSE, # whether x_axis labels should be rotated (useful for long names)
    show_x_var_frequencies = FALSE, # don't show x-var frequencies
    show_percentages = TRUE, # whether to show percentages in the bars
    min_percentage = 2, # hide percentages below this threshold
    file_path = NULL, # file path (with format) to save the plot
    res_ppi = 300, # resolution of the plot
    width_in = 6, height_in = 6, # width and height of the plot (in inches)
    ... # additional arguments passed to par()
    ) {
  # TODO:
  # if statements to check valid arguments

  # create table
  xtab <- table(col_var, x_var)

  # get percentages by levels of x_var
  if (type == "proportion") {
    xtab <- apply(X = xtab, MARGIN = 2, FUN = function(x) 100 * x / sum(x))
  }

  # open graphic device
  if (!is.null(file_path)) {
    file_format <- tolower(tools::file_ext(file_path)) # extract file format

    # open appropriate graphic device
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
      stop("File format not supported.")
    }
  } else {
    old_par <- par(no.readonly = TRUE) # save old par() settings if not saving plot
    on.exit(par(old_par), add = TRUE) # restore old par() settings
  }

  # prepare graphical parameters
  dots <- list(...)
  if (length(dots)) do.call(par, dots) # only call par() if something was passed
  par(xpd = TRUE)

  # get x-axis values
  x_vals <- barplot(xtab, plot = FALSE)

  # create barplot
  barplot(
    xtab,
    axisnames = FALSE, # hide axis names here
    col = color_palette, # colors
    xlab = "", # hide x-axis label here
    ylab = ifelse(type == "proportion", "Percentage", "Frequency"),
    legend.text = TRUE, # add legend
    args.legend = list( # additional settings for the legend
      title = col_var_name, # title of the legend
      bty = "n", ncol = 1, cex = 1,
      x = par("usr")[2], y = mean(par("usr")[3:4]),
      xjust = 0, yjust = 0.5
    )
  )

  # add x-axis labels
  if (rotate_x_var_labels) {
    text(
      labels = colnames(xtab), x = x_vals,
      y = par()$usr[3] - (par()$usr[4] - par()$usr[3]) * 0.035,
      srt = 45, cex = 0.8, adj = 1
    )
  } else {
    mtext(
      text = colnames(xtab),
      at = x_vals,
      side = 1,
      line = par("mgp")[2]
    )
  }

  # add x-axis name
  mtext(x_var_name, side = 1, line = par("mar")[1] - 1.5)

  # add percentages
  if (show_percentages && type == "frequency") {
    for (i_x in seq_len(ncol(xtab))) {
      x_level <- colnames(xtab)[i_x]
      for (i_y in seq_len(nrow(xtab))) {
        y_level <- rownames(xtab)[i_y]
        x_percent <- (xtab[y_level, x_level] / sum(xtab[, x_level])) * 100
        if (x_percent > min_percentage) {
          text(
            y = sum(xtab[1:i_y - 1, i_x]) + xtab[i_y, i_x] / 2,
            x = x_vals[i_x],
            adj = c(0.5, 0.5), cex = 0.5,
            labels = paste0(round(x_percent, 1), "%")
          )
        }
      }
    }
  } else if (show_percentages && type == "proportion") {
    for (i_x in seq_len(ncol(xtab))) {
      x_level <- colnames(xtab)[i_x]
      for (i_y in seq_len(nrow(xtab))) {
        y_level <- rownames(xtab)[i_y]
        if (xtab[y_level, x_level] > min_percentage) {
          text(
            y = sum(xtab[1:i_y - 1, i_x]) + xtab[i_y, i_x] / 2,
            x = x_vals[i_x],
            adj = c(0.5, 0.5), cex = 0.5,
            labels = paste0(round(xtab[y_level, x_level], 1), "%")
          )
        }
      }
    }
  }

  # add x_var frequencies
  if (show_x_var_frequencies) {
    mtext(
      side = 3,
      at = x_vals,
      text = table(x_var),
      line = 0.5
    )
    mtext(side = 3, text = paste(x_var_name, "frequency"), line = 2)
  }

  # close graphic device
  if (!is.null(file_path)) dev.off()
}



# Function to make dot plots ----------------------------------------------

make_dotplot <- function(
    x_var, # x-axis variable
    x_var_name = "", # name of the x-axis variable
    y_var, # y-axis variable
    y_var_name = "", # name of the y-axis variable
    expansion_factor = 5, # regulate dot size
    pt.pch = 16, # type of point
    pt.col = grey(0.5, 0.5), # color of the points
    show_legend = TRUE, # if FALSE, hide legend
    rotate_x_var_labels = FALSE, # whether x_axis labels should be rotated (useful for long names)
    file_path = NULL, # file path (with format) to save the plot
    res_ppi = 300, # resolution of the plot
    width_in = 6, height_in = 6, # width and height of the plot (in inches)
    ... # additional arguments passed to par()
    ) {
  # TODO:
  # if statements to check valid arguments

  # create table
  xtab <- table(y_var = y_var, x_var = x_var) |> as.data.frame()

  # get x and y coordinates
  xtab$x <- as.numeric(xtab$x_var)
  xtab$y <- as.numeric(xtab$y_var)

  # calculate dot size
  xtab$cex <- sqrt(xtab$Freq / max(xtab$Freq)) * expansion_factor

  # open graphic device
  if (!is.null(file_path)) {
    file_format <- tolower(tools::file_ext(file_path)) # extract file format

    # open appropriate graphic device
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
      stop("File format not supported.")
    }
  } else {
    old_par <- par(no.readonly = TRUE) # save old par() settings if not saving plot
    on.exit(par(old_par), add = TRUE) # restore old par() settings
  }

  # prepare graphical parameters
  dots <- list(...)
  if (length(dots)) do.call(par, dots) # only call par() if something was passed
  par(xpd = FALSE)

  # empty plot
  plot(NULL,
    xlim = c(1, nlevels(x_var)) + c(-0.5, 0.5),
    ylim = c(nlevels(y_var), 1) + c(0.5, -0.5),
    xaxs = "i", yaxs = "i", bty = "l",
    xlab = "", ylab = "", xaxt = "n", yaxt = "n"
  )

  # add grid
  abline(v = 1:nlevels(x_var), col = "gray70", lty = "dotted")
  abline(h = 1:nlevels(y_var), col = "gray70", lty = "dotted")

  # add points
  par(xpd = TRUE)
  points(x = xtab$x, y = xtab$y, cex = xtab$cex, pch = pt.pch, col = pt.col)

  # add axes
  axis(side = 2, at = 1:nlevels(y_var), labels = levels(y_var))
  mtext(y_var_name, side = 2, line = par("mar")[2] - 1.5, las = 0)
  axis(side = 1, at = 1:nlevels(x_var), labels = FALSE) # labels added later
  mtext(x_var_name, side = 1, line = par("mar")[1] - 1.5, las = 0)

  # add x-axis labels
  if (rotate_x_var_labels) {
    text(
      labels = levels(x_var), x = 1:nlevels(x_var),
      y = par()$usr[3] - (par()$usr[4] - par()$usr[3]) * 0.035,
      srt = 45, cex = 0.8, adj = 1
    )
  } else {
    mtext(
      text = levels(x_var),
      at = 1:nlevels(x_var),
      side = 1,
      line = par("mgp")[2]
    )
  }

  # add legend
  legend_labels <- c(pretty(xtab$Freq), 1) |>
    unique() |>
    sort()
  legend_pt.cx <- sqrt(legend_labels / max(xtab$Freq)) * expansion_factor
  legend(
    x = mean(c(par("usr")[2], grconvertX(1, from = "ndc", to = "user"))),
    y = mean(c(grconvertY(0, from = "ndc", to = "user"), grconvertY(1, from = "ndc", to = "user"))),
    # y = mean(par("usr")[3:4]),
    xjust = 0.5, yjust = 0.5, bty = "n", cex = 0.8,
    legend = legend_labels, pt.cex = legend_pt.cx, col = pt.col,
    pch = pt.pch, y.intersp = 2, x.intersp = 2, title = "N. observations"
  )

  # close graphic device
  if (!is.null(file_path)) dev.off()
}
