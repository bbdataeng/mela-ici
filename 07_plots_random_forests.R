# load("nonsync/07_plots_random_forests.RData") # run to restore working space

# Load data ---------------------------------------------------------------
load("nonsync/06_random_forests.RData") # load random forests data

# Load libraries ----------------------------------------------------------
library(ComplexHeatmap) # v2.26.0
library(grid)
library(paletteer)
library(cvms)
library(ggplot2)

# Plot settings -----------------------------------------------------------

resol <- 300 # resolution in ppi

# prepare colors for response (6 levels)
colors_response6 <- paletteer_c("grDevices::RdYlBu", nlevels(alldata$response_6levels))
names(colors_response6) <- levels(alldata$response_6levels)

# prepare colors for response (3 levels)
colors_response3 <- paletteer_c("grDevices::RdYlBu", nlevels(alldata$response_3levels))
names(colors_response3) <- levels(alldata$response_3levels)

# prepare colors for response (2 levels)
colors_response2 <- paletteer_c("grDevices::RdYlBu", nlevels(alldata$response_2levels))
names(colors_response2) <- levels(alldata$response_2levels)

# prepare colors for RF formulae
rf_formulae_colors <- paletteer_d("RColorBrewer::Paired", 8)
rf_formulae_colors <- rf_formulae_colors[!as.logical(seq_along(rf_formulae_colors) %% 2)]
names(rf_formulae_colors) <- paste0("RF", 1:4)
# add reference in to_do
to_do$rf_short <- substr(to_do$rf_formula, start = 1, stop = 3) |> factor()


# Make barplots of importance scores --------------------------------------

for (i in seq_len(nrow(to_do))) {
  png(
    file.path(to_do$folder[i], paste0(
      "vimp_barplot_", to_do$rf_formula[i], ".png"
    )),
    width = 8 * resol, height = 8 * resol, res = resol
  )
  par(mar = c(5.1, 12.5, 1.5, 0.5), las = 2)
  xx <- barplot(importances_df[i, ], horiz = TRUE, plot = FALSE)
  xlim <- range(importances_df[i, ], na.rm = TRUE)
  xlim <- xlim + c(-diff(xlim) * 0.1, diff(xlim) * 0.1)
  varnames <- ifelse(names(importances_df[i, ]) %in% names(cell_types_original),
    cell_types_original[names(importances_df[i, ])],
    names(importances_df[i, ])
  )
  plot(NULL,
    xlim = xlim, ylim = range(xx), axes = FALSE,
    xlab = "Variable Importance", ylab = "",
    main = to_do$rf_formula[i]
  )
  barplot(importances_df[i, ],
    horiz = TRUE, add = TRUE, names.arg = varnames,
    col = rf_formulae_colors[as.numeric(to_do$rf_short[i])]
  )
  text(
    x = as.numeric(importances_df[i, ]), y = xx, labels = importances_rank[[i]],
    pos = ifelse(importances_df[i, ] >= 0 | is.na(importances_df[i, ]), 4, 2)
  )
  dev.off()
  # make barplot of top 10 importance scores
  png(
    file.path(to_do$folder[i], paste0(
      "vimp_barplot_top10_", to_do$rf_formula[i], ".png"
    )),
    width = 6 * resol, height = 4 * resol, res = resol
  )
  par(mar = c(5.1, 12.5, 1.5, 0.5), las = 2)
  top10 <- importances_df[i, ] |>
    sort(decreasing = TRUE) |>
    head(10)
  xx <- barplot(top10, horiz = TRUE, plot = FALSE)
  xlim <- c(0, max(top10) * 1.1)
  varnames <- ifelse(names(top10) %in% names(cell_types_original),
    cell_types_original[names(top10)],
    names(top10)
  )
  plot(NULL,
    xlim = xlim, ylim = rev(range(xx)) + c(0.5, 0), axes = FALSE,
    xlab = "Variable Importance", ylab = "",
    main = to_do$rf_formula[i]
  )
  barplot(top10,
    horiz = TRUE, add = TRUE,
    names.arg = varnames,
    col = rf_formulae_colors[as.numeric(to_do$rf_short[i])]
  )
  text(x = as.numeric(top10), y = xx, labels = 1:10, pos = 4)
  dev.off()
}


# Compare variable importance across models -------------------------------

# data frame of variable importance ranks
vip_rank <- apply(importances_df, MARGIN = 1, FUN = function(x) {
  xx <- rank(x, na.last = TRUE)
  xx[is.na(x)] <- NA
  xx <- max(xx, na.rm = TRUE) + 1 - xx
  return(xx)
})
vip_rank <- vip_rank[order(rownames(vip_rank)), ] # order alphabetically
rownames(vip_rank) <- ifelse(
  rownames(vip_rank) %in% names(cell_types_original),
  cell_types_original[rownames(vip_rank)],
  rownames(vip_rank)
)

# prepare annotation dataframe
ann_df <- data.frame(
  to_do[, c("response_type", "rf_short")],
  row.names = colnames(vip_rank)
)
ann_df$response_levels <- ifelse(
  ann_df$response_type == "binary", 2, ifelse(
    ann_df$response_type == "ordinal3", 3, 6
  )
)
ann_df$rf_class <- ifelse(ann_df$response_levels == 2, "ranger", "ordfor")
ann_df <- ann_df[, c("rf_class", "response_levels", "rf_short")]
names(ann_df) <- c("RF Object Class", "N. Response Levels", "RF Formula")

# prepare annotation colors
palettes <- c(
  "ggsci::default_nejm", "ggthemes::Green", "RColorBrewer::Paired"
)
names(palettes) <- names(ann_df)
ann_colors <- lapply(names(ann_df), FUN = function(x) {
  xvar <- ann_df[, x]
  if (is.logical(xvar) | is.character(xvar)) xvar <- factor(xvar)
  if (is.factor(xvar)) {
    cols <- paletteer_d(palettes[x], nlevels(xvar))
    names(cols) <- levels(xvar)
  } else {
    cols <- paletteer_c(palettes[x], length(unique(xvar)))
    names(cols) <- sort(unique(xvar))
  }
  return(cols)
})
names(ann_colors) <- names(ann_df)
# change colors for rf formula as previously defined
ann_colors$`RF Formula` <- rf_formulae_colors

# create temporary heatmap and revert row order
ht <- Heatmap(vip_rank) |> draw()
new_order <- row_dend(ht) |> rev()

# create heatmap
png(file.path(output_folder, "heatmap_variable_importance.png"),
  width = 8 * resol, height = 6 * resol, res = resol
)
Heatmap(
  vip_rank,
  name = "Ranked Variable Importance\n(1 = most important)",
  col = paletteer_c("viridis::viridis", 150, -1),
  show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(
    show_legend = TRUE,
    df = ann_df,
    col = ann_colors,
    annotation_name_gp = gpar(fontface = "bold")
  ),
  cluster_rows = new_order,
  cluster_columns = FALSE,
  row_title = "Variables",
  column_title = "Random Forest models",
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    legend_width = unit(4, "cm") # adjust to taste
  )
) |> draw(
  merge_legend = TRUE, # all legends in one column
  heatmap_legend_side = "right", # keep everything on the right
  annotation_legend_side = "right"
)
dev.off()


# Compare accuracy metrics for binary RFs ---------------------------------

# see metrics of ranger (binary) random forests
names(accuracy_metrics[[1]])

# put together metrics of binary random forests
accuracy_ranger_wide <- cbind(
  to_do[
    to_do$response_type == "binary",
    c("rf_short", "technical_predictors", "age_and_gender")
  ],
  sapply(
    X = which(to_do$response_type == "binary"),
    FUN = function(x) accuracy_metrics[[x]]
  ) |>
    t() |> as.data.frame()
)

# reshape into long format
accuracy_ranger_long <- reshape(
  data = accuracy_ranger_wide,
  varying = c(
    "accuracy", "sensitivity", "specificity", "precision",
    "f1_score", "balanced_accuracy", "mcc"
  ),
  v.names = "metric",
  idvar = "rf_short",
  times = c(
    "accuracy", "sensitivity", "specificity", "precision",
    "f1_score", "balanced_accuracy", "mcc"
  ),
  timevar = "metric_type",
  direction = "long"
)
accuracy_ranger_long$metric_type <- factor( # metric type as factor
  accuracy_ranger_long$metric_type,
  levels = c(
    "accuracy", "sensitivity", "specificity", "precision",
    "f1_score", "balanced_accuracy", "mcc"
  )
)
levels(accuracy_ranger_long$metric_type) <- c(
  "Accuracy", "Sensitivity", "Specificity", "Precision",
  "F1 score", "Balanced\nAccuracy", "MCC"
)

# prepare coordinates for barplot
xx <- barplot(
  metric ~ rf_short + metric_type,
  data = accuracy_ranger_long,
  beside = TRUE, plot = FALSE
)
rownames(xx) <- levels(accuracy_ranger_long$rf_short)
colnames(xx) <- levels(accuracy_ranger_long$metric_type)
# create barplot
png(file.path(output_folder, "barplot_metrics_binary.png"),
  width = 9 * resol, height = 3.5 * resol, res = resol
)
par(las = 1, xpd = TRUE, tcl = -0.3, mar = c(4, 3.5, 3, 0.2), mgp = c(2.5, 0.8, 0))
ylabs <- seq(from = 0, to = 1, by = 0.25)
plot(NULL,
  xlim = range(xx), xaxt = "n",
  ylim = c(-0.1, 1), yaxs = "i", yaxt = "n",
  xlab = "", ylab = "Value", bty = "n"
)
axis(side = 2, at = ylabs)
segments(
  x0 = par("usr")[1], x1 = par("usr")[2],
  y0 = 0, y1 = 0, lty = 1, col = "black", lwd = 1
)
barplot(
  metric ~ rf_short + metric_type,
  data = accuracy_ranger_long,
  add = TRUE, axes = FALSE, axisnames = FALSE,
  beside = TRUE, col = rf_formulae_colors,
  legend = TRUE, args.legend = list(
    title = "Binary RF model", horiz = TRUE,
    title.font = 2, bty = "n", xpd = TRUE,
    x = mean(par("usr")[1:2]),
    y = grconvertY(0, from = "ndc", to = "user"),
    xjust = 0.5, yjust = 0
  )
)
for (x in seq_len(nrow(xx))) {
  for (y in seq_len(ncol(xx))) {
    yvalue <- accuracy_ranger_long$metric[
      as.numeric(accuracy_ranger_long$rf) == x &
        as.numeric(accuracy_ranger_long$metric_type) == y
    ]
    text(
      x = xx[x, y], y = yvalue, cex = 0.7, srt = 0,
      labels = round(yvalue, 2), pos = ifelse(yvalue >= 0, 3, 1)
    )
  }
}
text(
  x = apply(xx, 2, mean), y = -0.15, adj = c(0.5, 0.5),
  labels = levels(accuracy_ranger_long$metric_type)
)
dev.off()
# remove objects no longer needed
rm(x, y, yvalue, xx, ylabs)


# Compare 1d accuracy metrics for ordinal RFs -----------------------------

# see metrics of ordfor (ordinal) random forests
names(accuracy_metrics[[2]])

# put together monodimensional metrics
monodim_metrics_ordfor <- c( # names of 1d accuracy metrics
  "overall_accuracy", "macro_F1", "weighted_F1", "MAE",
  "adjacent_error_rate", "nonadjacent_error_rate", "QWK"
)
accuracy_ordfor_1d_wide <- lapply(
  X = accuracy_metrics[to_do$response_type != "binary"],
  FUN = function(x) x[monodim_metrics_ordfor]
) |>
  unlist() |> # unlist
  matrix( # create matrix
    ncol = length(monodim_metrics_ordfor), byrow = TRUE, dimnames = list(
      to_do$rf_formula[to_do$response_type != "binary"], # row names: ordinal rf models
      monodim_metrics_ordfor # column names: 1d accuracy metrics
    )
  )
accuracy_ordfor_1d_wide <- cbind( # bind with data from to_do
  to_do[
    to_do$response_type != "binary",
    c("rf_formula", "rf_short")
  ],
  accuracy_ordfor_1d_wide
)

# reshape into long format
accuracy_ordfor_long <- reshape(
  data = accuracy_ordfor_1d_wide,
  varying = monodim_metrics_ordfor,
  v.names = "metric",
  idvar = "rf_formula",
  times = monodim_metrics_ordfor,
  timevar = "metric_type",
  direction = "long"
)
accuracy_ordfor_long$metric_type <- factor( # metric type as factor
  accuracy_ordfor_long$metric_type,
  levels = monodim_metrics_ordfor
)
levels(accuracy_ordfor_long$metric_type) <- c(
  "Overall\nAccuracy", "Macro F1", "Weighted F1", "MAE",
  "Adjacent\nError Rate", "Non-Adjacent\nError Rate", "QWK"
)
accuracy_ordfor_long$rf_formula <- as.factor(accuracy_ordfor_long$rf_formula)

# define colors of binary rf models
ordfor_rf_colors <- paletteer_d("RColorBrewer::Paired", nlevels(accuracy_ordfor_long$rf_formula))
names(ordfor_rf_colors) <- levels(accuracy_ordfor_long$rf_formula)

# prepare coordinates for barplot
xx <- barplot(
  metric ~ rf_formula + metric_type,
  data = accuracy_ordfor_long,
  beside = TRUE, plot = FALSE
)
rownames(xx) <- levels(accuracy_ordfor_long$rf_formula)
colnames(xx) <- levels(accuracy_ordfor_long$metric_type)
# create barplot
png(file.path(output_folder, "barplot_metrics_ordinal.png"),
  width = 10 * resol, height = 4 * resol, res = resol
)
par(las = 1, xpd = TRUE, tcl = -0.3, mar = c(2, 3.5, 1.5, 0.2), mgp = c(2.5, 0.8, 0))
ylabs <- seq(from = 0, to = max(accuracy_ordfor_long$metric), by = 0.25)
plot(NULL,
  xlim = range(xx), xaxt = "n",
  ylim = range(accuracy_ordfor_long$metric), yaxt = "n",
  xlab = "", ylab = "Value", bty = "n"
)
axis(side = 2, at = ylabs)
segments(
  x0 = par("usr")[1], x1 = par("usr")[2],
  y0 = 0, y1 = 0, lty = 1, col = "black", lwd = 1
)
barplot(
  metric ~ rf_formula + metric_type,
  data = accuracy_ordfor_long,
  add = TRUE, axes = FALSE, axisnames = FALSE,
  beside = TRUE, col = ordfor_rf_colors
)
for (x in seq_len(nrow(xx))) {
  for (y in seq_len(ncol(xx))) {
    yvalue <- accuracy_ordfor_long$metric[
      as.numeric(accuracy_ordfor_long$rf_formula) == x &
        as.numeric(accuracy_ordfor_long$metric_type) == y
    ]
    text(
      x = xx[x, y], y = yvalue, cex = 0.5, srt = 0,
      labels = round(yvalue, 2), pos = ifelse(yvalue >= 0, 3, 1)
    )
  }
}
text(
  x = apply(xx, 2, mean), y = -0.3, adj = c(0.5, 0.5),
  labels = levels(accuracy_ordfor_long$metric_type)
)
legend(
  bty = "n", fill = ordfor_rf_colors,
  x = "topright",
  #y = grconvertY(0, from = "ndc", to = "user"),
  #xjust = 0.5, yjust = 0,
  cex = 0.8, ncol = 4,
  legend = names(ordfor_rf_colors), title = "Ordinal RF model (class ranger)",
  title.font = 2
)
dev.off()
# remove objects no longer needed
rm(x, y, yvalue, xx, ylabs)


# Compare 2d+ accuracy metrics for ordinal RFs ----------------------------

# see metrics of ordfor (ordinal) random forests
names(accuracy_metrics[[2]])

# put together monodimensional metrics
multidim_metrics_ordfor <- c( # names of multidimensional accuracy metrics
  "sensitivity_per_class", "precision_per_class", "f1_per_class"
)
accuracy_ordfor_multi_wide <- lapply(
  X = accuracy_metrics[to_do$response_type != "binary"],
  FUN = function(x) x[multidim_metrics_ordfor]
) |>
  do.call(what = rbind)
accuracy_ordfor_multi_wide <- cbind( # bind with data from to_do
  to_do[
    to_do$response_type != "binary",
    c("rf_formula", "rf_short")
  ],
  accuracy_ordfor_multi_wide
)
accuracy_ordfor_multi_wide

plotmat <- matrix(1:8, ncol = 4)
plotmat <- rbind(9:12, plotmat)
plotmat <- cbind(plotmat, c(0, 13:14))

png(file.path(output_folder, "barplot_metrics_perClass_ordinal.png"),
  height = 6 * resol, width = 10 * resol, res = resol
)
layout(plotmat, heights = c(1, 6, 6))
par(
  las = 2, mar = c(4, 3.5, 0.5, 0.5),
  mgp = c(2.5, 0.8, 0), tcl = -0.3, xpd = TRUE
)
for (i in seq_len(nrow(accuracy_ordfor_multi_wide))) {
  xxtype <- length(accuracy_ordfor_multi_wide[[i, 5]])
  xx <- accuracy_ordfor_multi_wide[i, ][3:5] |>
    unlist() |>
    matrix(nrow = xxtype)
  # xx[is.nan(xx) | is.na(xx)] <- 0
  dimnames(xx) <- list(
    names(accuracy_ordfor_multi_wide[[i, 3]]),
    gsub("_per_class", "", names(accuracy_ordfor_multi_wide)[3:5])
  )
  xxcoords <- barplot(xx,
    las = 2,
    beside = T, ylim = 0:1, space = c(0, 2), plot = FALSE
  )
  plot(NULL,
    xlim = range(xxcoords) + c(-1, 1), ylim = 0:1,
    ylab = "Value", xlab = "", xaxt = "n", bty = "n", yaxs = "i"
  )
  rect(
    xleft = par("usr")[1], xright = par("usr")[2],
    ybottom = par("usr")[3], ytop = par("usr")[4], col = "gray90", border = NA
  )
  barplot(xx,
    add = T,
    beside = T, ylim = c(-0.05, 1.05), las = 2,
    col = if (xxtype == 3) {
      colors_response3
    } else {
      colors_response6
    },
    ylab = "Value", names.arg = rep("", ncol(xx)),
    space = c(0, 2)
  )
  text(
    x = apply(xxcoords, 2, mean), y = -0.05, srt = 45, labels = colnames(xx),
    adj = 1
  )
}
# add column titles (RF formulas)
par(mar = c(0.5, 3.5, 0.5, 0.5), mgp = rep(0, 3))
for (xx in levels(accuracy_ordfor_multi_wide$rf_short)) {
  plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE, ann = FALSE)
  text(
    x = 0, y = 0,
    labels = xx, font = 2, cex = 1.5
  )
}

# add color legend
par(mar = rep(0.5, 4), mgp = rep(0, 3))
plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE, ann = FALSE)
legend(
  x = 0, y = 0, xjust = 0.5, yjust = 0.5,
  legend = names(colors_response3), fill = colors_response3,
  bty = "n", cex = 1.3, title = "Response\n(ordinal, 3 levels)", title.font = 2
)
plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), axes = FALSE, ann = FALSE)
legend(
  x = 0, y = 0, xjust = 0.5, yjust = 0.5,
  legend = names(colors_response6), fill = colors_response6,
  bty = "n", cex = 1.3, title = "Response\n(ordinal, 6 levels)", title.font = 2
)
dev.off()


# Plot confusion matrices -------------------------------------------------
for (i in seq_len(nrow(to_do))) {
  xx <- confusion_matrices[[i]] |> as.data.frame()
  xx <- subset(xx, Observed != "Sum" & Predicted != "Sum") |> droplevels()
  xx <- xx[seq_len(nrow(xx)) |> rep(times = xx$Freq), ]
  xcm <- confusion_matrix(targets = xx$Observed, predictions = xx$Predicted)
  plot_confusion_matrix(xcm,
    add_sums = TRUE, add_normalized = FALSE,
    class_order = rev(levels(xx$Observed)),
    add_zero_shading = F, palette = "Blues",
    intensity_by = "counts", rm_zero_text = FALSE,
    rm_zero_percentages = FALSE,
    sums_settings = sum_tile_settings(
      palette = "Oranges",
      tc_tile_border_color = "black"
    )
  ) |> ggsave(
    filename = file.path(to_do$folder[i], paste0(
      "CM_plot_", to_do$rf_formula[i], ".png"
    )),
    width = 3, height = 3, units = "in", dpi = 300
  )
}
rm(xx, xcm, i)

# Save image --------------------------------------------------------------
save.image("nonsync/07_plots_random_forests.RData")
