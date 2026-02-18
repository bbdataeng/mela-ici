# load("nonsync/06b_plots_rf_metadaset.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
library(caret)
library(ranger)
library(pROC)
library(grid)
library(paletteer)
library(cvms)
library(ggplot2)
library(gt)


# Prepare data ------------------------------------------------------------

load("nonsync/06a_rf_metadataset.RData") # load random forests data

# add RF short name in to_do
to_do$rf_short <- substr(to_do$rf_formula, start = 1, stop = 3) |> factor()

# keep only binary RFs in to_do
to_do <- subset(to_do, response_type == "binary")
rownames(to_do) <- to_do$rf_short


# Plot settings -----------------------------------------------------------

resol <- 300 # resolution in ppi

# prepare colors for RF formulae
rf_formulae_colors <- paletteer_d("RColorBrewer::Set2", nlevels(to_do$rf_short))
names(rf_formulae_colors) <- levels(to_do$rf_short)


# Visualization of RF formulae --------------------------------------------

formulae_data <- cbind(
  LM22 = rep(TRUE, nrow(to_do)),
  to_do[, c("age_and_gender", "hed_data", "technical_predictors")]
) |> as.matrix()
colnames(formulae_data) <- c(
  "LM22", "Age & Sex", "HED", "Technical\nPredictors"
)

# function to add squares in a plot
add_square <- function(x0, y0, ...) {
  x1 <- x0 + 1
  y1 <- y0 + 1
  rect(
    xleft = x0, xright = x1, ybottom = y0, ytop = y1,
    ...
  )
}

# draw heatmap
png(file.path(output_folder, "RF_formulae_heatmap.png"),
  width = 5 * resol, height = 3 * resol, res = resol
)
par(mar = c(0.5, 0.5, 3, 7))
plot(
  NULL,
  xlim = c(0.3, ncol(formulae_data) + 1), xaxs = "i",
  ylim = c(nrow(formulae_data), 0) + 1, yaxs = "i",
  xlab = "", ylab = "",
  bty = "n", axes = FALSE
)
rect(
  xleft = 0.4, xright = 0.9,
  ybottom = seq_len(nrow(formulae_data)), ytop = seq_len(nrow(formulae_data)) + 1,
  col = rf_formulae_colors, xpd = TRUE, border = NA
)
for (x in seq_len(ncol(formulae_data))) {
  for (y in seq_len(nrow(formulae_data))) {
    add_square(x, y, col = ifelse(
      formulae_data[y, x], "#4C72B0", "#D0D0D0"
    ), xpd = TRUE)
  }
}
text(
  x = 0.65, y = seq_len(nrow(formulae_data)) + 0.5,
  adj = 0.5, labels = rownames(formulae_data), xpd = TRUE
)
mtext(
  side = 3, line = 1, at = seq_len(ncol(formulae_data)) + 0.5,
  text = colnames(formulae_data), adj = 0.5, padj = 0.5
)
legend(
  x = mean(c(ncol(formulae_data) + 1, grconvertX(1, from = "ndc", to = "user"))),
  y = mean(par("usr")[3:4]), yjust = 0.5, xjust = 0.5,
  fill = c("#4C72B0", "#D0D0D0"), legend = c("Included", "Not Included"),
  bty = "n", xpd = TRUE
)
dev.off()



# ROC plots with AUC ------------------------------------------------------

# separately by RF model
for (i in seq_len(nrow(to_do))) {
  png(file.path(to_do$folder[i], paste0("ROC_plot_", to_do$rf_formula[i], ".png")),
    width = 3 * resol, height = 3 * resol, res = resol
  )
  plot(rocs_list[[i]],
    print.auc = TRUE,
    auc.polygon = TRUE,
    main = to_do$rf_short[i],
    xaxs = "i", yaxs = "i", las = 1,
    print.auc.x = 0.05, print.auc.y = 0.05, print.auc.cex = 1.3,
    print.auc.adj = c(1, 0),
    auc.polygon.col = rf_formulae_colors[i],
    identity.col = "black", identity.lty = "dashed"
  )
  dev.off()
}

# all together
png(file.path(output_folder, "ROC_plots.png"),
  width = 12 * resol, height = 2 * resol, res = resol
)
par(mfrow = c(1, 6), mar = c(4.1, 3.1, 3.1, 1.1))
for (i in seq_len(nrow(to_do))) {
  plot(rocs_list[[i]],
    print.auc = TRUE,
    auc.polygon = TRUE,
    main = to_do$rf_short[i],
    xaxs = "i", yaxs = "i", las = 1,
    print.auc.x = 0.05, print.auc.y = 0.05, print.auc.cex = 1.3,
    print.auc.adj = c(1, 0),
    auc.polygon.col = rf_formulae_colors[i],
    identity.col = "black", identity.lty = "dashed"
  )
}
dev.off()


# Barplot of median importance score --------------------------------------

# calculate median importances across RFs
median_imp <- apply(vimp_df, 2, median, na.rm = TRUE) |> sort()

# make barplot
png(
  file.path(output_folder, "median_importance_barplot.png"),
  width = 6 * resol, height = 8 * resol, res = resol
)
par(mar = c(3, 15.5, 0.5, 0.5), las = 1, mgp = c(2, 1, 0))
xx <- barplot(median_imp, horiz = TRUE, plot = FALSE)
varnames <- ifelse(names(median_imp) %in% names(cell_types_original),
  cell_types_original[names(median_imp)],
  names(median_imp)
)
plot(NULL,
  xlim = c(0, 100),
  ylim = range(xx), axes = FALSE,
  xlab = "Median Variable Importance", ylab = ""
)
grid()
barplot(median_imp,
  horiz = TRUE, add = TRUE, names.arg = varnames,
  col = "#4C72B0"
)
dev.off()


# Dotplot of importance scores --------------------------------------------

# make dotplot
png(
  file.path(output_folder, "median_importance_dotplot.png"),
  width = 6 * resol, height = 8 * resol, res = resol
)
par(mar = c(3, 15.5, 0.5, 0.5), las = 1, mgp = c(2, 1, 0), xpd = FALSE)
varnames <- ifelse(names(median_imp) %in% names(cell_types_original),
  cell_types_original[names(median_imp)],
  names(median_imp)
)
plot(NULL,
  xlim = c(1, 100),
  ylim = c(0, length(median_imp)) + 0.5, yaxt = "n", yaxs = "i",
  xlab = "Variable Importance", ylab = "",
  bty = "n"
)
grid()
segments(
  x0 = rep(par("usr")[1]), x1 = rep(par("usr")[2]),
  y0 = seq_along(median_imp), y1 = seq_along(median_imp),
  col = "lightgray", lty = "dotted"
)
segments(
  x0 = rep(0, length(median_imp)), x1 = median_imp,
  y0 = seq_along(median_imp), y1 = seq_along(median_imp),
  col = "#4C72B0", lwd = 5, lend = 0
)
set.seed(123456)
for (i in seq_along(median_imp)) {
  points(
    x = vimp_df[, names(median_imp)[i]],
    y = jitter(rep(i, nrow(vimp_df)), amount = 0.15),
    bg = rf_formulae_colors, pch = 21, cex = 1.1
  )
}
axis(side = 2, at = seq_along(varnames), labels = varnames, tick = FALSE)
legend(
  x = "bottomright", pch = 21, pt.bg = rf_formulae_colors,
  title = "Model", title.font = 1,
  legend = names(rf_formulae_colors), ncol = 2, xpd = TRUE
)
dev.off()


# Barplots of importance scores by RF model -------------------------------

# get inverse rank of variable importance scores
xxrank <- lapply(
  X = vimp_list, FUN = function(x) {
    xrank <- rank(x, na.last = TRUE)
    xrank[is.na(x)] <- NA
    xrank <- max(xrank, na.rm = TRUE) + 1 - xrank
    return(xrank)
  }
)

# for each RF
for (i in seq_len(nrow(to_do))) {
  # make barplot of all importances
  png(
    file.path(to_do$folder[i], paste0(
      "vimp_barplot_", to_do$rf_formula[i], ".png"
    )),
    width = 6 * resol, height = 8 * resol, res = resol
  )
  par(mar = c(3, 15.5, 1.5, 0.5), las = 1)
  xx <- barplot(vimp_df[i, ], horiz = TRUE, plot = FALSE)
  xlim <- range(vimp_df[i, ], na.rm = TRUE)
  xlim <- xlim + c(-diff(xlim) * 0.1, diff(xlim) * 0.1)
  varnames <- ifelse(names(vimp_df[i, ]) %in% names(cell_types_original),
    cell_types_original[names(vimp_df[i, ])],
    names(vimp_df[i, ])
  )
  plot(NULL,
    xlim = xlim, ylim = range(xx), axes = FALSE,
    xlab = "Variable Importance", ylab = "",
    main = to_do$rf_short[i]
  )
  barplot(vimp_df[i, ],
    horiz = TRUE, add = TRUE, names.arg = varnames,
    col = rf_formulae_colors[i]
  )
  text(
    x = as.numeric(vimp_df[i, ]), y = xx, labels = xxrank[[i]],
    pos = ifelse(vimp_df[i, ] >= 0 | is.na(vimp_df[i, ]), 4, 2)
  )
  dev.off()
  # make barplot of top 10 importance scores
  png(
    file.path(to_do$folder[i], paste0(
      "vimp_barplot_top10_", to_do$rf_formula[i], ".png"
    )),
    width = 6 * resol, height = 4 * resol, res = resol
  )
  par(mar = c(3, 15.5, 1.5, 0.5), las = 1)
  top10 <- vimp_df[i, ] |>
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
    main = to_do$rf_short[i]
  )
  barplot(top10,
    horiz = TRUE, add = TRUE,
    names.arg = varnames,
    col = rf_formulae_colors[i]
  )
  text(x = as.numeric(top10), y = xx, labels = 1:10, pos = 4)
  dev.off()
}


# Compare accuracy metrics for binary RFs ---------------------------------

# put metrics together
metrics_binary <- lapply(
  X = names(rf_cv_binary_list), FUN = function(x) {
    xauc <- rocs_list[[x]]$auc |> as.numeric()
    x1 <- confusion_matrices[[x]]$overall[c("Accuracy", "Kappa")]
    x2 <- confusion_matrices[[x]]$byClass
    return(c(
      AUC = xauc,
      x1, x2
    ))
  }
)
names(metrics_binary) <- names(rf_cv_binary_list)

# prepare a dataframe
accuracy_ranger_wide <- cbind(
  to_do[
    to_do$response_type == "binary",
    c("rf_short", "technical_predictors", "age_and_gender")
  ],
  sapply(
    X = to_do$rf_formula[to_do$response_type == "binary"],
    FUN = function(x) metrics_binary[[x]]
  ) |>
    t() |> as.data.frame()
)

# reshape into long format
accuracy_ranger_long <- reshape(
  data = accuracy_ranger_wide,
  varying = c(
    "AUC", "Accuracy", "Kappa", "Sensitivity", "Specificity",
    "Pos Pred Value", "Neg Pred Value", "Precision", "Recall",
    "F1", "Prevalence", "Detection Rate", "Detection Prevalence",
    "Balanced Accuracy"
  ),
  v.names = "metric",
  idvar = "rf_short",
  times = c(
    "AUC", "Accuracy", "Kappa", "Sensitivity", "Specificity",
    "Pos Pred Value", "Neg Pred Value", "Precision", "Recall",
    "F1", "Prevalence", "Detection Rate", "Detection Prevalence",
    "Balanced Accuracy"
  ),
  timevar = "metric_type",
  direction = "long"
)
accuracy_ranger_long$metric_type <- factor( # metric type as factor
  accuracy_ranger_long$metric_type,
  levels = c(
    "AUC", "Accuracy", "Kappa", "Sensitivity", "Specificity",
    "Pos Pred Value", "Neg Pred Value", "Precision", "Recall",
    "F1", "Prevalence", "Detection Rate", "Detection Prevalence",
    "Balanced Accuracy"
  )
)
accuracy_ranger_long <- subset(
  accuracy_ranger_long, metric_type %in% c(
    "Accuracy", "Sensitivity", "Specificity",
    "Pos Pred Value", "Neg Pred Value"
  )
) |> droplevels()
levels(accuracy_ranger_long$metric_type) <- c(
  "Accuracy", "Sensitivity", "Specificity", "PPV", "NPV"
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
  width = 6 * resol, height = 3.5 * resol, res = resol
)
par(las = 1, xpd = TRUE, tcl = -0.3, mar = c(2, 3.5, 4, 0.2), mgp = c(2.5, 0.8, 0))
ylabs <- seq(from = 0, to = 1, by = 0.25)
plot(NULL,
  xlim = range(xx), xaxt = "n",
  ylim = c(0, 1), yaxs = "i", yaxt = "n",
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
    title = "Random Forest model", horiz = TRUE,
    title.font = 2, bty = "n", xpd = TRUE,
    x = mean(par("usr")[1:2]),
    y = grconvertY(1, from = "ndc", to = "user"),
    xjust = 0.5, yjust = 1
  )
)
for (x in seq_len(nrow(xx))) {
  for (y in seq_len(ncol(xx))) {
    yvalue <- fakeyvalue <- accuracy_ranger_long$metric[
      as.numeric(accuracy_ranger_long$rf) == x &
        as.numeric(accuracy_ranger_long$metric_type) == y
    ]
    if (fakeyvalue < 0) fakeyvalue <- 0
    fakeyvalue <- fakeyvalue + 0.03
    text(
      labels = round(yvalue, 3), cex = 0.7,
      x = xx[x, y], y = fakeyvalue, srt = 90, adj = 0
    )
  }
}
mtext(
  text = levels(accuracy_ranger_long$metric_type),
  side = 1, at = apply(xx, 2, mean),
  line = 0.5
)

dev.off()
# remove objects no longer needed
rm(x, y, yvalue, xx, ylabs)


# Plot confusion matrices -------------------------------------------------
for (i in seq_len(nrow(to_do))) {
  xform <- to_do$rf_formula[i]
  xx <- confusion_matrices[[xform]]$table |> as.data.frame()
  xx <- xx[seq_len(nrow(xx)) |> rep(times = xx$Freq), ]
  xcm <- confusion_matrix(targets = xx$Reference, predictions = xx$Prediction)
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
rm(xx, xcm, i, xform)


# Table of sample sizes ---------------------------------------------------

tab2 <- lapply(to_do$rf_formula, FUN = function(x) {
  c(
    n_train = nrow(df_train_list[[x]]),
    tab_train = table(df_train_list[[x]]$response),
    n_test = nrow(df_test_list[[x]]),
    tab_test = table(df_test_list[[x]]$response)
  )
}) |>
  do.call(what = rbind) |>
  as.data.frame()
tab2 <- cbind(RF = to_do$rf_short, tab2)

# save table as excel
writexl::write_xlsx(
  tab2,
  file.path(output_folder, "sample_sizes_rawtable.xlsx"),
  format_headers = TRUE
)


# nicely format table and save as figure
xtab <- gt(tab2) |>
  tab_spanner(
    label = md("**Training/Validation data**"),
    columns = contains("train")
  ) |>
  tab_spanner(
    label = md("**Test data**"),
    columns = contains("test")
  ) |>
  cols_label(
    ends_with(".R") ~ "R",
    ends_with(".NR") ~ "NR",
    starts_with("n_") ~ "Total",
    matches("RF") ~ md("**Model**")
  ) |>
  cols_width(everything() ~ px(70)) |>
  cols_align("center", everything()) |>
  data_color(columns = "RF", method = "factor", palette = rf_formulae_colors)
gtsave(xtab, file.path(output_folder, "sample_sizes_table.png"))

# Save image --------------------------------------------------------------
save.image("nonsync/06b_plots_rf_metadaset.RData")
