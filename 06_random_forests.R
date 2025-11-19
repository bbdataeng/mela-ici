# load("nonsync/06_random_forests.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
library(ranger) # v0.17.0
library(ordinalForest) # v2.4-4
library(ComplexHeatmap) # v2.24.1
library(grid)
library(paletteer)
source("rf_functions.R")


# Prepare new directories -------------------------------------------------

# main output folder
output_folder <- "nonsync/06_random_forests"
if (!dir.exists(output_folder)) dir.create(output_folder)

# define combinations of rf models to be created
to_do <- expand.grid(
  response_type = c("binary", "ordinal3", "ordinal6"),
  technical_predictors = c(TRUE, FALSE),
  age_and_gender = c(TRUE, FALSE)
)

# create new folders
to_do$folder <- character(nrow(to_do))
for (i in seq_len(nrow(to_do))) {
  to_do$folder[i] <- file.path(
    output_folder, # main output folder
    paste0(to_do$response_type[i], "_response"), # binary or ordinal response
    ifelse(to_do$technical_predictors[i], # including or technical predictors or not
      "with_technical_preds", "wo_technical_preds"
    ),
    ifelse(to_do$age_and_gender[i], # including age and gender or not
      "with_age_gender", "wo_age_gender"
    )
  )
  if (!dir.exists(to_do$folder[i])) dir.create(to_do$folder[i], recursive = TRUE)
}

# plot resolution (pixel per inch)
resol <- 300


# Load and prepare data ---------------------------------------------------

metadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds") |> as.data.frame()
xdata <- readRDS("nonsync/01_clean_data/clean_cibersortx.rds") |> as.data.frame()
hed_data <- readRDS("nonsync/01_clean_data/clean_hed.rds") |> as.data.frame()

# which columns have NAs?
metadata |>
  apply(2, function(x) any(is.na(x))) |>
  which() # age, gender, enrichment_protocol
xdata |>
  apply(2, function(x) any(is.na(x))) |>
  which() # no missing values
hed_data |>
  apply(2, function(x) any(is.na(x))) |>
  which() # no missing values

# keep only necessary metadata columns
metadata <- metadata[, c(
  "accession", "response_6levels", "response_3levels", "response_2levels",
  "age", "gender", "enrichment_protocol", "dataset"
)]

# remove mean HED column
hed_data <- hed_data[, setdiff(names(hed_data), "HED_mean")]

# sanitize cell-type column names
cell_types_original <- colnames(xdata)
cell_types <- gsub(" ", "_", cell_types_original)
cell_types <- gsub("\\(|\\)", "", cell_types)
colnames(xdata) <- cell_types
names(cell_types_original) <- cell_types

# merge
alldata <- cbind(metadata, hed_data, xdata)
rm(metadata, hed_data, xdata)

# check response levels
str(alldata$response_6levels) # ordered factor (6 levels)
str(alldata$response_3levels) # ordered factor (3 levels)
str(alldata$response_2levels) # ordered factor (2 levels)
alldata$response_2levels <- factor( # response_2levels non-ordered
  alldata$response_2levels,
  ordered = FALSE
)

# define new dataframe for each row of to_do
to_do$df_all <- paste(
  "df", "all", gsub("inal", "", to_do$response_type),
  as.numeric(to_do$technical_predictors),
  as.numeric(to_do$age_and_gender),
  sep = "_"
)
for (i in seq_len(nrow(to_do))) {
  vars_to_include <- names(alldata)
  # exclude columns
  if (to_do$response_type[i] == "binary") { # select binary response
    vars_to_include <- setdiff(vars_to_include, c("response_6levels", "response_3levels"))
  } else if (to_do$response_type[i] == "ordinal3") { # select 3-level ordinal response
    vars_to_include <- setdiff(vars_to_include, c("response_6levels", "response_2levels"))
  } else { # select 6-level ordinal response
    vars_to_include <- setdiff(vars_to_include, c("response_3levels", "response_2levels"))
  }
  # exclude technical predictors when appropriate
  if (!to_do$technical_predictors[i]) {
    vars_to_include <- setdiff(
      vars_to_include, c("enrichment_protocol", "dataset")
    )
  }
  # exclude age and gender when appropriate
  if (!to_do$age_and_gender[i]) {
    vars_to_include <- setdiff(
      vars_to_include, c("age", "gender")
    )
  }
  # exclude incomplete cases
  xx <- alldata[, vars_to_include] |>
    na.exclude() |>
    droplevels()
  # uniform response variable name
  names(xx) <- gsub("response_.levels", "response", names(xx))
  # assign new object
  assign(x = to_do$df_all[i], value = xx)
}
rm(xx, vars_to_include)
ls(pattern = "^df") # new dataframes created



# Split dataset in train and test subsets ---------------------------------

# prepare new object names
to_do$df_train <- gsub("all", "train", to_do$df_all) # training data
to_do$df_test <- gsub("all", "test", to_do$df_all) # test data

# function for split stratified over the levels of a variable y
stratified_split <- function(y, p_train = 0.8) {
  stopifnot(is.factor(y))
  idx_train <- integer(0)
  for (lvl in levels(y)) {
    idx <- which(y == lvl)
    if (length(idx) == 0) next
    n_tr <- max(1, floor(length(idx) * p_train))
    idx_train <- c(idx_train, sample(idx, n_tr))
  }
  sort(unique(idx_train))
}

# get training data (~80%) and test data (~20%) for each case
set.seed(123)
for (i in seq_len(nrow(to_do))) {
  train_idx <- get(to_do$df_all[i])$response |> stratified_split(p_train = 0.8)
  assign(x = to_do$df_train[i], value = get(to_do$df_all[i])[train_idx, ])
  assign(x = to_do$df_test[i], value = get(to_do$df_all[i])[-train_idx, ])
}
ls(pattern = "^df") # new dataframes created

# save data in each folder
for (i in seq_len(nrow(to_do))) {
  write.csv(x = get(to_do$df_train[i]), file = file.path(
    to_do$folder[i], "data_train.csv"
  ), row.names = FALSE)
  write.csv(x = get(to_do$df_test[i]), file = file.path(
    to_do$folder[i], "data_test.csv"
  ), row.names = FALSE)
}

# define formulas
for (i in seq_len(nrow(to_do))) {
  predictors <- get(to_do$df_all[i]) |>
    names() |>
    setdiff(y = c("accession", "response")) |>
    paste(collapse = " + ")
  to_do$formula[i] <- paste0("response ~ ", predictors)
}

# save formula and number of observations in test/train data in each folder
for (i in seq_len(nrow(to_do))) {
  n_train <- nrow(get(to_do$df_train[i]))
  n_test <- nrow(get(to_do$df_test[i]))
  sink(file.path(to_do$folder[i], "formula_sampleSizes.txt"))
  cat("===== Formula =====", sep = "\n")
  cat(to_do$formula[i], "\n")
  cat("\n===== Excluded data (NAs present) =====", sep = "\n")
  cat(paste0(nrow(alldata) - n_train - n_test, " observations"), sep = "\n")
  cat("\n===== Training data =====", sep = "\n")
  cat(paste0(
    n_train, " observations (", round(n_train / (n_train + n_test), 3) * 100,
    "%)"
  ), sep = "\n")
  cat("\nResponse levels:")
  print(get(to_do$df_train[i])$response |> table())
  cat("\n===== Test data =====", sep = "\n")
  cat(paste0(
    n_test, " observations (", round(n_test / (n_train + n_test), 3) * 100,
    "%)"
  ), sep = "\n")
  cat("\nResponse levels:")
  print(get(to_do$df_test[i])$response |> table())
  sink()
}


# Create random forests ---------------------------------------------------

# prepare names of rf objects
to_do$rf <- gsub("df_all", "rf", to_do$df_all)

# get random forests on binary response
for (i in which(to_do$response_type == "binary")) {
  # get number of predictors
  n_predictors <- strsplit(to_do$formula[i], split = " \\+ ") |>
    unlist() |>
    length()
  # random forest
  ranger(
    formula = as.formula(to_do$formula[i]), # rf formula
    data = get(to_do$df_train[i]), # training data
    importance = "permutation",
    probability = TRUE,
    min.node.size = 5,
    na.action = "na.fail", # no NAs expected in the training data
    respect.unordered.factors = "partition",
    num.trees = 1000,
    mtry = floor(sqrt(n_predictors)),
    seed = 1 # set seed for reproducibility
  ) |> assign(x = to_do$rf[i]) # assign new object
}
ls(pattern = "^rf_bin") # new random forests objects

# get random forests on ordinal response
for (i in which(to_do$response_type != "binary")) {
  # get predictors
  predictors <- to_do$formula[i] |>
    gsub(pattern = "response ~ ", replacement = "") |>
    strsplit(split = " \\+ ") |>
    unlist()
  # set seed for reproducibility
  set.seed(1)
  # ordinal random forest
  ordfor(
    depvar = "response", # response variable name
    data = get(to_do$df_train[i])[, c("response", predictors)],
    nsets = 500,
    ntreeperdiv = 100,
    ntreefinal = 3000,
    mtry = floor(sqrt(length(predictors))),
    min.node.size = 10,
    perffunction = "probability", # uses Ranked Probability Score (RPS)
    importance = "rps",
    num.threads = 0 # 0 = auto
  ) |> assign(x = to_do$rf[i]) # assign new object
}
ls(pattern = "^rf_ord") # new random forests objects

# save rf objects
for (i in seq_len(nrow(to_do))) {
  saveRDS(
    object = get(to_do$rf[i]),
    file = file.path(to_do$folder[i], "rf.rds")
  )
}


# Get confusion matrices --------------------------------------------------
confusion_matrices <- lapply( # list of confusion matrices
  X = seq_len(nrow(to_do)),
  FUN = function(i) {
    get_confusion_matrix(
      rf_object = get(to_do$rf[i]),
      testdata = get(to_do$df_test[i]),
      show_sum = TRUE
    )
  }
)
names(confusion_matrices) <- to_do$rf

# Get accuracy metrics ----------------------------------------------------
accuracy_metrics <- lapply( # list of accuracy metrics
  X = seq_len(nrow(to_do)),
  FUN = function(i) {
    get_accuracy_metrics(
      rf_object = get(to_do$rf[i]),
      confusion_matrix = confusion_matrices[[i]],
      positive_level = "R"
    )
  }
)
names(accuracy_metrics) <- to_do$rf

# save confusion matrices and accuracy metrics in the respective folders
for (i in seq_len(nrow(to_do))) {
  sink(file.path(to_do$folder[i], "confusionMatrix_accuracyMetrics.txt"))
  cat("===== Confusion Matrix =====", sep = "\n")
  print(confusion_matrices[[i]])
  cat("\n===== Accuracy Metrics =====", sep = "\n")
  print(accuracy_metrics[[i]])
  sink()
}


# Get variable importance -------------------------------------------------
importances <- lapply( # list of variable importances
  X = seq_len(nrow(to_do)),
  FUN = function(i) get_importance(rf_object = get(to_do$rf[i]))
)
names(importances) <- to_do$rf

# create dataframe of variable importance scores
vars <- names(alldata) |> grepv(pattern = "response|accession", invert = TRUE)
for (i in seq_along(importances)) {
  importances[[i]] <- importances[[i]][vars]
  names(importances[[i]]) <- vars
}
importances_df <- do.call(rbind, importances)

# get ranked variable importance scores
importances_rank <- lapply(
  X = importances, FUN = function(x) {
    xrank <- rank(x, na.last = TRUE)
    xrank[is.na(x)] <- NA
    xrank <- max(xrank, na.rm = TRUE) + 1 - xrank
    return(xrank)
  }
)

# save data and barplot of importance scores
for (i in seq_len(nrow(to_do))) {
  # prepare ordered dataframe
  xx <- data.frame(
    variable = names(importances[[i]]),
    rank = rank(importances[[i]], na.last = TRUE),
    importance = importances[[i]]
  )
  xx$rank[is.na(xx$importance)] <- NA
  xx$rank <- max(xx$rank, na.rm = TRUE) + 1 - xx$rank
  xx <- xx[order(xx$rank, na.last = TRUE), ]
  rownames(xx) <- NULL
  # save dataframe
  write.table(xx, file.path(to_do$folder[i], "variable_importance.txt"),
    row.names = FALSE, sep = "\t"
  )
  # make barplot of all importance scores
  png(file.path(to_do$folder[i], "variable_importance_barplot.png"),
    width = 8 * resol, height = 8 * resol, res = resol
  )
  par(mar = c(5.1, 12.5, 0.5, 0.5), las = 2)
  xx <- barplot(importances_df[i, ], horiz = TRUE, plot = FALSE)
  xlim <- range(importances_df[i, ], na.rm = TRUE)
  xlim <- xlim + c(-diff(xlim) * 0.1, diff(xlim) * 0.1)
  varnames <- ifelse(names(importances_df[i, ]) %in% names(cell_types_original),
    cell_types_original[names(importances_df[i, ])],
    names(importances_df[i, ])
  )
  plot(NULL,
    xlim = xlim, ylim = range(xx), axes = FALSE,
    xlab = "Variable Importance", ylab = ""
  )
  barplot(importances_df[i, ],
    horiz = TRUE, add = TRUE, names.arg = varnames
  )
  text(
    x = as.numeric(importances_df[i, ]), y = xx, labels = importances_rank[[i]],
    pos = ifelse(importances_df[i, ] >= 0 | is.na(importances_df[i, ]), 4, 2)
  )
  dev.off()
  # make barplot of top 10 importance scores
  png(file.path(to_do$folder[i], "variable_importance_barplot_top10.png"),
    width = 6 * resol, height = 4 * resol, res = resol
  )
  par(mar = c(5.1, 12.5, 0.5, 0.5), las = 2)
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
    xlab = "Variable Importance", ylab = ""
  )
  barplot(top10,
    horiz = TRUE, add = TRUE,
    names.arg = varnames
  )
  text(x = as.numeric(top10), y = xx, labels = 1:10, pos = 4)
  dev.off()
}
rm(xlim, xx, top10, varnames, i)


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
  to_do[, c("response_type", "technical_predictors", "age_and_gender")],
  row.names = colnames(vip_rank)
)
ann_df$response_levels <- ifelse(
  ann_df$response_type == "binary", 2, ifelse(
    ann_df$response_type == "ordinal3", 3, 6
  )
)
ann_df$rf_class <- ifelse(ann_df$response_levels == 2, "ranger", "ordfor")
ann_df <- ann_df[, c("rf_class", "response_levels", "technical_predictors", "age_and_gender")]

# prepare annotation colors
palettes <- c(
  "ggsci::default_nejm", "ggthemes::Green",
  "ggsci::default_jama", "ggsci::default_jama"
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
  cluster_rows = TRUE,
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
    c("rf", "technical_predictors", "age_and_gender")
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
  idvar = "rf",
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
accuracy_ranger_long$rf <- as.factor(accuracy_ranger_long$rf)

# define colors of binary rf models
ranger_rf_colors <- paletteer_d("RColorBrewer::Paired", 2 * nlevels(accuracy_ranger_long$rf))
ranger_rf_colors <- ranger_rf_colors[!as.logical(seq_along(ranger_rf_colors) %% 2)]
names(ranger_rf_colors) <- levels(accuracy_ranger_long$rf)

# prepare coordinates for barplot
xx <- barplot(
  metric ~ rf + metric_type,
  data = accuracy_ranger_long,
  beside = TRUE, plot = FALSE
)
rownames(xx) <- levels(accuracy_ranger_long$rf)
colnames(xx) <- levels(accuracy_ranger_long$metric_type)
# create barplot
png(file.path(output_folder, "barplot_metrics_binary.png"),
  width = 10 * resol, height = 4 * resol, res = resol
)
par(las = 1, xpd = TRUE, tcl = -0.3, mar = c(6, 3.5, 0.5, 0.2), mgp = c(2.5, 0.8, 0))
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
  metric ~ rf + metric_type,
  data = accuracy_ranger_long,
  add = TRUE, axes = FALSE, axisnames = FALSE,
  beside = TRUE, col = ranger_rf_colors
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
  x = apply(xx, 2, mean), y = -0.2, adj = c(0.5, 0.5),
  labels = levels(accuracy_ranger_long$metric_type)
)
legend(
  bty = "n", horiz = TRUE, fill = ranger_rf_colors,
  x = mean(par("usr")[1:2]),
  y = grconvertY(0, from = "ndc", to = "user"),
  xjust = 0.5, yjust = 0, cex = 0.8,
  legend = c(
    "wo technical predictors,\nage and gender",
    "wo technical predictors",
    "wo age and gender",
    "all predictors"
  ), title = "Binary RF model (class ranger)", title.font = 2
)
dev.off()
# remove objects no longer needed
rm(x, y, yvalue, xx, ylabs)


# Compare accuracy metrics for ordinal RFs --------------------------------

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
      to_do$rf[to_do$response_type != "binary"], # row names: ordinal rf models
      monodim_metrics_ordfor # column names: 1d accuracy metrics
    )
  )
accuracy_ordfor_1d_wide <- cbind( # bind with data from to_do
  to_do[
    to_do$response_type != "binary",
    c("rf", "technical_predictors", "age_and_gender")
  ],
  accuracy_ordfor_1d_wide
)

# reshape into long format
accuracy_ordfor_long <- reshape(
  data = accuracy_ordfor_1d_wide,
  varying = monodim_metrics_ordfor,
  v.names = "metric",
  idvar = "rf",
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
accuracy_ordfor_long$rf <- factor(accuracy_ordfor_long$rf, levels = c(
  "rf_ord3_0_0", "rf_ord6_0_0",
  "rf_ord3_0_1", "rf_ord6_0_1",
  "rf_ord3_1_0", "rf_ord6_1_0",
  "rf_ord3_1_1", "rf_ord6_1_1"
))

# define colors of binary rf models
ordfor_rf_colors <- paletteer_d("RColorBrewer::Paired", nlevels(accuracy_ordfor_long$rf))
names(ordfor_rf_colors) <- levels(accuracy_ordfor_long$rf)

# prepare coordinates for barplot
xx <- barplot(
  metric ~ rf + metric_type,
  data = accuracy_ordfor_long,
  beside = TRUE, plot = FALSE
)
rownames(xx) <- levels(accuracy_ordfor_long$rf)
colnames(xx) <- levels(accuracy_ordfor_long$metric_type)
# create barplot
png(file.path(output_folder, "barplot_metrics_ordinal.png"),
  width = 10 * resol, height = 4 * resol, res = resol
)
par(las = 1, xpd = TRUE, tcl = -0.3, mar = c(6, 3.5, 0.5, 0.2), mgp = c(2.5, 0.8, 0))
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
  metric ~ rf + metric_type,
  data = accuracy_ordfor_long,
  add = TRUE, axes = FALSE, axisnames = FALSE,
  beside = TRUE, col = ordfor_rf_colors
)
for (x in seq_len(nrow(xx))) {
  for (y in seq_len(ncol(xx))) {
    yvalue <- accuracy_ordfor_long$metric[
      as.numeric(accuracy_ordfor_long$rf) == x &
        as.numeric(accuracy_ordfor_long$metric_type) == y
    ]
    text(
      x = xx[x, y], y = yvalue, cex = 0.5, srt = 0,
      labels = round(yvalue, 2), pos = ifelse(yvalue >= 0, 3, 1)
    )
  }
}
text(
  x = apply(xx, 2, mean), y = -0.4, adj = c(0.5, 0.5),
  labels = levels(accuracy_ordfor_long$metric_type)
)
legend(
  bty = "n", fill = ordfor_rf_colors,
  x = mean(par("usr")[1:2]),
  y = grconvertY(0, from = "ndc", to = "user"),
  xjust = 0.5, yjust = 0, cex = 0.8, ncol = 4,
  legend = c(
    "3 levels - wo technical, age and gender",
    "6 levels - wo technical, age and gender",
    "3 levels - wo technical",
    "6 levels - wo technical",
    "3 levels - wo age and gender",
    "6 levels - wo age and gender",
    "3 levels - all predictors",
    "6 levels - all predictors"
  ), title = "Ordinal RF model (class ranger)", title.font = 2
)
dev.off()
# remove objects no longer needed
rm(x, y, yvalue, xx, ylabs)


# Save image --------------------------------------------------------------

save.image("nonsync/06_random_forests.RData")
