# load("nonsync/07_rf_transcriptomics.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
library(ranger) # v0.17.0
library(ordinalForest) # v2.4-4
library(pROC) # v1.19.0.1
library(grid)
library(cvms)
library(paletteer)
library(ggplot2)
library(fields)
library(circlize)
source("rf_functions.R")


# Load and prepare data ---------------------------------------------------

# define output folder
output_folder <- "nonsync/07_rf_transcriptomics"
if (!dir.exists(output_folder)) dir.create(output_folder)

# load metadata
allmetadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds") |> as.data.frame()

# load VST matrix
vst <- read.delim("./nonsync/autogo_response_pre_icb/deseq_vst_data.txt",
  header = TRUE, row.names = 1
)

# inspect samples for which data are available
table(allmetadata$accession %in% colnames(vst)) |>
  addmargins() # 27/165 samples are not present in the VST data

# define samples
samples <- intersect(allmetadata$accession, colnames(vst))

# subset data
metadata <- subset(allmetadata, accession %in% samples)
vst <- vst[, metadata$accession]
stopifnot(identical(colnames(vst), metadata$accession)) # check that data match

# put data together as a data frame
xdata <- data.frame(
  response = metadata$response_2levels, # response (binary)
  t(vst) # vst data (trasponsed: genes in columns)
)
xdata$response <- factor(xdata$response, ordered = FALSE) # response as an unordered factor
table(xdata$response)

# check missing data
table(complete.cases(xdata)) # no missing data


# Split dataset in train and test subsets ---------------------------------

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

# get training data (~80%) and test data (~20%)
set.seed(123)
train_idx <- xdata$response |> stratified_split(p_train = 0.8)
train_data <- xdata[train_idx, ]
test_data <- xdata[-train_idx, ]

# save data
write.table(x = train_data, file = file.path(
  output_folder, "data_train.txt"
), sep = "\t")
write.table(x = test_data, file = file.path(
  output_folder, "data_test.txt"
), sep = "\t")


# save sample sizes as txt
n_train <- nrow(train_data)
n_test <- nrow(test_data)
sink(file.path(output_folder, "sampleSizes.txt"))
cat("===== Training data =====", sep = "\n")
cat(paste0(
  n_train, " observations (", round(n_train / (n_train + n_test), 3) * 100,
  "%)"
), sep = "\n")
cat("\nResponse levels:")
print(train_data$response |> table())
cat("\n===== Test data =====", sep = "\n")
cat(paste0(
  n_test, " observations (", round(n_test / (n_train + n_test), 3) * 100,
  "%)"
), sep = "\n")
cat("\nResponse levels:")
print(test_data$response |> table())
sink()




# Train random forest -----------------------------------------------------

# define formula
xformula <- paste0("response ~ '", paste(
  rownames(vst),
  collapse = "' + '"
), "'") |> as.formula()

# random forest
xrf <- ranger(
  y = train_data$response,
  x = train_data[, setdiff(colnames(train_data), "response")],
  importance = "permutation",
  probability = TRUE,
  min.node.size = 5,
  na.action = "na.fail", # no NAs expected in the training data
  respect.unordered.factors = "partition",
  num.trees = 1000,
  mtry = floor(sqrt(nrow(vst))),
  seed = 1 # set seed for reproducibility
)

# save rf object
saveRDS(xrf, file.path(output_folder, "rfObject.rds"))


# Examine prediction performance ------------------------------------------

# get confusion matrix
conf_matrix <- get_confusion_matrix(
  rf_object = xrf,
  testdata = test_data,
  show_sum = TRUE
)

# prediction performance metrics
ppm <- get_accuracy_metrics(
  rf_object = xrf,
  testdata = test_data,
  confusion_matrix = conf_matrix,
  positive_level = "R"
)

# save confusion matrices and accuracy metrics
sink(file.path(output_folder, "CM_accuracyMetrics.txt"))
cat("===== Confusion Matrix =====", sep = "\n")
print(conf_matrix)
cat("\n===== Accuracy Metrics =====", sep = "\n")
print(ppm)
sink()


# Get variable importance -------------------------------------------------

# get ordered importance scores
importances <- get_importance(rf_object = xrf)
any(is.na(importances)) # no NAs present
importances_rank <- seq_along(importances)
names(importances_rank) <- names(importances)

# save importance scores
# prepare ordered dataframe
xx <- data.frame(
  variable = names(importances),
  importance = importances,
  rank = seq_along(importances)
)
# save dataframe
write.table(xx, file.path(output_folder, "variable_importance.txt"),
  row.names = FALSE, sep = "\t"
)
rm(xx)


# Plot settings -----------------------------------------------------------

resol <- 300 # resolution in ppi

# prepare colors for response (binary)
colors_response <- paletteer_c("grDevices::RdYlBu", nlevels(xdata$response))
names(colors_response) <- levels(xdata$response)


# Barplot of top importance scores ----------------------------------------

ntop <- 50
plotdata <- head(importances, ntop) |> rev()
png(file.path(output_folder, paste0(
  "vimp_barplot_top", ntop, ".png"
)), width = 5 * resol, height = 10 * resol, res = resol)
par(mar = c(5, 8, 0, 0.5), las = 2)
barplot(plotdata,
  horiz = TRUE,
  col = "#1F78B4", border = NA
)
mtext(side = 1, text = "Importance", line = par("mar")[1] - 1, las = 0)
mtext(side = 2, text = paste0("Genes (ranked, top ", ntop, ")"), line = par("mar")[2] - 1, las = 0)
dev.off()


# Importance scores vs. LogFC ---------------------------------------------

# prepare data
plotdata <- data.frame(
  genes = names(importances),
  importance = importances,
  importance_rank = seq_along(importances)
) |> head(50)
# remove gene subfix (version)
plotdata$genes <- gsub("\\.[A-Z0-9]*", "", plotdata$genes)
# load and merge LogFC data
logfc <- read.delim("nonsync/autogo_response_pre_icb/NR_vs_R/DE_NR_vs_R_allres.tsv")
plotdata <- merge(plotdata, logfc, all.x = TRUE, sort = FALSE)
rownames(plotdata) <- plotdata$genes
plotdata <- plotdata[order(plotdata$importance_rank, decreasing = TRUE), ]
# prepare colors
xpal <- paletteer_c("grDevices::RdBu", 251, -1)
col_breaks <- seq(
  from = -max(abs(plotdata$log2FoldChange), na.rm = TRUE),
  to = max(abs(plotdata$log2FoldChange), na.rm = TRUE),
  length.out = length(xpal)
)
col_fun <- colorRamp2(breaks = col_breaks, colors = xpal)
# barplot
png(file.path(output_folder, paste0(
  "vimp_barplot_LogFC_top", ntop, ".png"
)), width = 5 * resol, height = 10 * resol, res = resol)
par(mar = c(5, 8, 0, 0.5), las = 2)
barplot(importances |> head(ntop) |> rev(),
  horiz = TRUE,
  col = ifelse(
    is.na(plotdata$log2FoldChange), "grey40", col_fun(plotdata$log2FoldChange)
  )
)
mtext(side = 1, text = "Importance", line = par("mar")[1] - 1, las = 0)
mtext(side = 2, text = paste0("Genes (ranked, top ", ntop, ")"), line = par("mar")[2] - 1, las = 0)
image.plot(
  legend.only = TRUE,
  horizontal = FALSE,
  zlim = range(col_breaks),
  col = col_fun(col_breaks),
  legend.mar = 4,
  legend.shrink = 0.7,
  smallplot = c(0.8, 0.85, 0.15, 0.35),
  legend.args = list(text = "Log2FC", side = 3, line = 0.5, las = 0)
)
dev.off()



# Barplot of prediction performance metrics -------------------------------

plotdata <- data.frame(
  metric_type = factor(names(ppm), levels = c(
    "AUC", "accuracy", "sensitivity", "specificity", "precision",
    "f1_score", "balanced_accuracy", "mcc"
  )),
  value = ppm,
  row.names = NULL
)
levels(plotdata$metric_type) <- c(
  "AUC", "Accuracy", "Sensitivity", "Specificity", "Precision",
  "F1 score", "Balanced\nAccuracy", "MCC"
)

# prepare coordinates for barplot
xx <- barplot(
  value ~ metric_type,
  data = plotdata, plot = FALSE
)
rownames(xx) <- levels(plotdata$metric_type)

# create barplot
png(file.path(output_folder, "barplot_metrics.png"),
  width = 5 * resol, height = 3.5 * resol, res = resol
)
par(las = 2, xpd = TRUE, tcl = -0.3, mar = c(5, 3.5, 0.8, 0.2), mgp = c(2.5, 0.8, 0))
xx <- barplot(
  value ~ metric_type,
  data = plotdata, ylim = 0:1,
  add = FALSE, axes = TRUE, axisnames = TRUE,
  col = "#1F78B4", ylab = "Value", xlab = "", border = NA
)
for (x in seq_len(nrow(plotdata))) {
  yvalue <- fakeyvalue <- plotdata$value[x]
  if (fakeyvalue < 0) fakeyvalue <- 0
  fakeyvalue <- fakeyvalue + 0.05
  text(
    labels = round(yvalue, 2), cex = 1,
    x = xx[x, ], y = fakeyvalue, srt = 90, adj = 0
  )
}
dev.off()
# remove objects no longer needed
rm(x, yvalue, xx)


# Plot confusion matrix ---------------------------------------------------

xx <- conf_matrix |> as.data.frame()
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
  filename = file.path(output_folder, "CM_plot.png"),
  width = 3, height = 3, units = "in", dpi = 300
)

# Save image --------------------------------------------------------------

save.image("nonsync/07_rf_transcriptomics.RData")
