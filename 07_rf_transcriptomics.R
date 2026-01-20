# load("nonsync/07_rf_transcriptomics.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(ranger) # v0.17.0
  library(caret) # v7.0-1
  library(MLmetrics) # v1.1.3
  library(pROC) # v1.19.0.1
  library(grid) # v4.5.2
  library(cvms) # v2.0.0
  library(paletteer)
  library(ggplot2)
  library(fields)
  library(circlize)
  source("rf_functions.R")
})


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

# get training data (~75%) and test data (~25%)
set.seed(100)
train_idx <- stratified_split(y = xdata$response, p_train = 0.75)
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

# tune hyperparameters
ctrl_binary <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  search = "random",
  summaryFunction = twoClassSummary, # summary for binary response
  savePredictions = "final"
)

# train binary RF with 5-fold cross-validation
set.seed(1) # set seed
rf_cv <- train(
  x = train_data[, -1], # variables to use for training
  y = train_data[, 1], # response variable
  method = "ranger",
  importance = "permutation",
  num.trees = 1000,
  trControl = ctrl_binary,
  tuneLength = 30,
  metric = "ROC"
)

# save rf object
saveRDS(rf_cv, file.path(output_folder, "rfObject.rds"))

# save summary of trained RF
sink(file.path(output_folder, "rf_cv_results.txt"))
print(rf_cv)
sink()

# save plot of hyperparameters tuning
resol <- 300
png(file.path(output_folder, "hyperpar_tuning.png"),
  width = 6 * resol, height = 4 * resol, res = resol
)
# plot(rf_cv) |> print()
ggplot(rf_cv) |> print()
# TODO: improve colors
dev.off()

# Examine prediction performance ------------------------------------------

# predicted classes (NR/R)
pred_classes <- predict(
  rf_cv,
  newdata = test_data[, -1] # test data set (without response variable)
)

# confusion matrix and performance metrics
conf_matrix <- confusionMatrix(
  data = pred_classes,
  reference = test_data[, 1], # true response in test data set
  positive = "R"
)

# save report of confusion matrix
sink(file.path(output_folder, "CM_metrics.txt"))
print(conf_matrix)
sink()

# predicted probabilities of NR and R
pred_probs <- predict(
  rf_cv,
  newdata = test_data[, -1], # test data set (without response variable)
  type = "prob"
)

# ROC
xroc <- roc(
  response = test_data[, 1], # true response in test data set
  predictor = pred_probs[, "R"], # probability of positive level
  quiet = FALSE
)

# ROC plot with AUC
png(file.path(output_folder, "ROC_plot.png"),
  width = 4 * resol, height = 4 * resol, res = resol
)
plot(xroc,
  print.auc = TRUE,
  auc.polygon = TRUE,
  xaxs = "i", yaxs = "i", las = 1
)
dev.off()


# Get variable importance -------------------------------------------------

# extract variable importance
vimp <- varImp(rf_cv)

# prepare dataframe of sorted variable importance scores
vimp_df <- vimp$importance
vimp_df <- vimp_df[order(vimp_df$Overall, decreasing = TRUE), , drop = FALSE]
vimp_df$Rank <- rank(vimp_df$Overall)
vimp_df$Rank <- abs( # invert rank
  vimp_df$Rank - max(vimp_df$Rank) - 1
)

# save dataframe
write.table(vimp_df, file.path(output_folder, "variable_importance.txt"),
  row.names = FALSE, sep = "\t"
)


# Plot settings -----------------------------------------------------------
resol <- 300 # resolution in ppi

# prepare colors for response (binary)
colors_response <- paletteer_c("grDevices::RdYlBu", nlevels(xdata$response))
names(colors_response) <- levels(xdata$response)


# Barplot of top importance scores ----------------------------------------
ntop <- 50
plotdata <- vimp_df$Overall
names(plotdata) <- rownames(vimp_df)
plotdata <- head(plotdata, ntop) |> rev()
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
plotdata2 <- data.frame(
  genes = rownames(vimp_df),
  importance = vimp_df$Overall,
  importance_rank = vimp_df$Rank
) |> head(ntop)
# remove gene subfix (version)
plotdata2$genes <- gsub("\\.[A-Z0-9]*", "", plotdata2$genes)
# load and merge LogFC data
logfc <- read.delim("nonsync/autogo_response_pre_icb/NR_vs_R/DE_NR_vs_R_allres.tsv")
plotdata2 <- merge(plotdata2, logfc, all.x = TRUE, sort = FALSE)
rownames(plotdata2) <- plotdata2$genes
plotdata2 <- plotdata2[order(plotdata2$importance_rank, decreasing = TRUE), ]
# prepare colors
xpal <- paletteer_c("grDevices::RdBu", 251, -1)
col_breaks <- seq(
  from = -max(abs(plotdata2$log2FoldChange), na.rm = TRUE),
  to = max(abs(plotdata2$log2FoldChange), na.rm = TRUE),
  length.out = length(xpal)
)
col_fun <- colorRamp2(breaks = col_breaks, colors = xpal)
# barplot
png(file.path(output_folder, paste0(
  "vimp_barplot_LogFC_top", ntop, ".png"
)), width = 5 * resol, height = 10 * resol, res = resol)
par(mar = c(5, 8, 0, 0.5), las = 2)
barplot(plotdata,
  horiz = TRUE,
  col = ifelse(
    is.na(plotdata2$log2FoldChange), "grey40", col_fun(plotdata2$log2FoldChange)
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

# put metrics together
metrics <- c(
  AUC = xroc$auc |> as.numeric(),
  conf_matrix$overall[c("Accuracy", "Kappa")],
  conf_matrix$byClass
)

# prepare metrics to plot
plotdata <- data.frame(
  metric_type = factor(names(metrics), levels = c(
    "AUC", "Accuracy", "Kappa", "Sensitivity", "Specificity",
    "Pos Pred Value", "Neg Pred Value", "Precision", "Recall",
    "F1", "Prevalence", "Detection Rate", "Detection Prevalence",
    "Balanced Accuracy"
  )),
  value = metrics,
  row.names = NULL
)
plotdata <- subset( # relevant subset of metrics
  plotdata, metric_type %in% c(
    "AUC", "Accuracy", "Sensitivity", "Specificity",
    "Pos Pred Value", "Neg Pred Value", "Precision"
  )
) |> droplevels()

# prepare coordinates for barplot
xx <- barplot(
  value ~ metric_type,
  data = plotdata, plot = FALSE
)
rownames(xx) <- levels(plotdata$metric_type)

# create barplot
png(file.path(output_folder, "barplot_metrics.png"),
  width = 4 * resol, height = 3.5 * resol, res = resol
)
par(las = 2, xpd = TRUE, tcl = -0.3, mar = c(7, 3.5, 0.8, 0.2), mgp = c(2.5, 0.8, 0))
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
xx <- conf_matrix$table |> as.data.frame()
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
  filename = file.path(output_folder, "CM_plot.png"),
  width = 3, height = 3, units = "in", dpi = 300
)

# Save image --------------------------------------------------------------

save.image("nonsync/07_rf_transcriptomics.RData")
