# load("nonsync/06a_rf_metadataset.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(ranger) # v0.17.0
  library(caret) # v7.0-1
  library(ordinalForest) # v2.4-4
  library(MLmetrics) # v1.1.3
  library(pROC) # v1.19.0.1
  library(grid) # v4.5.2
  library(ggplot2)
})


# Prepare new directories -------------------------------------------------

# main output folder
output_folder <- "nonsync/06_rf_metadataset"
if (!dir.exists(output_folder)) dir.create(output_folder)

# define combinations of rf models to be created
to_do <- expand.grid(
  response_type = c("binary", "ordinal3", "ordinal6"),
  age_and_gender = c(TRUE, FALSE),
  technical_predictors = c(FALSE, TRUE),
  hed_data = TRUE
)
# add RFs with only CIBERSORTx data
to_do <- rbind(
  to_do, data.frame(
    response_type = c("binary", "ordinal3", "ordinal6"),
    age_and_gender = FALSE,
    technical_predictors = FALSE,
    hed_data = FALSE
  )
)
# add RFs with CIBERSORTx + age and gender
to_do <- rbind(
  to_do, data.frame(
    response_type = c("binary", "ordinal3", "ordinal6"),
    age_and_gender = TRUE,
    technical_predictors = FALSE,
    hed_data = FALSE
  )
)
# define rf formulae names
to_do$rf_formula <- paste(
  paste0("RF", 1:6) |> rep(each = 3),
  c("bin", "ord3", "ord6"),
  sep = "_"
)

# create new folders
to_do$folder <- character(nrow(to_do))
for (i in seq_len(nrow(to_do))) {
  to_do$folder[i] <- file.path(
    output_folder, # main output folder
    to_do$rf_formula[i]
  )
  if (!dir.exists(to_do$folder[i])) dir.create(to_do$folder[i])
}


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

# create a list of dataframes with the full datasets for each row of to_do
df_all_list <- lapply(
  X = seq_len(nrow(to_do)), FUN = function(i) {
    # names of all variables
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
    # exclude HED when appropriate
    if (!to_do$hed_data[i]) vars_to_include <- grepv("^HED", vars_to_include, invert = TRUE)
    # exclude incomplete cases
    xx <- alldata[, vars_to_include] |>
      na.exclude() |>
      droplevels()
    # uniform response variable name
    names(xx) <- gsub("response_.levels", "response", names(xx))
    # return dataframe
    return(xx)
  }
)
names(df_all_list) <- to_do$rf_formula

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

# get training data (~75%) and test data (~25%) for each case
train_idx_list <- lapply(
  X = seq_len(nrow(to_do)), FUN = function(i) {
    set.seed(100) # set same seed each time
    stratified_split(
      y = df_all_list[[i]]$response,
      p_train = 0.75
    )
  }
)
df_train_list <- lapply(
  X = seq_len(nrow(to_do)), FUN = function(i) {
    df_all_list[[i]][train_idx_list[[i]], ]
  }
)
df_test_list <- lapply(
  X = seq_len(nrow(to_do)), FUN = function(i) {
    df_all_list[[i]][-train_idx_list[[i]], ]
  }
)
names(df_train_list) <- names(df_test_list) <- to_do$rf_formula

# save data in each folder
for (i in seq_len(nrow(to_do))) {
  write.csv(x = df_train_list[[i]], file = file.path(
    to_do$folder[i], paste0("data_train_", to_do$rf_formula[i], ".csv")
  ), row.names = FALSE)
  write.csv(x = df_test_list[[i]], file = file.path(
    to_do$folder[i], paste0("data_test_", to_do$rf_formula[i], ".csv")
  ), row.names = FALSE)
}

# define formulae
formulae_list <- lapply(
  X = seq_len(nrow(to_do)), FUN = function(i) {
    predictors <- names(df_all_list[[i]]) |>
      setdiff(y = c("accession", "response")) |>
      paste(collapse = " + ")
    return(paste0("response ~ ", predictors))
  }
)
names(formulae_list) <- to_do$rf_formula

# save formula and number of observations in test/train data in each folder
for (i in seq_len(nrow(to_do))) {
  n_train <- nrow(df_train_list[[i]])
  n_test <- nrow(df_test_list[[i]])
  sink(file.path(to_do$folder[i], paste0(
    "formula_sampleSizes_", to_do$rf_formula[i], ".txt"
  )))
  cat(paste0("===== ", to_do$rf_formula[i], " =====\n"), sep = "\n")
  cat("===== Formula =====", sep = "\n")
  cat(formulae_list[[i]], "\n")
  cat("\n===== Excluded data (NAs present) =====", sep = "\n")
  cat(paste0(nrow(alldata) - n_train - n_test, " observations"), sep = "\n")
  cat("\n===== Training data =====", sep = "\n")
  cat(paste0(
    n_train, " observations (", round(n_train / (n_train + n_test), 3) * 100,
    "%)"
  ), sep = "\n")
  cat("\nResponse levels:")
  print(df_train_list[[i]]$response |> table())
  cat("\n===== Test data =====", sep = "\n")
  cat(paste0(
    n_test, " observations (", round(n_test / (n_train + n_test), 3) * 100,
    "%)"
  ), sep = "\n")
  cat("\nResponse levels:")
  print(df_test_list[[i]]$response |> table())
  sink()
}



# Train binary random forests ---------------------------------------------

# tune hyperparameters for binary RFs of class ranger
ctrl_binary <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  search = "random",
  summaryFunction = twoClassSummary, # summary for binary response
  savePredictions = "final"
)
# train binary RFs with 5-fold cross-validation
rf_cv_binary_list <- lapply(
  X = which(to_do$response_type == "binary"), FUN = function(i) {
    set.seed(1) # set seed inside the loop
    train(
      form = as.formula(formulae_list[[i]]),
      data = df_train_list[[i]],
      method = "ranger",
      importance = "permutation",
      num.trees = 1000,
      trControl = ctrl_binary,
      tuneLength = 30,
      metric = "ROC"
    )
  }
)
names(rf_cv_binary_list) <- to_do$rf_formula[which(to_do$response_type == "binary")]

# save rf objects
for (i in seq_along(rf_cv_all)) {
  j <- which(to_do$rf_formula == names(rf_cv_binary_list)[i])
  saveRDS(
    object = rf_cv_binary_list[[i]],
    file = file.path(to_do$folder[j], paste0("rfObject_", to_do$rf_formula[j], ".rds"))
  )
}

# Train ordinal random forests --------------------------------------------

# ## ordinal random forests ##
# # tune hyperparameters for ordinal RFs of class ordfor
# ctrl_ordinal <- trainControl(
#   method = "cv",
#   number = 5,
#   classProbs = TRUE,
#   search = "random",
#   summaryFunction = multiClassSummary, # summary for multi-class response
#   savePredictions = "final"
# )
# # train ordinal RFs with 5-fold cross-validation
# rf_cv_ordinal_list <- lapply(
#   X = which(to_do$response_type != "binary"), FUN = function(i) {
#     set.seed(1) # set seed inside the loop
#     train(
#       form = as.formula(formulae_list[[i]]),
#       data = df_train_list[[i]],
#       method = "ordinalRF",
#       trControl = ctrl_ordinal,
#       tuneLength = 15
#     )
#   }
# )
# names(rf_cv_ordinal_list) <- to_do$rf_formula[which(to_do$response_type != "binary")]
#
# # put objects together in the same order as to_do
# rf_cv_all <- c(rf_cv_binary_list, rf_cv_ordinal_list)[to_do$rf_formula]
#
# # save rf objects
# for (i in seq_len(nrow(to_do))) {
#   saveRDS(
#     object = rf_cv_all[[i]],
#     file = file.path(to_do$folder[i], paste0("rfObject_", to_do$rf_formula[i], ".rds"))
#   )
# }


# Save tuned hyperparameters ----------------------------------------------
resol <- 300 # plot resolution in PPI

for (i in which(to_do$response_type == "binary")) {
  # save summary of trained RF
  sink(file.path(to_do$folder[i], paste0(
    "rf_cv_results_", to_do$rf_formula[i], ".txt"
  )))
  cat(paste0("===== ", to_do$rf_formula[i], " =====\n"), sep = "\n")
  print(rf_cv_binary_list[[to_do$rf_formula[i]]])
  sink()
  # save plot of hyperparameters tuning
  png(file.path(to_do$folder[i], paste0(
    "hyperpar_tuning_", to_do$rf_formula[i], ".png"
  )), width = 6 * resol, height = 4 * resol, res = resol)
  # plot(rf_cv_binary_list[[to_do$rf_formula[i]]]) |> print()
  ggplot(rf_cv_binary_list[[to_do$rf_formula[i]]]) |> print()
  # TODO: improve colors
  dev.off()
}


# Confusion Matrix and performance metrics --------------------------------
confusion_matrices <- lapply(
  X = which(to_do$response_type == "binary"), FUN = function(i) {
    xxformula <- to_do$rf_formula[i]
    pred_classes <- predict(
      rf_cv_binary_list[[xxformula]],
      newdata = df_test_list[[i]]
    )
    confmat <- confusionMatrix(
      data = pred_classes,
      reference = df_test_list[[i]]$response,
      positive = "R"
    )
    return(confmat)
  }
)
names(confusion_matrices) <- to_do$rf_formula[which(to_do$response_type == "binary")]

# save report of confusion matrix
for (i in seq_along(confusion_matrices)) {
  j <- which(to_do$rf_formula == names(confusion_matrices)[i])
  sink(file.path(to_do$folder[j], paste0(
    "CM_metrics_", to_do$rf_formula[j], ".txt"
  )))
  cat(paste0("===== ", to_do$rf_formula[j], " =====\n"), sep = "\n")
  print(confusion_matrices[[i]])
  sink()
}


# ROC and AUC -------------------------------------------------------------
rocs_list <- lapply(
  X = which(to_do$response_type == "binary"), FUN = function(i) {
    xxformula <- to_do$rf_formula[i]
    pred_probs <- predict(
      rf_cv_binary_list[[xxformula]],
      newdata = df_test_list[[i]],
      type = "prob"
    )
    xroc <- roc(
      response = df_test_list[[i]]$response,
      predictor = pred_probs[, "R"], # probability of positive level
      quiet = TRUE
    )
    return(xroc)
  }
)
names(rocs_list) <- to_do$rf_formula[which(to_do$response_type == "binary")]

# Variable importance -----------------------------------------------------
vimp_list <- lapply(
  rf_cv_binary_list, function(x) varImp(x)
)

# get names of all variables
allvars <- Reduce(
  f = union,
  x = lapply(vimp_list, function(x) rownames(x$importance))
) |>
  # adjust names where needed
  gsub(pattern = "^gender", replacement = "gender_") |>
  gsub(pattern = "^enrichment_protocol", replacement = "enrichm_prot_") |>
  gsub(pattern = "^dataset", replacement = "dataset_") |>
  # alphabetical order
  sort()

# keep same variable order in each RF
vimp_list <- lapply(
  X = rf_cv_binary_list, FUN = function(x) {
    vimp <- varImp(x)$importance
    xx <- as.vector(vimp$Overall)
    names(xx) <- rownames(vimp) |>
      gsub(pattern = "^gender", replacement = "gender_") |>
      gsub(pattern = "^enrichment_protocol", replacement = "enrichm_prot_") |>
      gsub(pattern = "^dataset", replacement = "dataset_")
    xx <- xx[allvars]
    names(xx) <- allvars
    return(xx)
  }
)

# create dataframe of variable importance scores
vimp_df <- do.call(rbind, vimp_list)

# save importance scores
for (i in seq_along(vimp_list)) {
  j <- which(to_do$rf_formula == names(vimp_list)[i])
  # prepare ordered dataframe
  xx <- data.frame(
    variable = names(vimp_list[[i]]),
    importance = round(vimp_list[[i]], 1),
    rank = rank(vimp_list[[i]], na.last = TRUE)
  )
  xx$rank[is.na(xx$importance)] <- NA
  xx <- xx[order(xx$rank, na.last = TRUE, decreasing = TRUE), ]
  rownames(xx) <- NULL
  # save dataframe
  write.table(xx, file.path(to_do$folder[j], paste0(
    "variable_importance_", to_do$rf_formula[j], ".txt"
  )),
  row.names = FALSE, sep = "\t"
  )
}


# Save image --------------------------------------------------------------
rm(xx, i, j, n_test, n_train)
save.image("nonsync/06a_rf_metadataset.RData")
