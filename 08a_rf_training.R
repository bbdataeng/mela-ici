# load("nonsync/08a_RF_metadataset.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(ranger) # v0.18.0
  library(caret) # v7.0-1
  library(ordinalForest) # v2.4-4
  library(MLmetrics) # v1.1.3
  library(pROC) # v1.19.0.1
  library(grid) # v4.6.0
  library(ggplot2) # v4.0.3
})

# Prepare output folders --------------------------------------------------

# main output folders
outdirs <- c(
  # output folder for complete dataset
  complete = "nonsync/08_RF/complete",
  # output folder for dataset without checkmate067
  nocheckmate067 = "nonsync/08_RF/nocheckmate067"
)

# create directories
for (i in seq_along(outdirs)) {
  if (!dir.exists(outdirs[[i]])) dir.create(outdirs[[i]], recursive = TRUE)
}

# prepare instructions
to_do_list <- lapply(outdirs, function(x) {
  # define combinations of RF models to be created
  to_do <- expand.grid(
    response_type = c("binary", "ordinal3", "ordinal4"),
    age_and_gender = c(TRUE, FALSE),
    technical_predictors = c(FALSE, TRUE),
    hed_data = TRUE
  )
  # add RFs with only CIBERSORTx data
  to_do <- rbind(
    to_do, data.frame(
      response_type = c("binary", "ordinal3", "ordinal4"),
      age_and_gender = FALSE,
      technical_predictors = FALSE,
      hed_data = FALSE
    )
  )
  # add RFs with CIBERSORTx + age and gender
  to_do <- rbind(
    to_do, data.frame(
      response_type = c("binary", "ordinal3", "ordinal4"),
      age_and_gender = TRUE,
      technical_predictors = FALSE,
      hed_data = FALSE
    )
  )
  # define rf formulae names
  to_do$rf_formula <- paste(
    paste0("RF", 1:6) |> rep(each = 3),
    c("bin", "ord3", "ord4"),
    sep = "_"
  )
  # create new folders
  to_do$folder <- character(nrow(to_do))
  for (i in seq_len(nrow(to_do))) {
    to_do$folder[i] <- file.path(
      x, # main output folder
      to_do$rf_formula[i]
    )
    if (!dir.exists(to_do$folder[i])) dir.create(to_do$folder[i])
  }
  to_do # return object
})


# Load and prepare data ---------------------------------------------------

# load complete data
metadata_complete <- readRDS("nonsync/04_clean_data/clean_metadata.rds")
xdata_complete <- readRDS("nonsync/04_clean_data/clean_cibersortx.rds")
hed_complete <- readRDS("nonsync/04_clean_data/clean_hed.rds")

# identify checkmate067 cohort
to_exclude <- which(
  metadata_complete$dataset == "Campbell-2023" &
    metadata_complete$enrichment_protocol == "targeted-mRNA-capture"
)

# keep only necessary metadata columns
metadata_complete <- metadata_complete[, c(
  "accession", "response_4levels", "response_3levels", "response_2levels",
  "age", "gender", "enrichment_protocol", "dataset"
)]
metadata_complete$cibersortx_Absolute_Score <- xdata_complete$Absolute_Score

# remove mean HED column
hed_complete <- hed_complete[, setdiff(names(hed_complete), "HED_mean")]

# sanitize cell-type column names
xdata_complete <- xdata_complete[, names(xdata_complete) != "Absolute_Score"]
cell_types_original <- colnames(xdata_complete)
cell_types <- gsub(" ", "_", cell_types_original)
cell_types <- gsub("\\(|\\)", "", cell_types)
colnames(xdata_complete) <- cell_types
names(cell_types_original) <- cell_types

# merge metadata, HED, and LM22 and create a list of data
alldata <- list(
  complete = cbind(metadata_complete, hed_complete, xdata_complete),
  nocheckmate067 = cbind(metadata_complete, hed_complete, xdata_complete)[-to_exclude, ]
)
rm(metadata_complete, hed_complete, xdata_complete, to_exclude)

# create a list of dataframes with the full datasets for each row of to_do
df_all_list <- lapply(seq_along(to_do_list), function(x) {
  dfs <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    # names of all variables
    vars_to_include <- names(alldata[[x]])
    # exclude columns
    if (to_do_list[[x]]$response_type[i] == "binary") { # select binary response
      vars_to_include <- setdiff(vars_to_include, c("response_4levels", "response_3levels"))
    } else if (to_do_list[[x]]$response_type[i] == "ordinal3") { # select 3-level ordinal response
      vars_to_include <- setdiff(vars_to_include, c("response_4levels", "response_2levels"))
    } else { # select 4-level ordinal response
      vars_to_include <- setdiff(vars_to_include, c("response_3levels", "response_2levels"))
    }
    # exclude technical predictors when appropriate
    if (!to_do_list[[x]]$technical_predictors[i]) {
      vars_to_include <- setdiff(
        vars_to_include, c("enrichment_protocol", "dataset")
      )
    }
    # exclude age and gender when appropriate
    if (!to_do_list[[x]]$age_and_gender[i]) {
      vars_to_include <- setdiff(
        vars_to_include, c("age", "gender")
      )
    }
    # exclude HED when appropriate
    if (!to_do_list[[x]]$hed_data[i]) vars_to_include <- grepv("^HED", vars_to_include, invert = TRUE)
    # exclude incomplete cases
    xx <- alldata[[x]][, vars_to_include] |>
      na.exclude() |>
      droplevels()
    # uniform response variable name
    names(xx) <- gsub("response_.levels", "response", names(xx))
    # return dataframe
    xx
  })
  names(dfs) <- to_do_list[[x]]$rf_formula
  dfs
})
names(df_all_list) <- names(to_do_list)

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

# get training data (~70%) and test data (~30%) for each case
train_idx_list <- lapply(seq_along(to_do_list), function(x) {
  indices <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    set.seed(100) # set same seed each time
    stratified_split(
      y = df_all_list[[x]][[i]]$response,
      p_train = 0.7
    )
  })
  names(indices) <- to_do_list[[x]]$rf_formula
  indices
})
names(train_idx_list) <- names(to_do_list)

df_train_list <- lapply(seq_along(to_do_list), function(x) {
  dfs <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    df_all_list[[x]][[i]][train_idx_list[[x]][[i]], ]
  })
  names(dfs) <- to_do_list[[x]]$rf_formula
  dfs
})
names(df_train_list) <- names(to_do_list)

df_test_list <- lapply(seq_along(to_do_list), function(x) {
  dfs <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    df_all_list[[x]][[i]][-train_idx_list[[x]][[i]], ]
  })
  names(dfs) <- to_do_list[[x]]$rf_formula
  dfs
})
names(df_test_list) <- names(to_do_list)

# save data in each folder
for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    write.csv(x = df_train_list[[x]][[i]], file = file.path(
      to_do_list[[x]]$folder[i], paste0("data_train_", to_do_list[[x]]$rf_formula[i], ".csv")
    ), row.names = FALSE)
    write.csv(x = df_test_list[[x]][[i]], file = file.path(
      to_do_list[[x]]$folder[i], paste0("data_test_", to_do_list[[x]]$rf_formula[i], ".csv")
    ), row.names = FALSE)
  }
}

# define formulae
formulae_list <- lapply(seq_along(to_do_list), function(x) {
  forms <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    predictors <- names(df_all_list[[x]][[i]]) |>
      setdiff(y = c("accession", "response")) |>
      paste(collapse = " + ")
    paste0("response ~ ", predictors)
  })
  names(forms) <- to_do_list[[x]]$rf_formula
  forms
})
names(formulae_list) <- names(to_do_list)

# save formula and number of observations in test/train data in each folder
for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    n_train <- nrow(df_train_list[[x]][[i]])
    n_test <- nrow(df_test_list[[x]][[i]])
    sink(file.path(to_do_list[[x]]$folder[i], paste0(
      "formula_sampleSizes_", to_do_list[[x]]$rf_formula[i], ".txt"
    )))
    cat(paste0("===== ", to_do_list[[x]]$rf_formula[i], " =====\n"), sep = "\n")
    cat("===== Formula =====", sep = "\n")
    cat(formulae_list[[x]][[i]], "\n")
    cat("\n===== Excluded data (NAs present) =====", sep = "\n")
    cat(paste0(nrow(alldata[[x]]) - n_train - n_test, " observations"), sep = "\n")
    cat("\n===== Training data =====", sep = "\n")
    cat(paste0(
      n_train, " observations (", round(n_train / (n_train + n_test), 3) * 100,
      "%)"
    ), sep = "\n")
    cat("\nResponse levels:")
    print(df_train_list[[x]][[i]]$response |> table())
    cat("\n===== Test data =====", sep = "\n")
    cat(paste0(
      n_test, " observations (", round(n_test / (n_train + n_test), 3) * 100,
      "%)"
    ), sep = "\n")
    cat("\nResponse levels:")
    print(df_test_list[[x]][[i]]$response |> table())
    sink()
  }
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
rf_cv_binary_list <- lapply(seq_along(to_do_list), function(x) {
  trains <- lapply(which(to_do_list[[x]]$response_type == "binary"), function(i) {
    set.seed(1) # set seed inside the loop
    train(
      form = as.formula(formulae_list[[x]][[i]]),
      data = df_train_list[[x]][[i]],
      method = "ranger",
      importance = "permutation",
      num.trees = 1000,
      trControl = ctrl_binary,
      tuneLength = 30,
      metric = "ROC"
    )
  })
  names(trains) <- to_do_list[[x]]$rf_formula[which(to_do_list[[x]]$response_type == "binary")]
  trains
})
names(rf_cv_binary_list) <- names(to_do_list)


# save rf objects
for (x in seq_along(to_do_list)) {
  for (i in seq_along(rf_cv_binary_list[[x]])) {
    j <- which(to_do_list[[x]]$rf_formula == names(rf_cv_binary_list[[x]])[i])
    saveRDS(
      object = rf_cv_binary_list[[x]][[i]],
      file = file.path(
        to_do_list[[x]]$folder[j],
        paste0("rfObject_", to_do_list[[x]]$rf_formula[j], ".rds")
      )
    )
  }
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
#   X = which(to_do_list[[x]]$response_type != "binary"), FUN = function(i) {
#     set.seed(1) # set seed inside the loop
#     train(
#       form = as.formula(formulae_list[[x]][[i]]),
#       data = df_train_list[[x]][[i]],
#       method = "ordinalRF",
#       trControl = ctrl_ordinal,
#       tuneLength = 15
#     )
#   }
# )
# names(rf_cv_ordinal_list) <- to_do_list[[x]]$rf_formula[which(to_do_list[[x]]$response_type != "binary")]
#
# # put objects together in the same order as to_do
# rf_cv_all <- c(rf_cv_binary_list, rf_cv_ordinal_list)[to_do_list[[x]]$rf_formula]
#
# # save rf objects
# for (i in seq_len(nrow(to_do))) {
#   saveRDS(
#     object = rf_cv_all[[i]],
#     file = file.path(to_do_list[[x]]$folder[i], paste0("rfObject_", to_do_list[[x]]$rf_formula[i], ".rds"))
#   )
# }


# Save tuned hyperparameters ----------------------------------------------
resol <- 300 # plot resolution in PPI

for (x in seq_along(to_do_list)) {
  for (i in which(to_do_list[[x]]$response_type == "binary")) {
    # save summary of trained RF
    sink(file.path(to_do_list[[x]]$folder[i], paste0(
      "rf_cv_results_", to_do_list[[x]]$rf_formula[i], ".txt"
    )))
    cat(paste0("===== ", to_do_list[[x]]$rf_formula[i], " =====\n"), sep = "\n")
    print(rf_cv_binary_list[[x]][[to_do_list[[x]]$rf_formula[i]]])
    sink()
    # save plot of hyperparameters tuning
    png(file.path(to_do_list[[x]]$folder[i], paste0(
      "hyperpar_tuning_", to_do_list[[x]]$rf_formula[i], ".png"
    )), width = 6 * resol, height = 4 * resol, res = resol)
    ggplot(rf_cv_binary_list[[x]][[to_do_list[[x]]$rf_formula[i]]]) |> print()
    dev.off()
  }
}


# Confusion Matrix and performance metrics --------------------------------
confusion_matrices <- lapply(seq_along(to_do_list), function(x) {
  cms <- lapply(which(to_do_list[[x]]$response_type == "binary"), function(i) {
    xxformula <- to_do_list[[x]]$rf_formula[i]
    pred_classes <- predict(
      rf_cv_binary_list[[x]][[xxformula]],
      newdata = df_test_list[[x]][[i]]
    )
    confusionMatrix(
      data = pred_classes,
      reference = df_test_list[[x]][[i]]$response,
      positive = "R"
    )
  })
  names(cms) <- to_do_list[[x]]$rf_formula[which(to_do_list[[x]]$response_type == "binary")]
  cms
})
names(confusion_matrices) <- names(to_do_list)

# save report of confusion matrices
for (x in seq_along(to_do_list)) {
  for (i in seq_along(confusion_matrices[[x]])) {
    j <- which(to_do_list[[x]]$rf_formula == names(confusion_matrices[[x]])[i])
    sink(file.path(to_do_list[[x]]$folder[j], paste0(
      "CM_metrics_", to_do_list[[x]]$rf_formula[j], ".txt"
    )))
    cat(paste0("===== ", to_do_list[[x]]$rf_formula[j], " =====\n"), sep = "\n")
    print(confusion_matrices[[x]][[i]])
    sink()
  }
}


# ROC and AUC -------------------------------------------------------------
rocs_list <- lapply(seq_along(to_do_list), function(x) {
  xroc <- lapply(which(to_do_list[[x]]$response_type == "binary"), function(i) {
    xxformula <- to_do_list[[x]]$rf_formula[i]
    pred_probs <- predict(
      rf_cv_binary_list[[x]][[xxformula]],
      newdata = df_test_list[[x]][[i]],
      type = "prob"
    )
    roc(
      response = df_test_list[[x]][[i]]$response,
      predictor = pred_probs[, "R"], # probability of positive level
      quiet = TRUE
    )
  })
  names(xroc) <- to_do_list[[x]]$rf_formula[which(to_do_list[[x]]$response_type == "binary")]
  xroc
})
names(rocs_list) <- names(to_do_list)


# Variable importance -----------------------------------------------------
vimp_list <- lapply(seq_along(to_do_list), function(x) {
  lapply(rf_cv_binary_list[[x]], varImp)
})
names(vimp_list) <- names(to_do_list)

# get names of all variables
allvars <- Reduce(
  f = union,
  x = lapply(vimp_list[[1]], function(x) rownames(x$importance))
) |>
  # adjust names where needed
  gsub(pattern = "^gender", replacement = "gender_") |>
  gsub(pattern = "^enrichment_protocol", replacement = "enrichm_prot_") |>
  gsub(pattern = "^dataset", replacement = "dataset_") |>
  # alphabetical order
  sort()

# keep same variable order in each RF
vimp_list <- lapply(rf_cv_binary_list, function(a) {
  lapply(a, function(b) {
    vimp <- varImp(b)$importance
    xx <- as.vector(vimp$Overall)
    names(xx) <- rownames(vimp) |>
      gsub(pattern = "^gender", replacement = "gender_") |>
      gsub(pattern = "^enrichment_protocol", replacement = "enrichm_prot_") |>
      gsub(pattern = "^dataset", replacement = "dataset_")
    xx <- xx[allvars]
    names(xx) <- allvars
    xx
  })
})

# create dataframes of variable importance scores
vimp_df_list <- lapply(vimp_list, do.call, what = rbind)

# save importance scores
for (x in seq_along(to_do_list)) {
  for (i in seq_along(vimp_list[[x]])) {
    j <- which(to_do_list[[x]]$rf_formula == names(vimp_list[[x]])[i])
    # prepare ordered dataframe
    xx <- data.frame(
      variable = names(vimp_list[[x]][[i]]),
      importance = round(vimp_list[[x]][[i]], 1),
      rank = rank(vimp_list[[x]][[i]], na.last = TRUE)
    )
    xx$rank[is.na(xx$importance)] <- NA
    xx <- xx[order(xx$rank, na.last = TRUE, decreasing = TRUE), ]
    rownames(xx) <- NULL
    # save dataframe
    write.table(xx, file.path(to_do_list[[x]]$folder[j], paste0(
      "variable_importance_", to_do_list[[x]]$rf_formula[j], ".txt"
    )),
    row.names = FALSE, sep = "\t"
    )
  }
}


# Save image --------------------------------------------------------------
rm(xx, i, j, n_test, n_train, x)
save.image("nonsync/08a_RF_metadataset.RData")
