# load("nonsync/06a_rf_metadataset.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(ranger) # v0.17.0
  library(caret) # v7.0-1
  library(ordinalForest) # v2.4-4
  library(pROC) # v1.19.0.1
  library(grid) # v4.5.2
  source("rf_functions.R")
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
# define rf formulae names
to_do$rf_formula <- paste(
  paste0("RF", 1:5) |> rep(each = 3),
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

# get training data (~80%) and test data (~20%) for each case
set.seed(123)
train_idx_list <- lapply(
  X = seq_len(nrow(to_do)), FUN = function(i) {
    stratified_split(
      y = df_all_list[[i]]$response,
      p_train = 0.8
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


# Create random forests ---------------------------------------------------
rf_list <- lapply(
  X = seq_len(nrow(to_do)), FUN = function(i) {
    # get predictors
    predictors <- formulae_list[[i]] |>
      gsub(pattern = "response ~ ", replacement = "") |>
      strsplit(split = " \\+ ") |>
      unlist()
    # binary random forest
    if (to_do$response_type[i] == "binary") {
      # train RF with ranger()
      xrf <- ranger(
        formula = as.formula(formulae_list[[i]]), # rf formula
        data = df_train_list[[i]], # training data
        importance = "permutation",
        probability = TRUE,
        min.node.size = 5,
        na.action = "na.fail", # no NAs expected in the training data
        respect.unordered.factors = "partition",
        num.trees = 1000,
        mtry = floor(sqrt(length(predictors))),
        seed = 1 # set seed for reproducibility
      )
    } else { # ordinal random forest
      # set seed for reproducibility
      set.seed(1)
      # train RF with ordfor()
      xrf <- ordfor(
        depvar = "response", # response variable name
        data = df_train_list[[i]][, c("response", predictors)],
        nsets = 500,
        ntreeperdiv = 100,
        ntreefinal = 3000,
        mtry = floor(sqrt(length(predictors))),
        min.node.size = 10,
        perffunction = "probability", # uses Ranked Probability Score (RPS)
        importance = "rps",
        num.threads = 0 # 0 = auto
      )
    }
    # return trained random forest
    return(xrf)
  }
)
names(rf_list) <- to_do$rf_formula

# save rf objects
for (i in seq_len(nrow(to_do))) {
  saveRDS(
    object = rf_list[[i]],
    file = file.path(to_do$folder[i], paste0("rfObject_", to_do$rf_formula[i], ".rds"))
  )
}


# Get confusion matrices --------------------------------------------------
confusion_matrices <- lapply( # list of confusion matrices
  X = seq_len(nrow(to_do)),
  FUN = function(i) {
    get_confusion_matrix(
      rf_object = rf_list[[i]],
      testdata = df_test_list[[i]],
      show_sum = TRUE
    )
  }
)
names(confusion_matrices) <- to_do$rf_formula

# Get accuracy metrics ----------------------------------------------------
accuracy_metrics <- lapply( # list of accuracy metrics
  X = seq_len(nrow(to_do)),
  FUN = function(i) {
    get_accuracy_metrics(
      rf_object = rf_list[[i]],
      testdata = df_test_list[[i]],
      confusion_matrix = confusion_matrices[[i]],
      positive_level = "R"
    )
  }
)
names(accuracy_metrics) <- to_do$rf_formula

# save confusion matrices and accuracy metrics in the respective folders
for (i in seq_len(nrow(to_do))) {
  sink(file.path(to_do$folder[i], paste0(
    "CM_accuracyMetrics_", to_do$rf_formula[i], ".txt"
  )))
  cat("===== Confusion Matrix =====", sep = "\n")
  print(confusion_matrices[[i]])
  cat("\n===== Accuracy Metrics =====", sep = "\n")
  print(accuracy_metrics[[i]])
  sink()
}


# Get variable importance -------------------------------------------------
importances <- lapply( # list of variable importances
  X = seq_len(nrow(to_do)),
  FUN = function(i) get_importance(rf_object = rf_list[[i]])
)
names(importances) <- to_do$rf_formula

# create dataframe of variable importance scores
vars <- names(alldata) |> grepv(pattern = "response|accession", invert = TRUE)
for (i in seq_along(importances)) {
  importances[[i]] <- importances[[i]][vars]
  names(importances[[i]]) <- vars
}
importances_df <- do.call(rbind, importances)

# prepare function for min-max normalization
minmax_norm <- function(x) {
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  return((x - xmin) / (xmax - xmin))
}

# get transformed ranked variable importance scores
importances_rank <- lapply(
  X = importances, FUN = function(x) {
    # get rank
    xrank <- rank(x, na.last = TRUE)
    xrank[is.na(x)] <- NA
    # min-max normalization
    xrank <- minmax_norm(xrank)
    return(xrank)
  }
)

# save importance scores
for (i in seq_len(nrow(to_do))) {
  # prepare ordered dataframe
  xx <- data.frame(
    variable = names(importances[[i]]),
    importance = importances[[i]],
    rank = rank(importances[[i]], na.last = TRUE)
  )
  xx$rank[is.na(xx$importance)] <- NA
  xx$tr_rank <- minmax_norm(xx$rank) |> round(3)
  xx <- xx[order(xx$rank, na.last = TRUE, decreasing = TRUE), ]
  rownames(xx) <- NULL
  # save dataframe
  write.table(xx, file.path(to_do$folder[i], paste0(
    "variable_importance_", to_do$rf_formula[i], ".txt"
  )),
  row.names = FALSE, sep = "\t"
  )
}


# Save image --------------------------------------------------------------
rm(xx, i, n_test, n_train, vars)
save.image("nonsync/06a_rf_metadataset.RData")
