# load("nonsync/06_random_forests.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
library(ranger) # v0.17.0
library(ordinalForest) # v2.4-4
source("rf_functions.R")
# library(ggplot2)
# library(irr) # v0.84.1
# library(pROC) # v1.19.0.1


# Prepare new directories -------------------------------------------------

# main output folder
output_folder <- "nonsync/06_random_forests"
if (!dir.exists(output_folder)) dir.create(output_folder)

# define combinations of rf models to be created
to_do <- expand.grid(
  response_type = c("binary", "ordinal"),
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
  "accession", "response", "response_group",
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
str(alldata$response) # ordered factor (6 levels)
str(alldata$response_group) # ordered factor (2 levels)
alldata$response_group <- factor( # response_group non-ordered
  alldata$response_group,
  ordered = FALSE
)

# define new dataframe for each row of to_do
to_do$df_all <- paste(
  "df", "all", substr(to_do$response_type, 1, 3),
  as.numeric(to_do$technical_predictors),
  as.numeric(to_do$age_and_gender),
  sep = "_"
)
for (i in seq_len(nrow(to_do))) {
  vars_to_include <- names(alldata)
  # exclude columns
  if (to_do$response_type[i] == "binary") { # select binary response
    vars_to_include <- setdiff(vars_to_include, c("response"))
  } else { # or select ordinal response
    vars_to_include <- setdiff(vars_to_include, c("response_group"))
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
  names(xx) <- gsub("response_group", "response", names(xx))
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
for (i in which(to_do$response_type == "ordinal")) {
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


# Save image --------------------------------------------------------------

save.image("nonsync/06_random_forests.RData")
