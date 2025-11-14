# Load libraries ----------------------------------------------------------
library(readxl)
library(writexl)
library(car)
library(gtools)
library(tidyverse)

# Prepare folder for cleaned data -----------------------------------------
output_folder <- "nonsync/01_clean_data"
if (!dir.exists(output_folder)) dir.create(output_folder)


# Load and rearrange data -------------------------------------------------

# load data
xdata <- read_excel("nonsync/00_raw_data/immunotherapy_datasets_metadata_cibersortx_extended.xlsx") |>
  as.data.frame() |> # convert to dataframe
  strings2factors(verbose = FALSE) # convert character columns to factors

metadata <- read_excel("nonsync/00_raw_data/immunotherapy_datasets_metadata.xlsx") |>
  as.data.frame() |> # convert to dataframe
  strings2factors(verbose = FALSE) # convert character columns to factors

# compare dimensions
dim(metadata) # 292 rows, 9 columns
dim(xdata) # 292 rows, 32 columns

# compare metadata
names(metadata)
names(xdata)
all(names(metadata) %in% names(xdata))
identical(metadata, xdata[, names(metadata)]) # identical metadata
metadata$response_group <- xdata$response_group # copy response_group

# extract CibersortX column names
cibersortx_cols <- setdiff(names(xdata), names(metadata))
cibersortx_cols

# combine data
xdata <- cbind(metadata, xdata[, cibersortx_cols])
rownames(xdata) <- xdata$accession
rm(metadata) # object no longer needed


# Extract patient_id column -----------------------------------------------

# utility function to check whether specific rows of a dataframe are identical
rows_identical <- function(df, idx) {
  if (length(idx) <= 1) stop("Only 1 index provided")
  v <- idx
  all(vapply(df, function(col) {
    ref <- col[v[1]]
    all(vapply(col[v[-1]], identical, logical(1L), ref))
  }, logical(1L)))
}

# check "sample_name"
table(table(xdata$sample_name)) # one sample_name is repeated twice
which(table(xdata$sample_name) == 2) # this is Pt27
subset(xdata, sample_name == "Pt27") # repeated rows with identical values except "accession"

# store original sample name
old_sample_name <- xdata$sample_name

# convert sample_name to character
xdata$sample_name <- as.character(xdata$sample_name)

# prepare new column patient_id with info about dataset
xdata$patient_id <- paste(xdata$sample_name, as.numeric(xdata$dataset), sep = "@")

# examine subset from dataset 1
# define columns to compare
columns_to_compare_lenient <- c("age", "gender")
columns_to_compare_strict <- c("treatment", "response", "response_group", "age", "gender", "enrichment_protocol", "dataset")

# loop over each dataset
for (xx_i in seq_along(levels(xdata$dataset))) {
  # subset data from that dataset
  xx <- subset(xdata, as.numeric(dataset) == xx_i)
  # xx = xx[order(xx$sample_name), ] # order data by sample_name
  # option 1: compare first substring before "-" or "_"
  xx_strsplit <- strsplit(xx$sample_name, "-|_") |> # split strings
    sapply(FUN = function(x) x[1]) |> # extract first bit of the strings
    unique() # get unique values

  # loop over each substring
  xx_to_check <- c()
  xx_matches_found <- FALSE
  for (xx_j in seq_along(xx_strsplit)) {
    xx_match_indices <- grep(paste0("^", xx_strsplit[xx_j], "(_|-|$)"), xx$sample_name)
    # if 2+ matches
    if (length(xx_match_indices) > 1) {
      # if rows are identical (strict rule)
      if (rows_identical(df = xx[, columns_to_compare_strict], idx = xx_match_indices)) {
        xx_matches_found <- TRUE
        xx$patient_id[xx_match_indices] <- paste0("reconstr_strict_", xx_strsplit[xx_j], "@", as.numeric(xx$dataset[xx_i]))
      } else if (rows_identical(df = xx[, columns_to_compare_lenient], idx = xx_match_indices)) {
        # if rows are identical (lenient rule)
        xx_matches_found <- TRUE
        xx$patient_id[xx_match_indices] <- paste0("reconstr_lenient_", xx_strsplit[xx_j], "@", as.numeric(xx$dataset[xx_i]))
      } else { # else, if rows differ
        # store for manual inspection
        xx_to_check <- c(xx_to_check, xx_match_indices)
      }
    }
  }

  if (!xx_matches_found) { # if no matches were found, try with first 2 bits of the strings
    xx_strsplit <- strsplit(xx$sample_name, "-|_") |> # split strings
      sapply(FUN = function(x) paste(x[1], x[2], sep = "_")) |> # extract first 2 bits of the strings
      unique() # get unique values
    xx_to_check <- c()
    for (xx_j in seq_along(xx_strsplit)) {
      xx_match_indices <- grep(paste0("^", xx_strsplit[xx_j], "(_|-)"), xx$sample_name)
      # if 2+ matches
      if (length(xx_match_indices) > 1) {
        # if rows are identical (strict rule)
        if (rows_identical(df = xx[, columns_to_compare_strict], idx = xx_match_indices)) {
          xx_matches_found <- TRUE
          xx$patient_id[xx_match_indices] <- paste0("reconstr_strict_", xx_strsplit[xx_j], "@", as.numeric(xx$dataset[xx_i]))
        } else if (rows_identical(df = xx[, columns_to_compare_lenient], idx = xx_match_indices)) {
          # if rows are identical (lenient rule)
          xx_matches_found <- TRUE
          xx$patient_id[xx_match_indices] <- paste0("reconstr_lenient_", xx_strsplit[xx_j], "@", as.numeric(xx$dataset[xx_i]))
        } else { # else, if rows differ
          # store for manual inspection
          xx_to_check <- c(xx_to_check, xx_match_indices)
        }
      }
    }
  }

  # save reconstructed changes
  if (xx_matches_found) {
    xdata$patient_id[as.numeric(xdata$dataset) == xx_i] <- xx$patient_id
  }

  # save rows to be checked manually
  if (xx_i == 1) {
    to_check_df <- xx[xx_to_check, ]
  } else {
    to_check_df <- rbind(to_check_df, xx[xx_to_check, ])
  }
}
# remove temporary objects
rm(list = ls(pattern = "^xx"))

# write data to check as excel
to_check_df <- to_check_df[
  order(to_check_df$sample_name),
  setdiff(names(to_check_df), c("original_sample_name", "patient_id"))
]
write_xlsx(to_check_df, file.path(output_folder, "to_check.xlsx"))

### manual corrections ###
# reconstruct patient "reconstr_manual_62@1" (identical data expect age increasing by 1)
xdata$patient_id[xdata$sample_name %in% to_check_df$sample_name] <-
  "reconstr_manual_62@1"
# merge reconstr_strict_MGH39@1 and reconstr_strict_39@1
xdata$patient_id[xdata$patient_id %in% c("reconstr_strict_MGH39@1", "reconstr_strict_39@1")] <-
  "reconstr_manual_MGH39_39@1"
# merge reconstr_lenient_208@1 and reconstr_strict_MGH208@1
xdata$patient_id[xdata$patient_id %in% c("reconstr_lenient_208@1", "reconstr_strict_MGH208@1")] <-
  "reconstr_manual_MGH208_208@1"
# merge reconstr_strict_42@1 and reconstr_strict_MGH42@1
xdata$patient_id[xdata$patient_id %in% c("reconstr_strict_42@1", "reconstr_strict_MGH42@1")] <-
  "reconstr_manual_MGH42_42@1"

# make sample names unique
xdata$sample_name <- make.unique(xdata$sample_name)

# reorder columns
xdata <- xdata[, c(
  "accession", "patient_id", "sample_name", "response", "response_group",
  "treatment", "biopsy_time", "age", "gender", "enrichment_protocol", "dataset",
  cibersortx_cols
)]
# reorder rows
xdata <- xdata[mixedorder(xdata$patient_id), ] # by patient_id
xdata <- xdata[order(xdata$dataset), ] # then by dataset
# save excel
# write_xlsx(xdata, file.path(output_folder, "reconstructed_patients.xlsx"))



# Inspect data ------------------------------------------------------------

# data structure
str(xdata)

# number of complete cases (rows)
table(complete.cases(xdata)) |> addmargins() # 156/292 complete cases (rows)

# check "accession"
table(table(xdata$accession)) # all unique values

# check "sample_name"
table(table(xdata$sample_name)) # now unique sample_names

# check "patient_id"
table(table(xdata$patient_id)) |> addmargins()
# 188 total patients
# 112 have only 1 biopsy
# 70 patients have 2 biopsies
# 6 patients have more than 2 biopsies (max = 10 biopsies)

# check "response"
table(xdata$response, useNA = "always") |> addmargins() # no NAs present, but level "UNK"
nlevels(xdata$response) # 8 levels

# check "response_group"
table(xdata$response_group, useNA = "always") |> addmargins() # no NAs present, but level "UNK"
nlevels(xdata$response_group) # 3 levels

# check "treatment"
table(xdata$treatment, useNA = "always") |> addmargins() # no NAs present
nlevels(xdata$treatment) # 10 levels

# check "biopsy_time"
table(xdata$biopsy_time, useNA = "always") |> addmargins() # no NAs present
nlevels(xdata$biopsy_time) # 6 levels

# check "age"
table(xdata$age) |> plot(main = "age") # age distribution
sum(is.na(xdata$age)) # 109 NAs present

# check "gender"
table(xdata$gender, useNA = "always") |> addmargins() # 56 F, 100 M, 136 NA

# check "enrichment_protocol"
table(xdata$enrichment_protocol, useNA = "always") |> addmargins() # no NAs present, but level "unspecified"
nlevels(xdata$enrichment_protocol) # 4 levels

# check "dataset"
table(xdata$dataset, useNA = "always") |> addmargins() # no NAs present
nlevels(xdata$dataset) # 5 levels
table(xdata$dataset) |> plot(las = 2, main = "dataset")


# Clean data --------------------------------------------------------------

### clean 7-level response ###
xx <- xdata$response
table(xx, useNA = "always") # have a look
xx[xx == "UNK"] <- NA # turn "UNK" to NA
xdata$response_7levels <- xx |>
  droplevels() |> # remove unused "UNK" level
  factor( # reorder levels from worst to best
    levels = c("PD", "NR", "SD", "PR", "PRCR", "R", "CR"),
    ordered = TRUE
  )
table(xdata$response_7levels, useNA = "always") # have a look

### create 6-level response ###
xx <- xdata$response_7levels
xx[xx == "PRCR" & !is.na(xx)] <- "R" # turn "PRCR" to "R"
xdata$response_6levels <- droplevels(xx)
table(xdata$response_6levels, useNA = "always") # have a look

### create 3-level response ###
xx <- xdata$response_7levels
xx[xx %in% c("PD", "NR") & !is.na(xx)] <- "NR" # level 1: non-responder
xx[xx %in% c("PR", "PRCR", "R", "CR") & !is.na(xx)] <- "R" # level 3: responder
xdata$response_3levels <- droplevels(xx)
table(xdata$response_3levels, useNA = "always") # have a look

### clean binary response ###
xx <- xdata$response_group
table(xx, useNA = "always") # have a look
xx[xx == "UNK"] <- NA # turn "UNK" to NA
xdata$response_2levels <- xx |>
  droplevels() |> # remove unused "UNK" level
  factor( # reorder levels from worst to best
    levels = c("NR", "R"),
    ordered = TRUE
  )
table(xdata$response_2levels, useNA = "always") # have a look

# reorder response columns
cols_to_keep <- c(
  "accession", "patient_id", "sample_name", "response_7levels",
  "response_6levels", "response_3levels", "response_2levels",
  "treatment", "biopsy_time", "age", "gender", "enrichment_protocol",
  "dataset", cibersortx_cols
)
xdata <- xdata[, cols_to_keep]

### clean treatment ###
table(xdata$treatment, useNA = "always") # have a look
old_treatment <- xdata$treatment # save non-cleaned data
# prepare cleaned treatment variable
xdata$treatment <- as.character(xdata$treatment) |> tolower() # lower case
xdata$treatment <- gsub("[[:space:]]+", "", xdata$treatment) # remove spaces
xdata$treatment <- gsub("-", "", xdata$treatment) # remove hyphens
xdata$treatment[ # clean "anti-PD-1"
  xdata$treatment %in% c("pembrolizumab", "nivolumab", "antipd1")
] <- "anti-PD-1"
xdata$treatment[ # clean "anti-PD-L1"
  xdata$treatment %in% c("antipdl1")
] <- "anti-PD-L1"
xdata$treatment[ # clean "anti-CTLA-4"
  xdata$treatment %in% c("antictla4")
] <- "anti-CTLA-4"
xdata$treatment[ # clean "anti-PD-1 + anti-CTLA-4"
  xdata$treatment %in% c("ipilimumab+pembrolizumab", "ipilimumab+nivolumab", "antipd1+antictla4")
] <-
  "anti-PD-1 + anti-CTLA-4"
xdata$treatment <- factor(xdata$treatment, levels = c( # transform to factor
  "anti-PD-1", "anti-PD-L1", "anti-CTLA-4", "anti-PD-1 + anti-CTLA-4"
))
table(xdata$treatment, useNA = "always") # have a look: quite unbalanced

### clean biopsy_time ###
table(xdata$biopsy_time, useNA = "always") # have a look
old_biopsy_time <- xdata$biopsy_time # save non-cleaned data
# clean text
xx <- tolower(trimws(xdata$biopsy_time))
xx <- gsub("[[:space:]]+", " ", xx)
# flag dabrafenib+trametinib therapy as a separate covariate
xdata$on_dabrafenib_trametinib <- grepl(
  "dabrafenib\\+trametinib", xx
) |> as.numeric()
# clean levels
xx[grepl("^pre", xx)] <- "PRE-ICB" # clean pre-ICB treatment
xx[grepl("^on|^early", xx)] <- "ON-ICB" # clean on-ICB treatment
# transform to factor and add to xdata
xdata$biopsy_time <- factor(
  xx,
  levels = c("PRE-ICB", "ON-ICB"), ordered = TRUE
)
table(xdata$biopsy_time, useNA = "always") # have a look

### clean gender ###
table(xdata$gender, useNA = "always") # no need to clean

### clean enrichment_protocol ###
table(xdata$enrichment_protocol, useNA = "always") # have a look
xx <- xdata$enrichment_protocol |> as.character()
xx[xx == "unspecified"] <- NA # change unspecified to NA
xx[xx == "poly-A selection"] <- "polyA-selection" # clean "polyA-selection"
xx[xx == "ribo-zero depletion"] <- "rRNA-depletion" # clean "rRNA-depletion"
xx[xx == "targeted mRNA capture"] <- "targeted-mRNA-capture" # clean "targeted-mRNA-capture"
xdata$enrichment_protocol <- as.factor(xx)
table(xdata$enrichment_protocol, useNA = "always") # have a look

### clean dataset ###
table(xdata$dataset, useNA = "always") # no need to clean

### clean patient_id ###
xdata$patient_id <- as.factor(xdata$patient_id)


# Manually exclude samples ------------------------------------------------

# exclude samples done on dabrafenib-trametinib treatment
xdata <- xdata[!xdata$on_dabrafenib_trametinib, ]
# remove "on_dabrafenib_trametinib" column
xdata <- xdata[, setdiff(names(xdata), "on_dabrafenib_trametinib")]

# exclude samples with unknown response
xdata <- xdata[!is.na(xdata$response_2levels), ] |> droplevels()

# exclude biopsies made during treatment
xdata <- subset(xdata, biopsy_time == "PRE-ICB") |> droplevels()

# List of preferred samples (highest number of reads)
preferred_samples <- c("SRR7344567", "SRR7344565", "SRR3184298")

# Keep only:
# - preferred samples if patient has multiple samples
# - all other samples for patients with only one
xdata <- xdata %>%
  group_by(patient_id) %>%
  filter(
    if (n_distinct(sample_name) > 1) {
      accession %in% preferred_samples
    } else {
      TRUE
    }
  ) %>%
  ungroup() %>%
  droplevels()

table(table(xdata$patient_id)) # only 1 biopsy per patient, as expected


# Add HED data ------------------------------------------------------------

hed_data <- read.table("nonsync/00_raw_data/HED.tsv", header = TRUE) |>
  strings2factors(verbose = FALSE) # convert character columns to factors
head(hed_data)

names(hed_data)
names(hed_data) <- c("accession", paste0("HED_locus", LETTERS[1:3]), "HED_mean")

xx <- merge(xdata, hed_data, by = "accession", all.x = TRUE)
all(complete.cases(xx[, grepv("^HED", names(xx))])) # no missing values, safe to merge
xdata <- xx
rownames(xdata) <- xdata$accession
rm(xx)


# Export xdata ---------------------------------------------------------

# extract HED data
hed_data <- xdata[, grepv("^HED", names(xdata))]

# extract metadata
metadata <- xdata[, setdiff(names(xdata), c(cibersortx_cols, names(hed_data)))]

# extract cibersortx data
cibersortx_data <- xdata[, cibersortx_cols]

# as RDS
saveRDS(metadata, file.path(output_folder, "clean_metadata.rds"))
saveRDS(hed_data, file.path(output_folder, "clean_hed.rds"))
saveRDS(cibersortx_data, file.path(output_folder, "clean_cibersortx.rds"))
saveRDS(xdata, file.path(output_folder, "clean_alldata.rds"))

# as CSV
write.csv(metadata, file.path(output_folder, "clean_metadata.csv"), row.names = FALSE)
write.csv(hed_data, file.path(output_folder, "clean_hed.csv"), row.names = TRUE)
write.csv(cibersortx_data, file.path(output_folder, "clean_cibersortx.csv"), row.names = TRUE)
write.csv(xdata, file.path(output_folder, "clean_alldata.csv"), row.names = FALSE)


# Prepare and export patient metadata -------------------------------------

patients <- data.frame(patient_id = levels(metadata$patient_id))
patients$n_samples <- tapply(X = metadata$sample_name, INDEX = metadata$patient_id, FUN = length)
patients$mean_age <- tapply(X = metadata$age, INDEX = metadata$patient_id, FUN = function(x) round(mean(x), 1))
patients$gender <- tapply(X = metadata$gender, INDEX = metadata$patient_id, FUN = function(x) unique(as.character(x)))
patients$dataset <- tapply(X = metadata$dataset, INDEX = metadata$patient_id, FUN = function(x) unique(as.character(x)))

patients <- patients[mixedorder(patients$patient_id), ]
patients <- patients[order(patients$dataset), ]
rownames(patients) <- NULL

saveRDS(patients, file.path(output_folder, "clean_patients_metadata.rds"))
write.csv(patients, file.path(output_folder, "clean_patients_metadata.csv"), row.names = FALSE)
