# Load libraries ----------------------------------------------------------
library(gtsummary)

# Prepare folder for tables -----------------------------------------------
output_folder <- "nonsync/02_tables"
if (!dir.exists(output_folder)) dir.create(output_folder)


# Load clean data ---------------------------------------------------------
metadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds")
xdata <- readRDS("nonsync/01_clean_data/clean_cibersortx.rds")

# create mock logical n_patients variable
metadata$n_patients <- TRUE
for (i in 2:nrow(metadata)) {
  metadata$n_patients[i] <- metadata$patient_id[i] != metadata$patient_id[i - 1]
}
nlevels(metadata$patient_id) == sum(metadata$n_patients) # it worked

# Summary tables ----------------------------------------------------------

## table with "response", split by dataset
cols_to_include <- c(
  "n_patients", "age", "gender", "response", "treatment",
  "biopsy_time", "enrichment_protocol"
)
tbl_summary(
  data = metadata,
  include = all_of(cols_to_include),
  by = "dataset",
  missing = "ifany",
  statistic = list(
    all_continuous() ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)",
    n_patients ~ "{n}"
  ),
) |>
  modify_header(all_stat_cols() ~ "**{level}**  \nN = {n} ({style_percent(p)}%)") |>
  bold_labels() |>
  modify_spanning_header(all_stat_cols() ~ "**Dataset**") |>
  italicize_levels() |>
  as_gt() |>
  gt::gtsave(filename = file.path(output_folder, "table_by_dataset_response.png"))


## table with "response", split by treatment
cols_to_include <- c(
  "age", "gender", "response",
  "biopsy_time", "enrichment_protocol", "dataset"
)
tbl_summary(
  data = metadata,
  include = all_of(cols_to_include),
  by = "treatment",
  missing = "ifany",
  statistic = list(
    all_continuous() ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)"
  ),
) |>
  modify_header(all_stat_cols() ~ "**{level}**  \nN = {n} ({style_percent(p)}%)") |>
  bold_labels() |>
  modify_spanning_header(all_stat_cols() ~ "**Treatment**") |>
  italicize_levels() |>
  as_gt() |>
  gt::gtsave(filename = file.path(output_folder, "table_by_treatment_response.png"))

## table split by "response"
cols_to_include <- c(
  "age", "gender", "treatment",
  "biopsy_time", "enrichment_protocol", "dataset"
)
tbl_summary(
  data = metadata,
  include = all_of(cols_to_include),
  by = "response",
  missing = "ifany",
  statistic = list(
    all_continuous() ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)"
  ),
) |>
  modify_header(all_stat_cols() ~ "**{level}**  \nN = {n} ({style_percent(p)}%)") |>
  bold_labels() |>
  modify_spanning_header(all_stat_cols() ~ "**Response**") |>
  italicize_levels() |>
  as_gt() |>
  gt::gtsave(filename = file.path(output_folder, "table_by_response.png"))

## table with "response_group", split by dataset
cols_to_include <- c(
  "n_patients", "age", "gender", "response_group", "treatment",
  "biopsy_time", "enrichment_protocol"
)
tbl_summary(
  data = metadata,
  include = all_of(cols_to_include),
  by = "dataset",
  missing = "ifany",
  statistic = list(
    all_continuous() ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)",
    n_patients ~ "{n}"
  ),
) |>
  modify_header(all_stat_cols() ~ "**{level}**  \nN = {n} ({style_percent(p)}%)") |>
  bold_labels() |>
  modify_spanning_header(all_stat_cols() ~ "**Dataset**") |>
  italicize_levels() |>
  as_gt() |>
  gt::gtsave(filename = file.path(output_folder, "table_by_dataset_response_group.png"))

## table with "response_group", split by treatment
cols_to_include <- c(
  "age", "gender", "response_group",
  "biopsy_time", "enrichment_protocol", "dataset"
)
tbl_summary(
  data = metadata,
  include = all_of(cols_to_include),
  by = "treatment",
  missing = "ifany",
  statistic = list(
    all_continuous() ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)"
  ),
) |>
  modify_header(all_stat_cols() ~ "**{level}**  \nN = {n} ({style_percent(p)}%)") |>
  bold_labels() |>
  modify_spanning_header(all_stat_cols() ~ "**Treatment**") |>
  italicize_levels() |>
  as_gt() |>
  gt::gtsave(filename = file.path(output_folder, "table_by_treatment_response_group.png"))

## table split by "response_group"
cols_to_include <- c(
  "age", "gender", "treatment",
  "biopsy_time", "enrichment_protocol", "dataset"
)
tbl_summary(
  data = metadata,
  include = all_of(cols_to_include),
  by = "response_group",
  missing = "ifany",
  statistic = list(
    all_continuous() ~ "{median} ({p25}, {p75})",
    all_categorical() ~ "{n} ({p}%)"
  ),
) |>
  modify_header(all_stat_cols() ~ "**{level}**  \nN = {n} ({style_percent(p)}%)") |>
  bold_labels() |>
  modify_spanning_header(all_stat_cols() ~ "**Response**") |>
  italicize_levels() |>
  as_gt() |>
  gt::gtsave(filename = file.path(output_folder, "table_by_response_group.png"))
