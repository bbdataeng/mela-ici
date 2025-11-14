# Load libraries ----------------------------------------------------------
library(gtsummary)

# Prepare folder for tables -----------------------------------------------
output_folder <- "nonsync/02_tables"
if (!dir.exists(output_folder)) dir.create(output_folder)


# Load clean data ---------------------------------------------------------
metadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds")
xdata <- readRDS("nonsync/01_clean_data/clean_cibersortx.rds")

# Summary tables ----------------------------------------------------------

# for each version of response variable
for (response in grepv("^response", names(metadata))) {
  ## table split by dataset
  cols_to_include <- c(
    "age", "gender", response, "treatment", "enrichment_protocol"
  )
  tbl_summary(
    data = metadata,
    include = all_of(cols_to_include),
    by = "dataset",
    missing = "ifany",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
  ) |>
    modify_header(all_stat_cols() ~ "**{level}**  \nN = {n} ({style_percent(p)}%)") |>
    bold_labels() |>
    modify_spanning_header(all_stat_cols() ~ "**Dataset**") |>
    italicize_levels() |>
    as_gt() |>
    gt::gtsave(filename = file.path(output_folder, paste0(
      "table_by_dataset_", response, ".png"
    )))


  ## table split by treatment
  cols_to_include <- c(
    "age", "gender", response, "enrichment_protocol", "dataset"
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
    gt::gtsave(filename = file.path(output_folder, paste0(
      "table_by_treatment_", response, ".png"
    )))

  ## table split by response
  cols_to_include <- c(
    "age", "gender", "treatment", "enrichment_protocol", "dataset"
  )
  tbl_summary(
    data = metadata,
    include = all_of(cols_to_include),
    by = response,
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
    gt::gtsave(filename = file.path(output_folder, paste0(
      "table_by_", response, ".png"
    )))
}
