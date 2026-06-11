# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(gtsummary) # v2.5.0
  library(stringr) # v1.6.0
})

# Prepare folder for tables -----------------------------------------------

# output folder for complete dataset
outdir_complete <- "nonsync/05_tables/complete"
if (!dir.exists(outdir_complete)) dir.create(outdir_complete, recursive = TRUE)

# output folder for dataset without checkmate067
outdir_nocheckmate067 <- "nonsync/05_tables/nocheckmate067"
if (!dir.exists(outdir_nocheckmate067)) dir.create(outdir_nocheckmate067)


# Prepare data ------------------------------------------------------------

# complete data
metadata_complete <- readRDS("nonsync/04_clean_data/clean_metadata.rds")

# exclude checkmate067
metadata_subset <- subset(
  metadata_complete,
  !(dataset == "Campbell-2023" & enrichment_protocol == "targeted-mRNA-capture")
) |> droplevels()


# Function to make tables -------------------------------------------------

make_table <- function(data, cols_to_include, by_var, path) {
  tbl_summary(
    data = data,
    include = all_of(cols_to_include),
    by = by_var,
    missing = "ifany",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
  ) |>
    modify_header(all_stat_cols() ~ "**{level}**  \nN = {n} ({style_percent(p)}%)") |>
    bold_labels() |>
    modify_spanning_header(all_stat_cols() ~ paste0("**", str_to_sentence(by_var), "**")) |>
    italicize_levels() |>
    as_gt() |>
    gt::gtsave(filename = path)
}


# Summary tables ----------------------------------------------------------

# for each version of response variable
for (response in grepv("^response", names(metadata_complete))) {
  ## table split by dataset
  cols_to_include <- c(
    "age", "gender", response, "treatment", "enrichment_protocol"
  )
  make_table(metadata_complete, cols_to_include, "dataset", file.path(
    outdir_complete, paste0("table_by_dataset_", response, ".docx")
  ))
  make_table(metadata_subset, cols_to_include, "dataset", file.path(
    outdir_nocheckmate067, paste0("table_by_dataset_", response, ".docx")
  ))

  ## table split by treatment
  cols_to_include <- c(
    "age", "gender", response, "enrichment_protocol", "dataset"
  )
  make_table(metadata_complete, cols_to_include, "treatment", file.path(
    outdir_complete, paste0("table_by_treatment_", response, ".docx")
  ))
  make_table(metadata_subset, cols_to_include, "treatment", file.path(
    outdir_nocheckmate067, paste0("table_by_treatment_", response, ".docx")
  ))

  ## table split by response
  cols_to_include <- c(
    "age", "gender", "treatment", "enrichment_protocol", "dataset"
  )
  make_table(metadata_complete, cols_to_include, response, file.path(
    outdir_complete, paste0("table_by_", response, ".docx")
  ))
  make_table(metadata_subset, cols_to_include, response, file.path(
    outdir_nocheckmate067, paste0("table_by_", response, ".docx")
  ))
}
