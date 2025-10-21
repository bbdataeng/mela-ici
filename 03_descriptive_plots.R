# Load libraries ----------------------------------------------------------
library(paletteer)
source("barplot_functions.R") # source function for making stacked barcharts

# Prepare folder for figures ----------------------------------------------
output_folder <- "nonsync/03_descriptive_plots"
if (!dir.exists(output_folder)) dir.create(output_folder)

# Load clean data ---------------------------------------------------------
metadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds")
xdata <- readRDS("nonsync/01_clean_data/clean_cibersortx.rds")


# Prepare plots -----------------------------------------------------------

# prepare colors for response
colors_response <- paletteer_d("ggsci::default_nejm", nlevels(metadata$response))
names(colors_response) <- levels(metadata$response)

# prepare colors for response_group
colors_response_group <- paletteer_d("ggsci::default_nejm", nlevels(metadata$response_group))
names(colors_response_group) <- levels(metadata$response_group)

# prepare colors for gender
colors_gender <- paletteer_d("RColorBrewer::Dark2", nlevels(metadata$gender))
names(colors_gender) <- levels(metadata$gender)

# prepare colors for treatment
colors_treatment <- paletteer_d("ggsci::default_jco", nlevels(metadata$treatment))
names(colors_treatment) <- levels(metadata$treatment)

# plot resolution
resol <- 300


# Prepare subfolders ------------------------------------------------------
response_folder <- file.path(output_folder, "response")
response_group_folder <- file.path(output_folder, "response_group")
gender_folder <- file.path(output_folder, "gender")
age_folder <- file.path(output_folder, "age")

for (xx in c(response_folder, response_group_folder, gender_folder, age_folder)) {
  if (!dir.exists(xx)) dir.create(xx)
}



# Plots response ~ treatment ----------------------------------------------

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = metadata$treatment, # x-axis variable
  x_var_name = "Treatment", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_treatment_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 6, # width and height of the plot (in inches)
  mar = c(8, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = metadata$treatment, # x-axis variable
  x_var_name = "Treatment", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_treatment_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 6, # width and height of the plot (in inches)
  mar = c(8, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# Plots response_group ~ treatment ----------------------------------------

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = metadata$treatment, # x-axis variable
  x_var_name = "Treatment", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_treatment_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 6, # width and height of the plot (in inches)
  mar = c(8, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = metadata$treatment, # x-axis variable
  x_var_name = "Treatment", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_treatment_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 6, # width and height of the plot (in inches)
  mar = c(8, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response ~ biopsy_time --------------------------------------------

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = metadata$biopsy_time, # x-axis variable
  x_var_name = "Biopsy Time", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = FALSE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_biopsy-time_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 5, # width and height of the plot (in inches)
  mar = c(3.5, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = metadata$biopsy_time, # x-axis variable
  x_var_name = "Biopsy Time", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = FALSE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_biopsy-time_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 5, # width and height of the plot (in inches)
  mar = c(3.5, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response_group ~ biopsy_time --------------------------------------

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = metadata$biopsy_time, # x-axis variable
  x_var_name = "Biopsy Time", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = FALSE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_biopsy-time_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 5, # width and height of the plot (in inches)
  mar = c(3.5, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = metadata$biopsy_time, # x-axis variable
  x_var_name = "Biopsy Time", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = FALSE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_biopsy-time_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 5, # width and height of the plot (in inches)
  mar = c(3.5, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response ~ gender -------------------------------------------------

gender_with_NA <- addNA(metadata$gender)
levels(gender_with_NA)
levels(gender_with_NA)[nlevels(gender_with_NA)] <- "Unknown"

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = gender_with_NA, # x-axis variable
  x_var_name = "Gender", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = FALSE, # don't rotate x_axis labels
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_gender_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 5, # width and height of the plot (in inches)
  mar = c(4, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = gender_with_NA, # x-axis variable
  x_var_name = "Gender", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = FALSE, # don't rotate x_axis labels
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_gender_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 5, # width and height of the plot (in inches)
  mar = c(4, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response_group ~ gender -------------------------------------------

gender_with_NA <- addNA(metadata$gender)
levels(gender_with_NA)
levels(gender_with_NA)[nlevels(gender_with_NA)] <- "Unknown"

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = gender_with_NA, # x-axis variable
  x_var_name = "Gender", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = FALSE, # don't rotate x_axis labels
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_gender_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 5, # width and height of the plot (in inches)
  mar = c(4, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = gender_with_NA, # x-axis variable
  x_var_name = "Gender", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = FALSE, # don't rotate x_axis labels
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_gender_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 4, height_in = 5, # width and height of the plot (in inches)
  mar = c(4, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response ~ enrichment_protocol ------------------------------------

enrich_prot_with_NA <- addNA(metadata$enrichment_protocol)
levels(enrich_prot_with_NA)[nlevels(enrich_prot_with_NA)] <- "Unknown"

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = enrich_prot_with_NA, # x-axis variable
  x_var_name = "Enrichment Protocol", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_enrichment-protocol_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 5, height_in = 5, # width and height of the plot (in inches)
  mar = c(8, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = enrich_prot_with_NA, # x-axis variable
  x_var_name = "Enrichment Protocol", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_enrichment-protocol_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 5, height_in = 5, # width and height of the plot (in inches)
  mar = c(8, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response_group ~ enrichment_protocol ------------------------------

enrich_prot_with_NA <- addNA(metadata$enrichment_protocol)
levels(enrich_prot_with_NA)[nlevels(enrich_prot_with_NA)] <- "Unknown"

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = enrich_prot_with_NA, # x-axis variable
  x_var_name = "Enrichment Protocol", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_enrichment-protocol_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 5, height_in = 5, # width and height of the plot (in inches)
  mar = c(8, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = enrich_prot_with_NA, # x-axis variable
  x_var_name = "Enrichment Protocol", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_enrichment-protocol_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 5, height_in = 5, # width and height of the plot (in inches)
  mar = c(8, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response ~ dataset ------------------------------------------------

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = metadata$dataset, # x-axis variable
  x_var_name = "Dataset", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_dataset_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 5, height_in = 5, # width and height of the plot (in inches)
  mar = c(6, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = metadata$dataset, # x-axis variable
  x_var_name = "Dataset", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_dataset_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 5, height_in = 5, # width and height of the plot (in inches)
  mar = c(6, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response_group ~ dataset ------------------------------------------

# stacked barchart showing absolute frequencies
make_stacked_barplot(
  x_var = metadata$dataset, # x-axis variable
  x_var_name = "Dataset", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_dataset_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 5, height_in = 5, # width and height of the plot (in inches)
  mar = c(6, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = metadata$dataset, # x-axis variable
  x_var_name = "Dataset", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_dataset_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 5, height_in = 5, # width and height of the plot (in inches)
  mar = c(6, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response ~ age ----------------------------------------------------

# dotplot with continuous age
xx <- table(Response = metadata$response, Age = metadata$age) |>
  as.data.frame()
xx$x <- as.numeric(as.character(xx$Age))
xx$y <- as.numeric(xx$Response)
expansion_factor <- 5
xx$cex <- (xx$Freq / max(xx$Freq)) |> sqrt() * expansion_factor
png(file.path(response_folder, "response_age_dotplot.png"),
  width = 6 * resol, height = 4 * resol, res = resol
)
par(mar = c(3.5, 5, 0.2, 6.5), mgp = c(3.5, 0.8, 0), tcl = -0.3, xpd = FALSE)
plot(NULL,
  xlim = range(metadata$age, na.rm = TRUE),
  ylim = c(nlevels(metadata$response), 1) + c(0.5, -0.5),
  yaxs = "i", bty = "l",
  xlab = "", ylab = "Response", yaxt = "n"
)
mtext(text = "Age", side = 1, line = par("mar")[1] - 1.5)
axis(
  side = 2, at = 1:nlevels(metadata$response),
  labels = levels(metadata$response), las = 1
)
grid()
abline(h = 1:nlevels(metadata$response), col = "gray70", lty = "dotted")
points(x = xx$x, y = xx$y, cex = xx$cex, pch = 16, col = adjustcolor("#023E8A", 0.8))
par(xpd = TRUE)
legend_labels <- c(0, 1, 2, 4, 6, 8, 10)
legend_pt.cx <- (legend_labels / max(xx$Freq)) |> sqrt() * expansion_factor
legend(
  x = par("usr")[2] + 2,
  y = mean(par("usr")[3:4]),
  xjust = 0, yjust = 0.5, bty = "n", cex = 0.8,
  legend = legend_labels, pt.cex = legend_pt.cx, col = adjustcolor("#023E8A", 0.8),
  pch = 16, y.intersp = 2, x.intersp = 2.5, title = "N. observations"
)
dev.off()

# stacked barchart showing absolute frequencies
bin_age <- cut(metadata$age, breaks = seq(0, 100, 5)) |>
  droplevels() |>
  addNA()
levels(bin_age)[nlevels(bin_age)] <- "Unknown"
make_stacked_barplot(
  x_var = bin_age, # x-axis variable
  x_var_name = "Age (binned)", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_aged_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 8, height_in = 6, # width and height of the plot (in inches)
  mar = c(4.5, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = bin_age, # x-axis variable
  x_var_name = "Age (binned)", # name of the x-axis variable
  col_var = metadata$response, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_folder, "response_age_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 8, height_in = 6, # width and height of the plot (in inches)
  mar = c(4.5, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots response_group ~ age ----------------------------------------------

# dotplot with continuous age
xx <- table(Response = metadata$response_group, Age = metadata$age) |>
  as.data.frame()
xx$x <- as.numeric(as.character(xx$Age))
xx$y <- as.numeric(xx$Response)
expansion_factor <- 5
xx$cex <- (xx$Freq / max(xx$Freq)) |> sqrt() * expansion_factor
png(file.path(response_group_folder, "response_age_dotplot.png"),
  width = 6 * resol, height = 3 * resol, res = resol
)
par(mar = c(3.5, 5, 0.2, 6.5), mgp = c(3.5, 0.8, 0), tcl = -0.3, xpd = FALSE)
plot(NULL,
  xlim = range(metadata$age, na.rm = TRUE),
  ylim = c(nlevels(metadata$response_group), 1) + c(0.5, -0.5),
  yaxs = "i", bty = "l",
  xlab = "", ylab = "Response", yaxt = "n"
)
mtext(text = "Age", side = 1, line = par("mar")[1] - 1.5)
axis(
  side = 2, at = 1:nlevels(metadata$response_group),
  labels = levels(metadata$response_group), las = 1
)
grid()
abline(h = 1:nlevels(metadata$response_group), col = "gray70", lty = "dotted")
points(x = xx$x, y = xx$y, cex = xx$cex, pch = 16, col = adjustcolor("#023E8A", 0.8))
par(xpd = TRUE)
legend_labels <- c(0, 1, 2, 4, 6, 8, 10)
legend_pt.cx <- (legend_labels / max(xx$Freq)) |> sqrt() * expansion_factor
legend(
  x = par("usr")[2] + 2,
  y = mean(par("usr")[3:4]),
  xjust = 0, yjust = 0.5, bty = "n", cex = 0.8,
  legend = legend_labels, pt.cex = legend_pt.cx, col = adjustcolor("#023E8A", 0.8),
  pch = 16, y.intersp = 2, x.intersp = 2.5, title = "N. observations"
)
dev.off()

# stacked barchart showing absolute frequencies
bin_age <- cut(metadata$age, breaks = seq(0, 100, 5)) |>
  droplevels() |>
  addNA()
levels(bin_age)[nlevels(bin_age)] <- "Unknown"
make_stacked_barplot(
  x_var = bin_age, # x-axis variable
  x_var_name = "Age (binned)", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "frequency",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_aged_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 8, height_in = 6, # width and height of the plot (in inches)
  mar = c(4.5, 4, 1, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = bin_age, # x-axis variable
  x_var_name = "Age (binned)", # name of the x-axis variable
  col_var = metadata$response_group, # color-axis variable
  col_var_name = "Response", # name of the color-axis variable
  type = "proportion",
  color_palette = colors_response_group, # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(response_group_folder, "response_age_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 8, height_in = 6, # width and height of the plot (in inches)
  mar = c(4.5, 4, 4, 5), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)



# Plots gender ~ age ------------------------------------------------------

# dotplot with continuous age
xx <- table(Gender = gender_with_NA, Age = metadata$age) |>
  as.data.frame()
xx$x <- as.numeric(as.character(xx$Age))
xx$y <- as.numeric(xx$Gender)
expansion_factor <- 5
xx$cex <- (xx$Freq / max(xx$Freq)) |> sqrt() * expansion_factor
png(file.path(gender_folder, "gender_age_dotplot.png"),
  width = 6 * resol, height = 3.5 * resol, res = resol
)
par(mar = c(3.5, 5, 0.2, 6.5), mgp = c(3.5, 0.8, 0), tcl = -0.3, xpd = FALSE)
plot(NULL,
  xlim = range(metadata$age, na.rm = TRUE),
  ylim = c(nlevels(gender_with_NA), 1) + c(0.5, -0.5),
  yaxs = "i", bty = "l",
  xlab = "", ylab = "Gender", yaxt = "n"
)
mtext(text = "Age", side = 1, line = par("mar")[1] - 1.5)
axis(
  side = 2, at = 1:nlevels(gender_with_NA),
  labels = levels(gender_with_NA), las = 1
)
abline(h = 1:nlevels(gender_with_NA), col = "gray70", lty = "dotted")
abline(v = seq(0, 100, 10), col = "gray70", lty = "dotted")
points(x = xx$x, y = xx$y, cex = xx$cex, pch = 16, col = adjustcolor("#023E8A", 0.8))
par(xpd = TRUE)
legend_labels <- c(0, 1, 2, 4, 6, 8, 10)
legend_pt.cx <- (legend_labels / max(xx$Freq)) |> sqrt() * expansion_factor
legend(
  x = par("usr")[2] + 2,
  y = mean(par("usr")[3:4]),
  xjust = 0, yjust = 0.5, bty = "n", cex = 0.8,
  legend = legend_labels, pt.cex = legend_pt.cx, col = adjustcolor("#023E8A", 0.8),
  pch = 16, y.intersp = 2, x.intersp = 2.5, title = "N. observations"
)
dev.off()

# stacked barchart showing absolute frequencies
bin_age <- cut(metadata$age, breaks = seq(0, 100, 5)) |>
  droplevels() |>
  addNA()
levels(bin_age)[nlevels(bin_age)] <- "Unknown"
make_stacked_barplot(
  x_var = bin_age, # x-axis variable
  x_var_name = "Age (binned)", # name of the x-axis variable
  col_var = gender_with_NA, # color-axis variable
  col_var_name = "Gender", # name of the color-axis variable
  type = "frequency",
  color_palette = c(colors_gender, "grey"), # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = FALSE, # don't show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(gender_folder, "gender_aged_stacked_frequencies.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 8, height_in = 6, # width and height of the plot (in inches)
  mar = c(4.5, 4, 1, 6), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)

# stacked barchart showing proportions
make_stacked_barplot(
  x_var = bin_age, # x-axis variable
  x_var_name = "Age (binned)", # name of the x-axis variable
  col_var = gender_with_NA, # color-axis variable
  col_var_name = "Gender", # name of the color-axis variable
  type = "proportion",
  color_palette = c(colors_gender, "grey"), # vector of colors matching levels of col_var
  rotate_x_var_labels = TRUE, # rotate x_axis labels (useful for long names)
  show_x_var_frequencies = TRUE, # show x-var frequencies
  show_percentages = TRUE, # whether to show percentages in the bars
  min_percentage = 2, # hide percentages below this threshold
  file_path = file.path(gender_folder, "gender_age_stacked_proportions.png"), # file path (with format) to save the plot
  res_ppi = resol, # resolution of the plot
  width_in = 8, height_in = 6, # width and height of the plot (in inches)
  mar = c(4.5, 4, 4, 6), mgp = c(2.5, 0.8, 0), tcl = -0.3, las = 1 # additional par() parameters
)


# Plots age ~ response ----------------------------------------------------

# age ~ response + gender + dataset
age_NAs_as_0 <- metadata$age
set.seed(1)
age_NAs_as_0[is.na(age_NAs_as_0)] <- rep(0, sum(is.na(age_NAs_as_0))) |> jitter(amount = 2)
png(file.path(age_folder, "age_response_gender_dataset_boxplot.png"),
  width = 6.5 * resol, height = 5 * resol, res = resol
)
par(mar = c(3, 3, 1, 8), mgp = c(2, 0.8, 0), tcl = -0.3, xpd = FALSE, las = 1)
plot(NULL,
  xlim = c(0, nlevels(metadata$response)) + 0.5, xaxt = "n",
  ylim = c(0, max(metadata$age, na.rm = TRUE)), yaxt = "n",
  xlab = "Response", ylab = "Age", bty = "n",
)
axis(side = 1, at = 1:nlevels(metadata$response), labels = levels(metadata$response))
axis(side = 2, at = seq(10, 90, 10))
axis(side = 2, at = 0, labels = "NA")
grid()
par(xpd = TRUE)
rect(
  xleft = par("usr")[1], xright = par("usr")[2],
  ybottom = par("usr")[3], ytop = -par("usr")[3],
  col = adjustcolor("red", 0.1)
)
boxplot(age ~ response,
  data = metadata, add = TRUE, col = "white", lwd = 1.5,
  outline = FALSE, axes = FALSE
)
set.seed(12345)
points(
  x = as.numeric(metadata$response) |> jitter(amount = 0.2),
  y = age_NAs_as_0,
  pch = 20 + as.numeric(metadata$dataset),
  cex = 1.2, lwd = 1.5,
  col = adjustcolor(c(colors_gender, "grey"), .8)[as.numeric(gender_with_NA)]
)
bottom_coord <- grconvertY(0, from = "ndc", to = "user")
top_coord <- grconvertY(1, from = "ndc", to = "user")
right_coord <- grconvertX(1, from = "ndc", to = "user")
legend(
  x = mean(c(par("usr")[2], right_coord)),
  y = mean(c(bottom_coord, top_coord)),
  xjust = 0.5, yjust = 0, bty = "n",
  title = "Gender",
  legend = levels(gender_with_NA),
  fill = adjustcolor(c(colors_gender, "grey"), .8)
)
legend(
  x = mean(c(par("usr")[2], right_coord)),
  y = mean(c(bottom_coord, top_coord)),
  xjust = 0.5, yjust = 1, bty = "n",
  title = "Dataset",
  legend = levels(metadata$dataset),
  pch = 20 + 1:nlevels(metadata$dataset),
  pt.cex = 1.2, pt.lwd = 1.5
)
dev.off()

# age ~ response + treatment + dataset
age_NAs_as_0 <- metadata$age
set.seed(1)
age_NAs_as_0[is.na(age_NAs_as_0)] <- rep(0, sum(is.na(age_NAs_as_0))) |> jitter(amount = 2)
png(file.path(age_folder, "age_response_treatment_dataset_boxplot.png"),
  width = 6.5 * resol, height = 5 * resol, res = resol
)
par(mar = c(3, 3, 1, 9), mgp = c(2, 0.8, 0), tcl = -0.3, xpd = FALSE, las = 1)
plot(NULL,
  xlim = c(0, nlevels(metadata$response)) + 0.5, xaxt = "n",
  ylim = c(0, max(metadata$age, na.rm = TRUE)), yaxt = "n",
  xlab = "Response", ylab = "Age", bty = "n",
)
axis(side = 1, at = 1:nlevels(metadata$response), labels = levels(metadata$response))
axis(side = 2, at = seq(10, 90, 10))
axis(side = 2, at = 0, labels = "NA")
grid()
par(xpd = TRUE)
rect(
  xleft = par("usr")[1], xright = par("usr")[2],
  ybottom = par("usr")[3], ytop = -par("usr")[3],
  col = adjustcolor("red", 0.1)
)
boxplot(age ~ response,
  data = metadata, add = TRUE, col = "white", lwd = 1.5,
  outline = FALSE, axes = FALSE
)
set.seed(12345)
points(
  x = as.numeric(metadata$response) |> jitter(amount = 0.2),
  y = age_NAs_as_0,
  pch = 20 + as.numeric(metadata$dataset),
  cex = 1.2, lwd = 1.5,
  col = adjustcolor(colors_treatment, .8)[as.numeric(metadata$treatment)]
)
bottom_coord <- grconvertY(0, from = "ndc", to = "user")
top_coord <- grconvertY(1, from = "ndc", to = "user")
right_coord <- grconvertX(1, from = "ndc", to = "user")
legend(
  x = mean(c(par("usr")[2], right_coord)),
  y = mean(c(bottom_coord, top_coord)),
  xjust = 0.5, yjust = 0, bty = "n", cex = 0.8, title.cex = 1,
  title = "Treatment",
  legend = levels(metadata$treatment),
  fill = adjustcolor(colors_treatment, .8)
)
legend(
  x = mean(c(par("usr")[2], right_coord)),
  y = mean(c(bottom_coord, top_coord)),
  xjust = 0.5, yjust = 1, bty = "n",
  title = "Dataset",
  legend = levels(metadata$dataset),
  pch = 20 + 1:nlevels(metadata$dataset),
  pt.cex = 1.2, pt.lwd = 1.5
)
dev.off()


# Plots age ~ dataset -----------------------------------------------------
age_NAs_as_0 <- metadata$age
set.seed(1)
age_NAs_as_0[is.na(age_NAs_as_0)] <- rep(0, sum(is.na(age_NAs_as_0))) |> jitter(amount = 2)
png(file.path(age_folder, "age_dataset_response_gender_boxplot.png"),
  width = 6.5 * resol, height = 5 * resol, res = resol
)
par(mar = c(6, 3, 1, 8), mgp = c(2, 0.8, 0), tcl = -0.3, xpd = FALSE, las = 1)
plot(NULL,
  xlim = c(0, nlevels(metadata$dataset)) + 0.5, xaxt = "n",
  ylim = c(0, max(metadata$age, na.rm = TRUE)), yaxt = "n",
  xlab = "", ylab = "Age", bty = "n",
)
axis(side = 1, at = 1:nlevels(metadata$dataset), labels = FALSE)
axis(side = 2, at = seq(10, 90, 10))
axis(side = 2, at = 0, labels = "NA")
grid()
par(xpd = TRUE)
text(
  labels = levels(metadata$dataset), x = 1:nlevels(metadata$dataset),
  y = par()$usr[3] - (par()$usr[4] - par()$usr[3]) * 0.035,
  srt = 45, cex = 1, adj = 1
)
mtext(side = 1, text = "Dataset", line = par("mar")[1] - 1)
rect(
  xleft = par("usr")[1], xright = par("usr")[2],
  ybottom = par("usr")[3], ytop = -par("usr")[3],
  col = adjustcolor("red", 0.1)
)
boxplot(age ~ dataset,
  data = metadata, add = TRUE, col = "white", lwd = 1.5,
  outline = FALSE, axes = FALSE
)
set.seed(12345)
points(
  x = as.numeric(metadata$dataset) |> jitter(amount = 0.2),
  y = age_NAs_as_0,
  pch = 20 + as.numeric(gender_with_NA),
  cex = 1.2, lwd = 1.5,
  col = adjustcolor(colors_response, .8)[as.numeric(metadata$response)]
)
bottom_coord <- grconvertY(0, from = "ndc", to = "user")
top_coord <- grconvertY(1, from = "ndc", to = "user")
right_coord <- grconvertX(1, from = "ndc", to = "user")
legend(
  x = mean(c(par("usr")[2], right_coord)),
  y = mean(c(bottom_coord, top_coord)),
  xjust = 0.5, yjust = 0, bty = "n",
  title = "Response",
  legend = levels(metadata$response),
  fill = adjustcolor(colors_response, .8)
)
legend(
  x = mean(c(par("usr")[2], right_coord)),
  y = mean(c(bottom_coord, top_coord)),
  xjust = 0.5, yjust = 1, bty = "n",
  title = "Gender",
  legend = levels(gender_with_NA),
  pch = 20 + 1:nlevels(gender_with_NA),
  pt.cex = 1.2, pt.lwd = 1.5
)
dev.off()
