# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(paletteer) # v1.7.0
  library(corrplot) # v0.95
  library(writexl) # v1.5.4
  library(car) # v3.1-5
  source("plotPCA.R")
})


# Prepare output folders --------------------------------------------------

# output folders
outdirs <- list(
  # output folder for complete dataset
  complete = "nonsync/06_PCA_multicollinearity/complete",
  # output folder for dataset without checkmate067
  nocheckmate067 = "nonsync/06_PCA_multicollinearity/nocheckmate067"
)

# create directories
for (i in seq_along(outdirs)) {
  if (!dir.exists(outdirs[[i]])) dir.create(outdirs[[i]], recursive = TRUE)
}


# Prepare data ------------------------------------------------------------

# load complete data
metadata_complete <- readRDS("nonsync/04_clean_data/clean_metadata.rds")
xdata_complete <- readRDS("nonsync/04_clean_data/clean_cibersortx.rds")
hed_complete <- readRDS("nonsync/04_clean_data/clean_hed.rds")

# move Absolute Score to the metadata
metadata_complete$cibersortx_Absolute_Score <- xdata_complete$Absolute_Score
xdata_complete <- xdata_complete[, names(xdata_complete) != "Absolute_Score"]

# add a level "Unknown" to unknown gender
metadata_complete$gender <- addNA(metadata_complete$gender)
levels(metadata_complete$gender)[nlevels(metadata_complete$gender)] <- "Unknown"

# add a level "Unknown" to unknown response
metadata_complete$response_4levels <- addNA(metadata_complete$response_4levels)
levels(metadata_complete$response_4levels)[nlevels(metadata_complete$response_4levels)] <- "Unknown"
metadata_complete$response_3levels <- addNA(metadata_complete$response_3levels)
levels(metadata_complete$response_3levels)[nlevels(metadata_complete$response_3levels)] <- "Unknown"

# add a level "Unknown" to unknown enrichment protocol
metadata_complete$enrichment_protocol <- addNA(metadata_complete$enrichment_protocol)
levels(metadata_complete$enrichment_protocol)[nlevels(metadata_complete$enrichment_protocol)] <- "Unknown"

# sanitize cell types names
cell_types <- colnames(xdata_complete)
cell_types_sanitized <- gsub(" ", "_", cell_types)
cell_types_sanitized <- gsub("\\(|\\)", "", cell_types_sanitized)

# put all data together
alldata_complete <- cbind(metadata_complete, hed_complete, xdata_complete)

# exclude checkmate067 and create lists of data
to_exclude <- which(
  alldata_complete$dataset == "Campbell-2023" &
    alldata_complete$enrichment_protocol == "targeted-mRNA-capture"
)
alldata <- list(
  complete = alldata_complete,
  nocheckmate067 = alldata_complete[-to_exclude, ] |> droplevels()
)
metadata <- list(
  complete = metadata_complete,
  nocheckmate067 = metadata_complete[-to_exclude, ] |> droplevels()
)
xdata <- list(
  complete = xdata_complete,
  nocheckmate067 = xdata_complete[-to_exclude, ] |> droplevels()
)
heddata <- list(
  complete = hed_complete,
  nocheckmate067 = hed_complete[-to_exclude, ] |> droplevels()
)
rm(metadata_complete, hed_complete, xdata_complete, to_exclude)


# Plot settings -----------------------------------------------------------
resol <- 300
transparency_colors <- 0.8

# prepare colors for response (4 levels)
colors_response4 <- paletteer_c("grDevices::RdYlBu", nlevels(metadata[[i]]$response_4levels) - 1) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_response4) <- levels(metadata[[i]]$response_4levels)[-nlevels(metadata[[i]]$response_4levels)]
colors_response4["Unknown"] <- adjustcolor("white", transparency_colors) # white for NA

# prepare colors for response (3 levels)
colors_response3 <- paletteer_c("grDevices::RdYlBu", nlevels(metadata[[i]]$response_3levels) - 1) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_response3) <- levels(metadata[[i]]$response_3levels)[-nlevels(metadata[[i]]$response_3levels)]
colors_response3["Unknown"] <- adjustcolor("white", transparency_colors) # white for NA

# prepare colors for response (2 levels)
colors_response2 <- paletteer_d("ggsci::default_jama", nlevels(metadata[[i]]$response_2levels)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_response2) <- levels(metadata[[i]]$response_2levels)

# prepare colors for sex
colors_gender <- paletteer_d("RColorBrewer::Dark2", nlevels(metadata[[i]]$gender)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_gender) <- levels(metadata[[i]]$gender)
colors_gender["Unknown"] <- adjustcolor("white", transparency_colors) # white for NA

# prepare colors for enrichment protocol
colors_enrichment_protocol <- paletteer_d("ggsci::default_igv", nlevels(metadata[[i]]$enrichment_protocol)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_enrichment_protocol) <- levels(metadata[[i]]$enrichment_protocol)
colors_enrichment_protocol["Unknown"] <- adjustcolor("white", transparency_colors) # white for NA

# prepare colors for dataset
colors_dataset <- c("#C8DE7B", "#8B0000", "#E07A5F", "#FB6F92", "#E0BE36", "#00A0D1") |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_dataset) <- levels(metadata[[i]]$dataset)

# prepare colors for treatment
colors_treatment <- paletteer_d("ggsci::default_jco", nlevels(metadata[[i]]$treatment)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_treatment) <- levels(metadata[[i]]$treatment)


# PCA ---------------------------------------------------------------------

# transform and scale data
trdata <- lapply(xdata, function(x) {
  # log-transform data
  xx <- as.matrix(x)
  xx <- log(x + 0.001) # log transform adding a small epsilon to avoid log(0)

  # scale data to a mean of 0 and sd of 1 for each column (i.e, cell type)
  xx <- scale(xx)
  stopifnot(all(is.finite(xx))) # no "bad" values
  xx # return transformed and scaled data matrix
})

# run PCA
pca <- lapply(trdata, prcomp, center = FALSE, scale. = FALSE)

# prepare folders for PCA results
pca_dirs <- lapply(outdirs, file.path, "PCA")
for (i in seq_along(pca_dirs)) {
  if (!dir.exists(pca_dirs[[i]])) dir.create(pca_dirs[[i]])
}

# export pca results
for (i in seq_along(pca)) {
  saveRDS(pca[[i]], file.path(pca_dirs[[i]], "results_pca.rds"))
}

# plot PCs
pcs_to_plot <- paste0("PC", 1:3)
combinations <- combn(pcs_to_plot, 2) |>
  as.data.frame() |>
  as.list()
for (i in seq_along(pca)) {
  for (combination in combinations) {
    plot_pca_factor(
      pca_res = pca[[i]], # output of prcomp()
      PCs_to_plot = combination, # PCs to plot on X and Y axis respectively
      metadata = metadata[[i]][names(metadata[[i]]) != "cibersortx_Absolute_Score"], # dataframe of metadata[[i]] (rows matching PCA data)
      vars_to_use = c( # which variables (columns of metadata[[i]]) should be used?
        "response_4levels", "response_3levels", "response_2levels", "treatment", "gender", "enrichment_protocol", "dataset"
      ),
      colors_vars_list = list( # list of colors for the variables
        response_4levels = colors_response4,
        response_3levels = colors_response3,
        response_2levels = colors_response2,
        treatment = colors_treatment,
        gender = colors_gender,
        enrichment_protocol = colors_enrichment_protocol,
        dataset = colors_dataset
      ),
      par_mar = c(2.5, 2.5, 2.5, 8), # graphical parameters par("mar")
      dark_theme = FALSE, # set to true to have dark theme
      pt.cex = 1.3, # point size
      equal_axis_scale = FALSE, # set to TRUE to force same scaling of the X and Y axes
      legend_cex = 0.7, # relative dimension of the legend
      legend_ncols = 1, # columns to split the legend in
      output_folder = pca_dirs[[i]], # folder to save the plots
      file_format = "png", # file format
      res_ppi = 300, # resolution (pixels per inch)
      width_in = 6, height_in = 4 # width and height of the plot in inches
    )
  }
}

# scree plot
for (i in seq_along(pca)) {
  percentVar <- round(pca[[i]]$sdev^2 / sum(pca[[i]]$sdev^2) * 100, 1)
  names(percentVar) <- colnames(pca[[i]]$x)
  png(file.path(pca_dirs[[i]], "scree_plot.png"),
    width = 5 * resol, height = 4 * resol, res = resol
  )
  par(mar = c(3.5, 3, 0.5, 0.1), mgp = c(2, 0.8, 0), tcl = -0.3)
  barplot(percentVar, las = 2, ylab = "% of total variance explained")
  dev.off()
}


# PCA loadings ------------------------------------------------------------

# export table of loadings
for (i in seq_along(pca)) {
  write.csv(
    as.data.frame(round(pca[[i]]$rotation, 3)),
    file.path(outdirs[[i]], "loadings_pca.csv")
  )
}

# function to plot PCA loadings as arrows
arrows_pca_loadings <- function(
  PCx = "PC1", PCy = "PC2",
  text.cex = 0.5, col.text = "black", loadings_data, ...
) {
  rot_x <- loadings_data[, PCx]
  rot_y <- loadings_data[, PCy]
  arrows(
    x0 = rep(0, length(rot_x)), y0 = rep(0, length(rot_y)),
    x1 = rot_x, y1 = rot_y, ...
  )
  angles <- atan2(rot_y, rot_x)
  angles <- ifelse(angles < 0, angles + 2 * pi, angles)
  position_text <- ifelse(
    angles >= pi / 4 & angles <= pi * (3 / 4), 3, ifelse(
      angles >= pi * (3 / 4) & angles <= pi * (5 / 4), 2, ifelse(
        angles >= pi * (5 / 4) & angles <= pi * (7 / 4), 1, 4
      )
    )
  )
  text(
    x = rot_x, y = rot_y, labels = rownames(loadings_data),
    pos = position_text, cex = text.cex, col = col.text
  )
}

# function to plot rotation
plot_rotation <- function(PCx = "PC1", PCy = "PC2", loadings_data) {
  par(las = 1, mar = rep(3, 4), mgp = c(2, 0.7, 0), tcl = -0.3, xpd = TRUE)
  plot(NULL,
    xlim = c(-1, 1), ylim = c(-1, 1),
    xlab = PCx, ylab = PCy,
    xaxs = "i", yaxs = "i", asp = 1, bty = "l",
    main = paste0("Rotation ", PCy, " ~ ", PCx),
  )
  colors_celltypes <- paletteer_d("ggsci::default_igv", nrow(loadings_data))
  arrows_pca_loadings(
    PCx = PCx, PCy = PCy, loadings_data = loadings_data,
    text.cex = 0.7, col.text = colors_celltypes, lwd = 1, length = 0.03, col = colors_celltypes
  )
}

# plot rotation with arrows for PCA loadings
for (i in seq_along(pca)) {
  for (combination in combinations) {
    png(
      file.path(outdirs[[i]], paste0(
        "PCA_rotation_", combination[1], "_", combination[2], ".png"
      )),
      width = 5 * resol, height = 5 * resol, res = resol
    )
    plot_rotation(
      PCx = combination[1], PCy = combination[2], loadings_data = pca[[i]]$rotation
    )
    dev.off()
  }
}


# VIF ---------------------------------------------------------------------

for (i in seq_along(alldata)) {
  # prepare data for modelling (only binary response and CIBERSORTx data)
  vifdata <- cbind(response = alldata[[i]]$response_2levels, xdata[[i]])
  names(vifdata) <- c("response", cell_types_sanitized)
  vifdata$response <- as.numeric(vifdata$response)

  # fit linear model
  xformula <- paste0("response ~ ", paste(
    cell_types_sanitized,
    collapse = " + "
  )) |> as.formula()
  xml <- lm(xformula, data = vifdata)

  # evaluate VIFs
  vifs <- vif(xml)
  range(vifs) |> round(2)

  # save VIFs as txt
  sink(file.path(outdirs[[i]], "vifs_cibersortx.txt"))
  print(round(vifs, 3))
  sink()
}


# Correlation between all numeric variables -------------------------------

for (i in seq_along(alldata)) {
  numeric_variables <- names(alldata[[i]])[sapply(alldata[[i]], is.numeric)]
  cormat_all <- cor(alldata[[i]][, numeric_variables],
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  # make correlation plot
  png(file.path(outdirs[[i]], "corrplot_all_numeric.png"),
    width = 8.5 * resol, height = 6.5 * resol, res = resol
  )
  corrplot(
    cormat_all,
    order = "original",
    col = paletteer_c("grDevices::RdBu", 150, direction = -1),
    tl.col = "black", tl.cex = 0.8, tl.srt = 45,
    method = "circle",
    type = "full",
    is.corr = TRUE,
    title = "Spearman correlation matrix",
    diag = TRUE,
    outline = FALSE,
    mar = c(0, 0, 1.5, 0)
  )
  dev.off()

  # export correlation data
  cormat_triangle <- cormat_all
  cormat_triangle[upper.tri(cormat_triangle, diag = TRUE)] <-
    NA # assign NA to upper triangle, including the diagonal
  cor_df <- as.data.frame.table(cormat_triangle)
  names(cor_df) <- c("Var1", "Var2", "spearman_correlation")
  cor_df <- na.exclude(cor_df)
  cor_df <- cor_df[order(abs(cor_df$spearman_correlation), decreasing = TRUE), ]
  cor_df$spearman_correlation <- round(cor_df$spearman_correlation, 3)
  write_xlsx(cor_df, file.path(outdirs[[i]], "correlations_cibersortx_data.xlsx"))

  # create scatterplots of strongly correlated data
  xx_dir <- file.path(outdirs[[i]], "scatterplots_correlations")
  if (!dir.exists(xx_dir)) dir.create(xx_dir)
  strong_cor_df <- subset(cor_df, abs(spearman_correlation) > 0.5)
  for (j in seq_len(nrow(strong_cor_df))) {
    x_cell <- strong_cor_df$Var1[j] |> as.character()
    y_cell <- strong_cor_df$Var2[j] |> as.character()
    png(file.path(xx_dir, paste0(
      "scatter.", j, ".", x_cell, "-", y_cell, ".png"
    )), height = 4.5 * resol, width = 4 * resol, res = resol)
    plot(
      x = rank(alldata[[i]][, x_cell]), y = rank(alldata[[i]][, y_cell]),
      las = 1, xlab = paste(x_cell, "(Rank)"), ylab = paste(y_cell, "(Rank)"), asp = 1,
      main = paste0("Spearman correlation = ", strong_cor_df$spearman_correlation[j]),
      pch = 21, bg = colors_dataset[as.numeric(metadata[[i]]$dataset)]
    )
    dev.off()
  }
}


# Correlation between cibersortx columns ----------------------------------

for (i in seq_along(alldata)) {
  xx <- c(names(xdata[[i]]), "cibersortx_Absolute_Score")
  cormat_cibersortx <- cor(alldata[[i]][, xx],
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  png(file.path(outdirs[[i]], "corrplot_cibersortx.png"),
    width = 8 * resol, height = 6 * resol, res = resol
  )
  corrplot(
    cormat_cibersortx,
    order = "original",
    col = paletteer_c("grDevices::RdBu", 150, direction = -1),
    tl.col = "black", tl.cex = 0.8, tl.srt = 45,
    method = "circle",
    type = "full",
    is.corr = TRUE,
    title = "Spearman correlation matrix",
    diag = TRUE,
    outline = FALSE,
    mar = c(0, 0, 1.5, 0)
  )
  dev.off()
}
