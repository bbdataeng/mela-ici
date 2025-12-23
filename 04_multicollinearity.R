# Load libraries ----------------------------------------------------------
library(paletteer)
library(corrplot)
library(writexl)
source("plotPCA.R")

# Prepare folder for figures ----------------------------------------------
output_folder <- "nonsync/04_multicollinearity"
if (!dir.exists(output_folder)) dir.create(output_folder)

# Load clean data ---------------------------------------------------------
metadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds") |> as.data.frame()
xdata <- readRDS("nonsync/01_clean_data/clean_cibersortx.rds") |> as.data.frame()
hed_data <- readRDS("nonsync/01_clean_data/clean_hed.rds") |> as.data.frame()

# add a level "Unknown" to unknown gender
metadata$gender <- addNA(metadata$gender)
levels(metadata$gender)[nlevels(metadata$gender)] <- "Unknown"

# add a level "Unknown" to unknown enrichment protocol
metadata$enrichment_protocol <- addNA(metadata$enrichment_protocol)
levels(metadata$enrichment_protocol)[nlevels(metadata$enrichment_protocol)] <- "Unknown"

# sanitize cell types names
cell_types <- colnames(xdata)
cell_types_sanitized <- gsub(" ", "_", cell_types)
cell_types_sanitized <- gsub("\\(|\\)", "", cell_types_sanitized)

# put all data together
alldata <- cbind(metadata, hed_data, xdata)


# Plot settings -----------------------------------------------------------
resol <- 300
transparency_colors <- 0.8

# prepare colors for response (6 levels)
colors_response6 <- paletteer_c("grDevices::RdYlBu", nlevels(metadata$response_6levels)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_response6) <- levels(metadata$response_6levels)

# prepare colors for response (3 levels)
colors_response3 <- paletteer_c("grDevices::RdYlBu", nlevels(metadata$response_3levels)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_response3) <- levels(metadata$response_3levels)

# prepare colors for response (2 levels)
colors_response2 <- paletteer_d("ggsci::default_jama", nlevels(metadata$response_2levels)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_response2) <- levels(metadata$response_2levels)

# prepare colors for sex
colors_gender <- paletteer_d("RColorBrewer::Dark2", nlevels(metadata$gender)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_gender) <- levels(metadata$gender)
colors_gender["Unknown"] <- adjustcolor("white", transparency_colors) # white for NA

# prepare colors for enrichment protocol
colors_enrichment_protocol <- paletteer_d("ggsci::default_igv", nlevels(metadata$enrichment_protocol)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_enrichment_protocol) <- levels(metadata$enrichment_protocol)
colors_enrichment_protocol["Unknown"] <- adjustcolor("white", transparency_colors) # white for NA

# prepare colors for dataset
colors_dataset <- paletteer_d("ggsci::default_nejm", nlevels(metadata$dataset)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_dataset) <- levels(metadata$dataset)

# prepare colors for treatment
colors_treatment <- paletteer_d("ggsci::default_jco", nlevels(metadata$treatment)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_treatment) <- levels(metadata$treatment)


# PCA ---------------------------------------------------------------------

# log-transform data
trdata <- as.matrix(xdata)
trdata <- log(trdata + 0.001) # log transform adding a small epsilon to avoid log(0)

# scale data to a mean of 0 and sd of 1 for each column (i.e, cell type)
trdata <- scale(trdata)
apply(trdata, 2, mean) |> round(3) # each column has mean ~ 0
apply(trdata, 2, sd) |> round(3) # each column has sd ~ 1
stopifnot(all(is.finite(trdata))) # no "bad" values

# run PCA
pca <- prcomp(trdata, center = FALSE, scale. = FALSE)

# prepare folder for PCA results
pca_folder <- file.path(output_folder, "PCA")
if (!dir.exists(pca_folder)) dir.create(pca_folder)

# export pca results
saveRDS(pca, file.path(pca_folder, "results_pca.rds"))

# plot PCs
pcs_to_plot <- paste0("PC", 1:3)
combinations <- combn(pcs_to_plot, 2) |>
  as.data.frame() |>
  as.list()
for (combination in combinations) {
  plot_pca_factor(
    pca_res = pca, # output of prcomp()
    PCs_to_plot = combination, # PCs to plot on X and Y axis respectively
    metadata = metadata, # dataframe of metadata (rows matching PCA data)
    vars_to_use = c( # which variables (columns of metadata) should be used?
      "response_6levels", "response_3levels", "response_2levels", "treatment", "gender", "enrichment_protocol", "dataset"
    ),
    colors_vars_list = list( # list of colors for the variables
      response_6levels = colors_response6,
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
    output_folder = pca_folder, # folder to save the plots
    file_format = "png", # file format
    res_ppi = 300, # resolution (pixels per inch)
    width_in = 6, height_in = 4 # width and height of the plot in inches
  )
}

# scree plot
percentVar <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
names(percentVar) <- colnames(pca$x)
png(file.path(pca_folder, "scree_plot.png"),
  width = 5 * resol, height = 4 * resol, res = resol
)
par(mar = c(3.5, 3, 0.5, 0.1), mgp = c(2, 0.8, 0), tcl = -0.3)
barplot(percentVar, las = 2, ylab = "% of total variance explained")
dev.off()


# PCA loadings ------------------------------------------------------------

set_dark_theme <- function() {
  par(
    bg = "black", # background of plot area
    fg = "white", # foreground: axis ticks, box
    col = "white", # default plotting color
    col.axis = "white", # axis tick labels
    col.lab = "white", # axis titles
    col.main = "white", # main title
    col.sub = "white" # subtitle (if used)
  )
}
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

plot_rotation <- function(PCx = "PC1", PCy = "PC2", loadings_data) {
  # set_dark_theme()
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

for (combination in combinations) {
  png(
    file.path(output_folder, paste0(
      "PCA_rotation_", combination[1], "_", combination[2], ".png"
    )),
    width = 5 * resol, height = 5 * resol, res = resol
  )
  plot_rotation(
    PCx = combination[1], PCy = combination[2], loadings_data = pca$rotation
  )
  dev.off()
}


# Correlation between all numeric variables -------------------------------

numeric_variables <- names(alldata)[sapply(alldata, is.numeric)]
cormat_all <- cor(alldata[, numeric_variables],
  method = "spearman",
  use = "pairwise.complete.obs"
)
png(file.path(output_folder, "corrplot_all_numeric.png"),
  width = 8 * resol, height = 6 * resol, res = resol
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

# Inspect collinearity between predictors ---------------------------------

# define threshold for cutting tree
hclust_height_thr <- 0.5

# correlation matrix of cell types
all(complete.cases(xdata)) # no missing values
cormat <- cor(xdata, method = "spearman") # calculate Spearman correlation

# cluster correlated cell types
hc <- hclust(as.dist(1 - abs(cormat)), method = "average")

# cut tree below a defined threshold
grp <- cutree(hc, h = hclust_height_thr) # all merges below threshold
table(grp)
unique(grp) |> length() # number of clusters

# get representative cell types for clustered groups
get_representatives <- function(grp, cormat) {
  split_names <- split(colnames(cormat), grp)
  reps <- lapply(split_names, function(vars) {
    sub <- abs(cormat[vars, vars, drop = FALSE])
    vars[which.max(rowMeans(sub))]
  })
  unlist(reps)
}
reps <- get_representatives(grp, cormat)
reps

# plot dendrogram
labels_marked <- ifelse( # mark representatives with a star
  cell_types %in% reps, paste0(cell_types, " *"), cell_types
)
ord <- hc$order
xx <- grp[ord] |> as.vector()
xx_groups <- names(table(xx))[table(xx) > 1]
xx <- unique(xx)
which_rects_indexes <- which(xx %in% xx_groups)
png(file.path(output_folder, "dendrogram_correlated_cibersortx.png"),
  width = 6 * resol, height = 5 * resol, res = resol
)
par(mar = c(3, 4.1, 2, 0), las = 1, xpd = TRUE, lwd = 1)
plot(
  hc,
  labels = labels_marked, xlab = "", cex = 0.7,
  main = "Dendrogram of immune cell types"
)
par(lwd = 2)
rect.hclust(hc, h = hclust_height_thr, which = which_rects_indexes, border = "#023E8A")
text(
  x = mean(par("usr")[1:2]),
  y = grconvertY(0, from = "ndc", to = "user"),
  labels = paste0("distance = 1 - |Spearman ρ|\nTree cut at height < ", hclust_height_thr),
  pos = 3
)
dev.off()

# plot correlation matrix
png(file.path(output_folder, "corrplot_cibersortx.png"),
  width = 8 * resol, height = 6 * resol, res = resol
)
corrplot(
  cormat[ord, ord],
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

# correlation matrix after selecting variables
cormat_selected <- cor(xdata[, as.vector(reps)], method = "spearman")
new_ord <- cell_types[ord]
new_ord <- new_ord[new_ord %in% reps]
png(file.path(output_folder, "corrplot_selected_cibersortx.png"),
  width = 8 * resol, height = 6 * resol, res = resol
)
corrplot(
  cormat_selected[new_ord, new_ord],
  order = "original",
  col = paletteer_c("grDevices::RdBu", 150, direction = -1),
  tl.col = "black", tl.cex = 0.8, tl.srt = 45,
  method = "circle",
  type = "full",
  is.corr = TRUE,
  title = "Spearman correlation matrix of selected variables",
  diag = TRUE,
  outline = FALSE,
  mar = c(0, 0, 1.5, 0)
)
dev.off()


# Export selected variables -----------------------------------------------

# summary table
cluster_members <- split(names(grp), grp)
summary_list <- lapply(names(cluster_members), function(g) {
  data.frame(
    cluster = as.integer(g),
    representative = reps[g],
    members = paste(cluster_members[[g]], collapse = ", "),
    row.names = NULL
  )
})
cluster_summary <- do.call(rbind, summary_list)
cluster_summary
write_xlsx(
  cluster_summary,
  file.path(output_folder, "clusters_correlated_cibersortx.xlsx")
)

# txt file with selected variables
writeLines(as.vector(reps), file.path(output_folder, "selected_cibersortx.txt"))


# Export correlation data -------------------------------------------------

cormat_triangle <- cormat
cormat_triangle[upper.tri(cormat_triangle, diag = TRUE)] <-
  NA # assign NA to upper triangle, including the diagonal
cor_df <- as.data.frame.table(cormat_triangle)
names(cor_df) <- c("Var1", "Var2", "spearman_correlation")
cor_df <- na.exclude(cor_df)
cor_df <- cor_df[order(abs(cor_df$spearman_correlation), decreasing = TRUE), ]
cor_df$spearman_correlation <- round(cor_df$spearman_correlation, 3)
write_xlsx(cor_df, file.path(output_folder, "correlations_cibersortx_data.xlsx"))


# Create scatterplots of strongly correlated variables --------------------

# create scatterplots of strongly correlated variables
xx_dir <- file.path(output_folder, "scatterplots_correlations")
if (!dir.exists(xx_dir)) dir.create(xx_dir)
strong_cor_df <- subset(cor_df, abs(spearman_correlation) > 0.5)
for (i in seq_len(nrow(strong_cor_df))) {
  x_cell <- strong_cor_df$Var1[i] |> as.character()
  y_cell <- strong_cor_df$Var2[i] |> as.character()
  png(file.path(xx_dir, paste0(
    "scatter.", i, ".", x_cell, "-", y_cell, ".png"
  )), height = 4.5 * resol, width = 4 * resol, res = resol)
  plot(
    x = rank(xdata[, x_cell]), y = rank(xdata[, y_cell]),
    las = 1, xlab = paste(x_cell, "(Rank)"), ylab = paste(y_cell, "(Rank)"), asp = 1,
    main = paste0("Spearman correlation = ", strong_cor_df$spearman_correlation[i]),
    pch = 21, bg = colors_dataset[as.numeric(metadata$dataset)]
  )
  dev.off()
}
