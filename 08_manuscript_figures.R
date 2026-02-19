library(paletteer)
library(ggplot2)
library(ggrepel)
library(ggsankey)


# Plot settings -----------------------------------------------------------

# folder for plots
plots_folder <- "nonsync/plots_manuscript"
if (!dir.exists(plots_folder)) dir.create(plots_folder)

# resolution
resol <- 300

# Volcano plot NR vs. R ---------------------------------------------------
source("my_volcano_plot.R")

dea_folder <- "nonsync/autogo_response_pre_icb/NR_vs_R"
dir.exists(dea_folder)

# load DEA results
volcanodata <- read.table(
  file.path(dea_folder, "DE_NR_vs_R_allres.tsv"),
  header = TRUE,
  row.names = 1
)
# order by lowest adjusted p value
volcanodata <- volcanodata[order(volcanodata$padj), ]
# define thresholds
thr_lfc <- 1
thr_padj <- 0.05
# genes to label
sig_genes <- volcanodata |>
  subset(abs(log2FoldChange) > thr_lfc & padj < thr_padj) |>
  rownames()
n_top_genes_to_label <- 30
genes_to_label <- c( #
  head(sig_genes, n_top_genes_to_label), # top n genes by pvalue
  "C4B", "CRP", "MBL2", "KLRC1", "ITIH4" # other custom genes
) |> unique()
stopifnot(all(genes_to_label %in% sig_genes))
# create volcano plot
set.seed(123)
my_volcano_plot(
  dea_res = volcanodata, # data.frame with differential expression results; must contain columns 'log2FoldChange' and 'padj'; rownames = gene IDs
  main = paste0( # main title of the plot (character string)
    "Volcano plot of NR vs R\n",
    "with ABS(Log2FC) > ", thr_lfc,
    " and adjusted P-value < ", thr_padj
  ),
  cutoff_log2FC = thr_lfc, # numeric threshold for |log2FoldChange|
  cutoff_padj = thr_padj, # numeric adjusted p-value threshold
  genes_to_label = genes_to_label, # optional character vector of gene IDs to label on the plot
  col_up_genes = "red", # color for up-regulated genes (any valid R color)
  col_down_genes = "blue", # color for down-regulated genes (any valid R color)
  col_other_genes = "grey40", # color for non-significant genes (any valid R color)
  file_path = file.path(plots_folder, "volcano_NR_vs_R.png"), # optional file path with extension
  xlim_range = c( # optional numeric vector of length 2 to set x-axis limits with clipping at boundaries
    -max(abs(volcanodata$log2FoldChange)), max(abs(volcanodata$log2FoldChange))
  ),
  width_in = 6, # numeric width (in inches) used when saving the plot
  height_in = 6, # numeric height (in inches) used when saving the plot
  show_legend = TRUE # logical; if TRUE show legend, if FALSE hide legend
)


# Load data ---------------------------------------------------------------
metadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds") |> as.data.frame()
xdata <- readRDS("nonsync/01_clean_data/clean_cibersortx.rds") |> as.data.frame()
hed_data <- readRDS("nonsync/01_clean_data/clean_hed.rds") |> as.data.frame()
pca_res <- readRDS("nonsync/04_multicollinearity/PCA/results_pca.rds")


# Plot settings -----------------------------------------------------------

# folder for plots
plots_folder <- "nonsync/plots_manuscript"
if (!dir.exists(plots_folder)) dir.create(plots_folder)

# resolution
resol <- 300

# colors transparency
transparency_colors <- 0.8

# prepare colors for response (7 levels)
colors_response7 <- paletteer_c("viridis::viridis", nlevels(metadata$response_7levels)) |>
  adjustcolor(alpha.f = transparency_colors)
names(colors_response7) <- levels(metadata$response_7levels)

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

colors_vars_list <- list( # list of colors for the variables
  response_6levels = colors_response6,
  response_3levels = colors_response3,
  response_2levels = colors_response2,
  treatment = colors_treatment,
  gender = colors_gender,
  enrichment_protocol = colors_enrichment_protocol,
  dataset = colors_dataset
)


# PCA plots for PC1 vs. PC2 -----------------------------------------------

# function for drawing arrows corresponding to PCA loadings
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

# calculate % of variation in each PC
percent_var <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)
names(percent_var) <- colnames(pca_res$x)

# define symmetric axis range
xrange <- c(
  -max(abs(pca_res$x[, c("PC1", "PC2")])),
  max(abs(pca_res$x[, c("PC1", "PC2")]))
)

# prepare rotation data
rotdata <- pca_res$rotation[, c("PC1", "PC2")] |> as.data.frame()
rotdata$length <- sqrt(rotdata$PC1^2 + rotdata$PC2^2)
rotdata <- rotdata[order(rotdata$length, decreasing = FALSE), ]
rotdata <- tail(rotdata, 10) |> as.matrix()
colors_celltypes <- paletteer_d("ggsci::default_igv", nrow(rotdata))

# open graphic device
png(
  filename = file.path(plots_folder, "PCA_PC1PC2.png"),
  width = resol * 15, height = resol * 3.5, res = resol
)

# set graphical parameters
par(
  mfrow = c(1, 4),
  mar = c(4, 4, 4, 4), mgp = c(2.5, 0.6, 0),
  tcl = -0.4, las = 1, xpd = FALSE
)

# plot response
plot(NULL,
  xlim = xrange, ylim = xrange, asp = 1,
  bty = "l", ylab = paste0(
    "PC2", ": ", percent_var["PC2"],
    "% of total variation"
  ),
  xlab = paste0(
    "PC1", ": ", percent_var["PC1"],
    "% of total variation"
  ),
  main = "PCA by response"
)
abline(h = 0, col = "gray40", lty = "dotted")
abline(v = 0, col = "gray40", lty = "dotted")
points(
  x = pca_res$x[, "PC1"], y = pca_res$x[, "PC2"],
  pch = 21, cex = 1.5,
  bg = colors_vars_list[["response_2levels"]][as.numeric(
    metadata[, "response_2levels"]
  )]
)
legend(
  x = "bottomleft", legend = levels(metadata[, "response_2levels"]),
  pt.bg = colors_vars_list[["response_2levels"]], pch = 21, pt.cex = 1.5,
  border = NA, cex = 0.8, ncol = 2, xpd = TRUE
)

# plot enrichment protocol
plot(NULL,
  xlim = xrange, ylim = xrange, asp = 1,
  bty = "l", ylab = paste0(
    "PC2", ": ", percent_var["PC2"],
    "% of total variation"
  ),
  xlab = paste0(
    "PC1", ": ", percent_var["PC1"],
    "% of total variation"
  ),
  main = "PCA by enrichment protocol"
)
abline(h = 0, col = "gray40", lty = "dotted")
abline(v = 0, col = "gray40", lty = "dotted")
points(
  x = pca_res$x[, "PC1"], y = pca_res$x[, "PC2"],
  pch = 21, cex = 1.5,
  bg = colors_vars_list[["enrichment_protocol"]][as.numeric(
    metadata[, "enrichment_protocol"]
  )]
)
legend(
  x = "bottomright", legend = names(colors_vars_list[["enrichment_protocol"]]),
  pt.bg = colors_vars_list[["enrichment_protocol"]], pch = 21, pt.cex = 1.5,
  border = NA, cex = 0.8, ncol = 1, xpd = TRUE
)

# plot dataset
plot(NULL,
  xlim = xrange, ylim = xrange, asp = 1,
  bty = "l", ylab = paste0(
    "PC2", ": ", percent_var["PC2"],
    "% of total variation"
  ),
  xlab = paste0(
    "PC1", ": ", percent_var["PC1"],
    "% of total variation"
  ),
  main = "PCA by dataset"
)
abline(h = 0, col = "gray40", lty = "dotted")
abline(v = 0, col = "gray40", lty = "dotted")
points(
  x = pca_res$x[, "PC1"], y = pca_res$x[, "PC2"],
  pch = 21, cex = 1.5,
  bg = colors_vars_list[["dataset"]][as.numeric(
    metadata[, "dataset"]
  )]
)
legend(
  x = "bottomleft", legend = levels(metadata[, "dataset"]),
  pt.bg = colors_vars_list[["dataset"]], pch = 21, pt.cex = 1.5,
  border = NA, cex = 0.8, ncol = 1, xpd = TRUE
)

# plot rotation
par(las = 1)
plot(NULL,
  xlim = c(-0.3, 0.7), ylim = c(-0.5, 0.5),
  xlab = "Loading on PC1", ylab = "Loading on PC2",
  xaxs = "i", yaxs = "i", asp = 1, bty = "l", # axes = FALSE,
  main = paste0("PCA Loading"),
)
abline(h = 0, col = "gray40", lty = "dotted")
abline(v = 0, col = "gray40", lty = "dotted")
par(xpd = TRUE)
arrows_pca_loadings(
  loadings_data = rotdata,
  text.cex = 1, col.text = colors_celltypes, lwd = 2,
  length = 0.05, col = colors_celltypes
)


# close graphic device
dev.off()

# PCA plots for PC1 vs. PC3 -----------------------------------------------

# define symmetric axis range
xrange <- c(
  -max(abs(pca_res$x[, c("PC1", "PC3")])),
  max(abs(pca_res$x[, c("PC1", "PC3")]))
)

# prepare rotation data
rotdata <- pca_res$rotation[, c("PC1", "PC3")] |> as.data.frame()
rotdata$length <- sqrt(rotdata$PC1^2 + rotdata$PC3^2)
rotdata <- rotdata[order(rotdata$length, decreasing = FALSE), ]
rotdata <- tail(rotdata, 10) |> as.matrix()
colors_celltypes <- paletteer_d("ggsci::default_igv", nrow(rotdata))

# open graphic device
png(
  filename = file.path(plots_folder, "PCA_PC1PC3.png"),
  width = resol * 8, height = resol * 8, res = resol
)

# set graphical parameters
par(
  mfrow = c(2, 2),
  mar = c(2.5, 2.5, 2.5, 0.2), mgp = c(1.5, 0.4, 0),
  tcl = -0.2, las = 1, xpd = FALSE
)

# plot response
plot(NULL,
  xlim = xrange, ylim = xrange, asp = 1,
  bty = "l", ylab = paste0(
    "PC3", ": ", percent_var["PC3"],
    "% of total variation"
  ),
  xlab = paste0(
    "PC1", ": ", percent_var["PC1"],
    "% of total variation"
  ),
  main = "PCA by response"
)
abline(h = 0, col = "gray40", lty = "dotted")
abline(v = 0, col = "gray40", lty = "dotted")
points(
  x = pca_res$x[, "PC1"], y = pca_res$x[, "PC3"],
  pch = 21, cex = 1.5,
  bg = colors_vars_list[["response_2levels"]][as.numeric(
    metadata[, "response_2levels"]
  )]
)
legend(
  x = "bottomleft", legend = levels(metadata[, "response_2levels"]),
  pt.bg = colors_vars_list[["response_2levels"]], pch = 21, pt.cex = 1.5,
  border = NA, cex = 0.8, ncol = 2, xpd = TRUE
)

# plot enrichment protocol
plot(NULL,
  xlim = xrange, ylim = xrange, asp = 1,
  bty = "l", ylab = paste0(
    "PC3", ": ", percent_var["PC3"],
    "% of total variation"
  ),
  xlab = paste0(
    "PC1", ": ", percent_var["PC1"],
    "% of total variation"
  ),
  main = "PCA by enrichment protocol"
)
abline(h = 0, col = "gray40", lty = "dotted")
abline(v = 0, col = "gray40", lty = "dotted")
points(
  x = pca_res$x[, "PC1"], y = pca_res$x[, "PC3"],
  pch = 21, cex = 1.5,
  bg = colors_vars_list[["enrichment_protocol"]][as.numeric(
    metadata[, "enrichment_protocol"]
  )]
)
legend(
  x = "bottomleft", legend = names(colors_vars_list[["enrichment_protocol"]]),
  pt.bg = colors_vars_list[["enrichment_protocol"]], pch = 21, pt.cex = 1.5,
  border = NA, cex = 0.8, ncol = 1, xpd = TRUE
)

# plot dataset
plot(NULL,
  xlim = xrange, ylim = xrange, asp = 1,
  bty = "l", ylab = paste0(
    "PC3", ": ", percent_var["PC3"],
    "% of total variation"
  ),
  xlab = paste0(
    "PC1", ": ", percent_var["PC1"],
    "% of total variation"
  ),
  main = "PCA by dataset"
)
abline(h = 0, col = "gray40", lty = "dotted")
abline(v = 0, col = "gray40", lty = "dotted")
points(
  x = pca_res$x[, "PC1"], y = pca_res$x[, "PC3"],
  pch = 21, cex = 1.5,
  bg = colors_vars_list[["dataset"]][as.numeric(
    metadata[, "dataset"]
  )]
)
legend(
  x = "bottomleft", legend = levels(metadata[, "dataset"]),
  pt.bg = colors_vars_list[["dataset"]], pch = 21, pt.cex = 1.5,
  border = NA, cex = 0.8, ncol = 1, xpd = TRUE
)

# plot rotation
par(las = 1)
plot(NULL,
  xlim = c(-0.4, 0.8), ylim = c(-0.5, 0.5),
  xlab = "Loading on PC1", ylab = "Loading on PC3",
  xaxs = "i", yaxs = "i", asp = 1, bty = "l",
  main = paste0("PCA Loading"),
)
abline(h = 0, col = "gray40", lty = "dotted")
abline(v = 0, col = "gray40", lty = "dotted")
par(xpd = TRUE)
arrows_pca_loadings(
  PCx = "PC1", PCy = "PC3",
  loadings_data = rotdata,
  text.cex = 1, col.text = colors_celltypes, lwd = 2,
  length = 0.05, col = colors_celltypes
)


# close graphic device
dev.off()


# Sankey plot of response classes collapsing ------------------------------

# prepare data and colors
xx <- metadata[, c("dataset", grepv("response", names(metadata)))]
names(xx)
long_data <- make_long(
  xx, dataset, response_7levels, response_6levels, response_3levels, response_2levels
)
levels(long_data$x)
levels(long_data$x) <- levels(long_data$next_x) <- c(
  "Dataset", "Response (7)", "Response (6)", "Response (3)", "Response (2)"
)
long_data$node_chr <- long_data$node
long_data$node <- factor(long_data$node, levels = c(
  levels(xx$dataset), levels(xx$response_7levels)
))
pal_node <- c(colors_response7, colors_dataset)
long_data$fill_node <- unname(pal_node[long_data$node_chr])
leg_dataset <- data.frame(key = names(colors_dataset))
leg_response <- data.frame(key = names(colors_response7))

# create plot
pl <- ggplot(long_data, aes(
  x = x, next_x = next_x,
  node = node, next_node = next_node,
  label = node
)) +
  geom_sankey(
    aes(fill = I(fill_node)),
    flow.alpha = 0.6,
    show.legend = FALSE,
    position = "identity",
    type = "sankey"
  ) +
  geom_sankey_label(size = 3, color = "black", fill = "white") +
  # invisible legend layers
  geom_point(
    data = leg_dataset,
    aes(x = 1, y = 1, fill = key),
    inherit.aes = FALSE,
    alpha = 0,
    shape = 22
  ) +
  geom_point(
    data = leg_response,
    aes(x = 1, y = 1, color = key),
    inherit.aes = FALSE,
    alpha = 0
  ) +
  scale_fill_manual(
    name   = "Dataset",
    values = colors_dataset,
    breaks = names(colors_dataset)
  ) +
  scale_color_manual(
    name   = "Response Class",
    values = colors_response7,
    breaks = names(colors_response7)
  ) +
  guides(
    fill  = guide_legend(ncol = 1, override.aes = list(alpha = 1, shape = 22, size = 6, colour = NA)),
    color = guide_legend(ncol = 1, override.aes = list(alpha = 1, shape = 15, size = 6, stroke = 0))
  ) +
  # theme
  theme_void() +
  theme(
    axis.text.x = element_text(color = "black", size = 10),
    legend.position = "right",
    legend.box = "vertical",
    plot.background = element_rect(fill = "white")
  )

# save plot
ggsave(
  plot = pl, filename = file.path(plots_folder, "sankey_response_classes.png"),
  width = 6, height = 4, units = "in", dpi = 300
)
