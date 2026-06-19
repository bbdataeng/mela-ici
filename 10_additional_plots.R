# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(paletteer) # v1.7.0
  library(ggplot2) # v4.0.3
  library(ggrepel) # v0.9.8
  library(ggsankey) # v0.0.99999
})

# Prepare output folders --------------------------------------------------

# main output folders
outdirs <- c(
  # output folder for complete dataset
  complete = "nonsync/10_additional_plots/complete",
  # output folder for dataset without checkmate067
  nocheckmate067 = "nonsync/10_additional_plots/nocheckmate067"
)

# create directories
for (i in seq_along(outdirs)) {
  if (!dir.exists(outdirs[[i]])) dir.create(outdirs[[i]], recursive = TRUE)
}


# Load data ---------------------------------------------------------------

metadata_complete <- readRDS("nonsync/04_clean_data/clean_metadata.rds")
xdata_complete <- readRDS("nonsync/04_clean_data/clean_cibersortx.rds")
hed_complete <- readRDS("nonsync/04_clean_data/clean_hed.rds")
pca_res <- list(
  complete = readRDS("nonsync/06_PCA_multicollinearity/complete/PCA/results_pca.rds"),
  nocheckmate067 = readRDS("nonsync/06_PCA_multicollinearity/nocheckmate067/PCA/results_pca.rds")
)

# clean anatomical site
metadata_complete$anatomical_location_s1 <- as.factor(metadata_complete$anatomical_location_s1) |> addNA()
levels(metadata_complete$anatomical_location_s1)[nlevels(metadata_complete$anatomical_location_s1)] <- "Unknown"

# identify checkmate067 cohort
to_exclude <- which(
  metadata_complete$dataset == "Campbell-2023" &
    metadata_complete$enrichment_protocol == "targeted-mRNA-capture"
)

# create lists of data
alldata <- list(
  complete = cbind(metadata_complete, hed_complete, xdata_complete),
  nocheckmate067 = cbind(metadata_complete, hed_complete, xdata_complete)[-to_exclude, ]
)
rm(metadata_complete, hed_complete, xdata_complete, to_exclude)


# Plot settings and colors ------------------------------------------------

# resolution in pixels per inch
resol <- 500

# prepare colors for response (4 levels)
colors_response4 <- paletteer_c("viridis::viridis", nlevels(alldata[[1]]$response_4levels))
names(colors_response4) <- levels(alldata[[1]]$response_4levels)

# prepare colors for response (2 levels)
colors_response2 <- c("NR" = "#D5B3FF", "R" = "#FBD960")

# prepare colors for enrichment protocol
alldata <- lapply(alldata, function(x) {
  x$enrichment_protocol <- addNA(x$enrichment_protocol)
  levels(x$enrichment_protocol)[nlevels(x$enrichment_protocol)] <- "Unknown"
  x
})
colors_enrichment_protocol <- c("#6A8532", "#6A4C93", "#D95D39", "#C9ADA7")
names(colors_enrichment_protocol) <- levels(alldata[[1]]$enrichment_protocol)

# prepare colors for dataset
colors_dataset <- c(
  "Hugo-2016" = "#E0BE36",
  "Riaz-2017" = "#00A0D1",
  "Auslander-2018" = "#C8DE7B",
  "Gide-2019" = "#FB6F92",
  "Du-2021" = "#E07A5F",
  "Campbell-2023" = "#8B0000"
)
stopifnot(all(names(colors_dataset) %in% alldata[[1]]$dataset))
alldata <- lapply(alldata, function(x) {
  x$dataset <- factor(x$dataset, levels = names(colors_dataset))
  x
})

# prepare colors for treatment
colors_treatment <- paletteer_d("ggsci::default_jco", nlevels(alldata[[1]]$treatment))
names(colors_treatment) <- levels(alldata[[1]]$treatment)

# prepare colors for anatomical site
colors_anatomy <- paletteer_d("ggsci::default_nejm", nlevels(alldata[[1]]$anatomical_location_s1))
names(colors_anatomy) <- levels(alldata[[1]]$anatomical_location_s1)
colors_anatomy["Unknown"] <- "white"


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
  file_path = file.path(outdirs[x], "volcano_NR_vs_R.png"), # optional file path with extension
  xlim_range = c( # optional numeric vector of length 2 to set x-axis limits with clipping at boundaries
    -max(abs(volcanodata$log2FoldChange)), max(abs(volcanodata$log2FoldChange))
  ),
  width_in = 6, # numeric width (in inches) used when saving the plot
  height_in = 6, # numeric height (in inches) used when saving the plot
  show_legend = TRUE # logical; if TRUE show legend, if FALSE hide legend
)


# Prepare PCA plots -------------------------------------------------------

# function for PCA plot
make_pca_plot <- function(
  pca_res, PCx = "PC1", PCy = "PC2", by_var, xlim, ylim,
  main, legend_pos = "bottomleft", legend_ncol = 1, cols
) {
  # calculate % of variation in each PC
  percent_var <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)
  names(percent_var) <- colnames(pca_res$x)
  # empty plot
  plot(NULL,
    xlim = xlim, ylim = ylim, asp = 1, bty = "l",
    ylab = paste0(
      PCy, ": ", percent_var[PCy], "% of total variation"
    ),
    xlab = paste0(
      PCx, ": ", percent_var[PCx], "% of total variation"
    ),
    main = main
  )
  # grid lines at 0
  abline(h = 0, col = "gray40", lty = "dotted")
  abline(v = 0, col = "gray40", lty = "dotted")
  # add points
  points(
    x = pca_res$x[, PCx], y = pca_res$x[, PCy],
    pch = 21, cex = 1.5, bg = cols[as.numeric(by_var)]
  )
  # add legend
  legend(
    x = legend_pos, legend = levels(by_var),
    pt.bg = cols, pch = 21, pt.cex = 1.5,
    border = NA, cex = 0.8, ncol = legend_ncol, xpd = TRUE
  )
}

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


# Plots PC2~PC1 -----------------------------------------------------------

for (x in names(outdirs)) {
  # prepare rotation data and colors
  rotdata <- pca_res[[x]]$rotation[, c("PC1", "PC2")] |> as.data.frame()
  rotdata$length <- sqrt(rotdata$PC1^2 + rotdata$PC2^2)
  rotdata <- rotdata[order(rotdata$length, decreasing = FALSE), ]
  rotdata <- tail(rotdata, 10) |> as.matrix()
  colors_celltypes <- paletteer_d("ggsci::default_igv", nrow(rotdata))

  # define symmetric axis range
  xrange <- c(
    -max(abs(pca_res[[x]]$x[, c("PC1", "PC2")])),
    max(abs(pca_res[[x]]$x[, c("PC1", "PC2")]))
  ) * 1.2

  # open graphic device
  png(
    filename = file.path(outdirs[x], "PCA_PC1PC2.png"),
    width = resol * 18, height = resol * 3.5, res = resol
  )

  # set graphical parameters
  par(
    mfrow = c(1, 5),
    mar = c(4, 4, 4, 4), mgp = c(2.5, 0.6, 0),
    tcl = -0.4, las = 1, xpd = FALSE, font.main = 1
  )

  # plot response
  make_pca_plot(
    pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC2",
    by_var = alldata[[x]]$response_2levels,
    xlim = xrange, ylim = xrange + 0.7, main = "Response",
    legend_ncol = 2, legend_pos = "topright", cols = colors_response2
  )

  # plot treatment
  make_pca_plot(
    pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC2",
    by_var = alldata[[x]]$treatment,
    xlim = xrange, ylim = xrange + 0.7, main = "Treatment",
    legend_ncol = 2, legend_pos = "topright", cols = colors_treatment
  )

  # plot enrichment protocol
  make_pca_plot(
    pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC2",
    by_var = alldata[[x]]$enrichment_protocol,
    xlim = xrange, ylim = xrange + 0.7,
    main = "Enrichment Protocol", legend_pos = "topright",
    legend_ncol = 2, cols = colors_enrichment_protocol
  )

  # plot dataset
  make_pca_plot(
    pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC2",
    by_var = alldata[[x]]$dataset,
    xlim = xrange, ylim = xrange + 0.7, legend_ncol = 2,
    main = "Dataset", legend_pos = "topright", cols = colors_dataset
  )

  # plot rotation
  par(las = 1)
  plot(NULL,
    xlim = c(-0.3, 0.7), ylim = c(-0.5, 0.5),
    xlab = "Loading on PC1", ylab = "Loading on PC2",
    xaxs = "i", yaxs = "i", asp = 1, bty = "l",
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
}


# Plots for anatomical location -------------------------------------------

for (x in names(outdirs)) {
  # open graphic device
  png(
    filename = file.path(outdirs[x], "PCA_PC1PC2PC3_anatomical_location.png"),
    width = resol * 3.5, height = resol * 18, res = resol
  )


  # set graphical parameters
  par(
    mfrow = c(5, 1),
    mar = c(4, 4, 4, 4), mgp = c(2.5, 0.6, 0),
    tcl = -0.4, las = 1, xpd = FALSE, font.main = 1
  )

  # plot anatomical location on PC2~PC1
  xrange <- c(
    -max(abs(pca_res[[x]]$x[, c("PC1", "PC2")])),
    max(abs(pca_res[[x]]$x[, c("PC1", "PC2")]))
  ) * 1.2
  make_pca_plot(
    pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC2",
    by_var = alldata[[x]]$anatomical_location_s1,
    xlim = xrange, ylim = c(-5.3, 6.7), main = "Anatomical location",
    legend_ncol = 2, legend_pos = "topright", cols = colors_anatomy
  )

  # plot anatomical location on PC3~PC1
  xrange <- c(
    -max(abs(pca_res[[x]]$x[, c("PC1", "PC3")])),
    max(abs(pca_res[[x]]$x[, c("PC1", "PC3")]))
  ) * 1.2
  make_pca_plot(
    pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC3",
    by_var = alldata[[x]]$anatomical_location_s1,
    xlim = xrange, ylim = c(-5.3, 6.7), main = "Anatomical location",
    legend_ncol = 2, legend_pos = "bottomleft", cols = colors_anatomy
  )

  # close graphic device
  dev.off()
}


# Plots PC3~PC1 -----------------------------------------------------------

for (x in names(outdirs)) {
  # define symmetric axis range
  xrange <- c(
    -max(abs(pca_res[[x]]$x[, c("PC1", "PC3")])),
    max(abs(pca_res[[x]]$x[, c("PC1", "PC3")]))
  ) * 1.2

  # prepare rotation data and colors
  rotdata <- pca_res[[x]]$rotation[, c("PC1", "PC3")] |> as.data.frame()
  rotdata$length <- sqrt(rotdata$PC1^2 + rotdata$PC3^2)
  rotdata <- rotdata[order(rotdata$length, decreasing = FALSE), ]
  rotdata <- tail(rotdata, 10) |> as.matrix()
  colors_celltypes <- paletteer_d("ggsci::default_igv", nrow(rotdata))

  # open graphic device
  png(
    filename = file.path(outdirs[x], "PCA_PC1PC3.png"),
    width = resol * 14.4, height = resol * 3.5, res = resol
  )

  # set graphical parameters
  par(
    mfrow = c(1, 4),
    mar = c(4, 4, 4, 4), mgp = c(2.5, 0.6, 0),
    tcl = -0.4, las = 1, xpd = FALSE, font.main = 1
  )

  # plot response
  make_pca_plot(
    pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC3",
    by_var = alldata[[x]]$response_2levels,
    xlim = xrange, ylim = c(-5.3, 6.7), main = "Response",
    legend_ncol = 2, legend_pos = "bottomleft", cols = colors_response2
  )

  # plot treatment
  make_pca_plot(
    pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC3",
    by_var = alldata[[x]]$treatment,
    xlim = xrange, ylim = c(-5.3, 6.7), main = "Treatment",
    legend_ncol = 2, legend_pos = "bottomleft", cols = colors_treatment
  )

  # # plot enrichment protocol
  # make_pca_plot(
  #   pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC3",
  #   by_var = alldata[[x]]$enrichment_protocol,
  #   xlim = xrange, ylim = c(-5.3, 6.7),
  #   main = "Enrichment Protocol", legend_pos = "bottomleft",
  #   legend_ncol = 2, cols = colors_enrichment_protocol
  # )

  # plot dataset
  make_pca_plot(
    pca_res = pca_res[[x]], PCx = "PC1", PCy = "PC3",
    by_var = alldata[[x]]$dataset,
    xlim = xrange, ylim = c(-5.3, 6.7), legend_ncol = 2,
    main = "Dataset", legend_pos = "bottomleft", cols = colors_dataset
  )

  # plot rotation
  par(las = 1)
  plot(NULL,
    xlim = c(-0.4, 1), ylim = c(-0.5, 0.5),
    xlab = "Loading on PC1", ylab = "Loading on PC3",
    asp = 1, bty = "l", main = paste0("Loadings"),
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
}


# Sankey plot of response classes collapsing ------------------------------

for (x in names(outdirs)) {
  # prepare data and colors
  xx <- alldata[[x]][, c("dataset", grepv("response", names(alldata[[x]])))]
  names(xx)
  long_data <- make_long(
    xx, dataset, response_4levels, response_3levels, response_2levels
  )
  levels(long_data$x)
  levels(long_data$x) <- levels(long_data$next_x) <- c(
    "Dataset", "RECIST", "3-classes", "Binary"
  )
  long_data$node_chr <- long_data$node
  long_data$node <- factor(long_data$node, levels = c(
    levels(xx$dataset), "PD", "NR", "SD", "PR", "R", "CR"
  ))
  pal_node <- c(colors_response4, colors_response2, colors_dataset)
  long_data$fill_node <- unname(pal_node[long_data$node_chr])
  leg_dataset <- data.frame(key = names(colors_dataset))
  leg_response <- data.frame(key = c(names(colors_response4), names(colors_response2)))

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
      name = "Dataset",
      values = colors_dataset,
      breaks = names(colors_dataset)
    ) +
    scale_color_manual(
      name = "Response Class",
      values = c(colors_response4, colors_response2),
      breaks = c(names(colors_response4), names(colors_response2))
    ) +
    guides(
      fill = guide_legend(ncol = 1, override.aes = list(alpha = 1, shape = 22, size = 6, colour = NA)),
      color = guide_legend(ncol = 1, override.aes = list(alpha = 1, shape = 15, size = 6, stroke = 0))
    ) +
    # theme
    theme_void() +
    theme(
      axis.text.x = element_text(color = "black", size = 10, angle = 0),
      legend.position = "right",
      legend.box = "vertical",
      plot.background = element_rect(fill = "white")
    )

  # save plot
  ggsave(
    plot = pl, filename = file.path(outdirs[x], "sankey_response_classes.png"),
    width = 5, height = 4, units = "in", dpi = 300
  )
}
