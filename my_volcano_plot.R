library(ggplot2)
library(ggrepel)

my_volcano_wLabels <- function(x, main = "", cutoff_log2FC = 1, cutoff_padj = 0.05,
                               deg_list = NULL, label_genes,
                               col_up_genes, col_down_genes, col_other_genes,
                               file_path = NULL, report_cutoffs = F, descr = "",
                               xlim_range = NULL) {
  # Data preparation
  x$neg_log10_padj <- -log10(x$padj)
  x$category <- "Other"
  x$color <- col_other_genes

  if (!is.null(deg_list)) {
    # options for relevant DEGs
    up_down_genes <- x[deg_list, ]
    up_down_genes$category[up_down_genes$log2FoldChange > 0] <-
      "Up-regulated"
    up_down_genes$category[up_down_genes$log2FoldChange < 0] <-
      "Down-regulated"
    up_down_genes$color[up_down_genes$category == "Up-regulated"] <-
      col_up_genes
    up_down_genes$color[up_down_genes$category == "Down-regulated"] <-
      col_down_genes
    x[deg_list, ] <- up_down_genes
  } else {
    x$category[x$padj < cutoff_padj & x$log2FoldChange > cutoff_log2FC] <- "Up-regulated"
    x$category[x$padj < cutoff_padj & x$log2FoldChange < -cutoff_log2FC] <- "Down-regulated"
    x$color[x$category == "Up-regulated"] <- col_up_genes
    x$color[x$category == "Down-regulated"] <- col_down_genes
  }

  # Create the base ggplot
  p <- ggplot(x, aes(x = log2FoldChange, y = neg_log10_padj, color = category, alpha = 0.3)) +
    geom_point(size = 2) +
    scale_color_manual(values = c(
      "Up-regulated" = col_up_genes,
      "Down-regulated" = col_down_genes,
      "Other" = col_other_genes
    )) +
    scale_fill_manual(values = c(
      "Up-regulated" = col_up_genes,
      "Down-regulated" = col_down_genes,
      "Other" = col_other_genes
    )) +
    scale_alpha_identity() + # Use the alpha column for transparency
    geom_hline(yintercept = -log10(cutoff_padj), linetype = "dashed") +
    geom_vline(xintercept = c(-cutoff_log2FC, cutoff_log2FC), linetype = "dashed") +
    labs(
      title = ifelse(report_cutoffs,
        yes = paste0(
          main, " (Cutoffs: |Log2FC| > ", cutoff_log2FC,
          " and Adjusted p-value < ", cutoff_padj, ")"
        ),
        no = ifelse((is.null(descr) | descr == ""), yes = main,
          no = paste0(main, " (", descr, ")")
        )
      ),
      x = expression(Log[2] ~ "Fold Change"),
      y = expression(-Log[10] ~ "Adjusted p-value"),
      color = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom", # Place the legend at the bottom
      legend.box = "vertical", # Stack legend items vertically
      plot.margin = margin(1, 1, 2, 1, unit = "lines") # Extra space for annotation
    )


  # Apply x-axis limits if specified
  if (!is.null(xlim_range)) {
    p <- p + xlim(xlim_range)
  }

  # Add labels with ggrepel
  labeled_data <- x[rownames(x) %in% label_genes, ]
  labeled_data$gene <- rownames(labeled_data)
  # black border for labelled points
  p <- p + geom_point(
    data = labeled_data,
    aes(
      x = log2FoldChange,
      y = neg_log10_padj,
      fill = category
    ),
    shape = 21,
    colour = "black",
    stroke = 0.6,
    size = 2.8,
    alpha = 1,
    show.legend = FALSE
  )
  # add labels
  p <- p + geom_label_repel(
    data = labeled_data,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    show.legend = FALSE,
    min.segment.length = 0,
    alpha = 1,
    segment.alpha = 1,
    label.size = 0.25, # spessore del bordo del box
    label.r = unit(0.15, "lines"), # angoli leggermente arrotondati
    fill = "white", # colore di sfondo del box
    colour = "black", # colore del testo
    arrow = arrow(length = unit(0.01, "npc"), type = "open")
  )


  # Save to file if file_path is provided
  if (!is.null(file_path)) {
    ggsave(file_path, plot = p, width = 8, height = 8)
  } else {
    print(p)
  }
}
