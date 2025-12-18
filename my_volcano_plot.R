library(ggplot2)
library(ggrepel)

my_volcano_plot <- function(x, main = "", cutoff_log2FC = 1,
                            cutoff_padj = 0.05, col_up_genes, col_down_genes,
                            col_other_genes, file_path = NULL) {
  x$neg_log10_padj <- -log10(x$padj)
  x$col <- col_other_genes
  x$padj_ok <- x$padj < cutoff_padj
  x$LFC_up <- x$log2FoldChange > cutoff_log2FC
  x$LFC_down <- -x$log2FoldChange > cutoff_log2FC
  x[x$padj_ok & x$LFC_up, "col"] <- col_up_genes
  x[x$padj_ok & x$LFC_down, "col"] <- col_down_genes
  # create volcano plot
  if (!is.null(file_path)) {
    pdf(file_path, width = 8, height = 8)
  }
  layout(mat = matrix(1:2, ncol = 1), heights = c(1, 3))
  # plot with title, legend and subtitle
  par(mar = c(3, 0.1, 4.1, 0.1))
  plot(NULL, axes = F, xlim = c(-1, 1), ylim = c(-1, 1), ylab = "", xlab = "", main = main, )
  mtext(
    side = 1, line = 0.5,
    text = paste0("cut-offs: |LFC| > ", cutoff_log2FC, " and p-adj < ", cutoff_padj)
  )
  legend(
    x = "center", pch = 20, col = c(col_down_genes, col_up_genes, col_other_genes),
    legend = c(
      paste0("down-regulated genes (", sum(x$col == col_down_genes), ")"),
      paste0("up-regulated genes (", sum(x$col == col_up_genes), ")"),
      paste0("other genes (", sum(x$col == col_other_genes), ")")
    ),
    horiz = T, bty = "n", y.intersp = 4
  )
  # actual volcano plot
  par(mar = c(4, 4.1, 0.1, 4.1))
  plot(NULL,
    main = "", las = 1,
    xlim = c(-max(abs(x$log2FoldChange)), max(abs(x$log2FoldChange))),
    ylim = c(0, max(x$neg_log10_padj + 3)),
    # xaxt = "n", yaxt="n",
    xlab = expression(Log[2] ~ "Fold Change"),
    ylab = expression(-Log[10] ~ "Adjusted p-value"),
  )
  xlabs <- sort(unique(c(
    pretty(x$log2FoldChange),
    cutoff_log2FC, -cutoff_log2FC
  )))
  yats <- sort(unique(c(seq(from = 0, to = 40, by = 10), -log10(cutoff_padj))))
  ylabs <- 10^(-yats)
  # axis(side=1, at=xlabs, labels = xlabs)
  # axis(side=2, at=yats, labels = ylabs, las=1)
  grid()
  abline(h = -log10(cutoff_padj), col = "black", lty = 2, lwd = 2)
  abline(v = c(cutoff_log2FC, -cutoff_log2FC), col = "black", lty = 2, lwd = 2)
  points(
    x = x$log2FoldChange, y = x$neg_log10_padj,
    pch = 20, # Circle points
    cex = 0.7, # Point size
    col = x$col
  )
  if (!is.null(file_path)) {
    dev.off()
  }
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
}


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

  # Set alpha levels for transparency
  x$alpha <- ifelse(x$category == "Other", 0.3, 0.7)

  # Create the base ggplot
  p <- ggplot(x, aes(x = log2FoldChange, y = neg_log10_padj, color = category, alpha = alpha)) +
    geom_point(size = 2) +
    scale_color_manual(values = c(
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
  p <- p + geom_text_repel(
    data = labeled_data,
    aes(label = gene),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
    max.overlaps = Inf
  )

  # Save to file if file_path is provided
  if (!is.null(file_path)) {
    ggsave(file_path, plot = p, width = 8, height = 8)
  } else {
    print(p)
  }
}
