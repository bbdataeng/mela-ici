library(ggplot2)
library(ggrepel)

my_volcano_plot <- function(
  dea_res, # data.frame with differential expression results; must contain columns 'log2FoldChange' and 'padj'; rownames = gene IDs
  main = "", # main title of the plot (character string)
  cutoff_log2FC = 1, # numeric threshold for |log2FoldChange|
  cutoff_padj = 0.05, # numeric adjusted p-value threshold
  deg_list = NULL, # optional character vector specifying genes to classify as DE instead of using cutoffs
  genes_to_label = NULL, # optional character vector of gene IDs to label on the plot
  col_up_genes, # color for up-regulated genes (any valid R color)
  col_down_genes, # color for down-regulated genes (any valid R color)
  col_other_genes, # color for non-significant genes (any valid R color)
  file_path = NULL, # optional file path with extension
  report_cutoffs = FALSE, # if TRUE, include cutoff values in the plot title
  descr = "", # optional short description appended to the title (ignored if report_cutoffs = TRUE)
  xlim_range = NULL, # optional numeric vector of length 2 to set x-axis limits with clipping at boundaries
  width_in = 8, # numeric width (in inches) used when saving the plot
  height_in = 8, # numeric height (in inches) used when saving the plot
  show_legend = TRUE # logical; if TRUE show legend, if FALSE hide legend
) {
  # basic input checks
  req_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(req_cols, colnames(dea_res))
  if (length(missing_cols) > 0) {
    stop("dea_res is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  if (is.null(rownames(dea_res))) {
    stop("dea_res must have rownames (gene IDs).")
  }

  ### data preparation ###
  # handle padj edge cases: NA and 0
  padj <- dea_res$padj
  padj_safe <- padj
  padj_safe[is.na(padj_safe)] <- NA_real_
  padj_safe[!is.na(padj_safe) & padj_safe <= 0] <- .Machine$double.xmin

  dea_res$neg_log10_padj <- -log10(padj_safe)
  dea_res$category <- "Non-significant"

  ### X clipping (optional) ###
  # if xlim_range is provided, clamp log2FC to the boundaries and mark clamped points
  dea_res$x_plot <- dea_res$log2FoldChange
  dea_res$clipped <- FALSE
  if (!is.null(xlim_range)) {
    if (length(xlim_range) != 2 || any(!is.finite(xlim_range))) {
      stop("xlim_range must be a numeric vector of length 2 with finite values (e.g., c(-6, 6)).")
    }
    x_low <- min(xlim_range)
    x_high <- max(xlim_range)

    dea_res$clipped <- !is.na(dea_res$log2FoldChange) & (dea_res$log2FoldChange < x_low | dea_res$log2FoldChange > x_high)
    dea_res$x_plot <- pmin(pmax(dea_res$log2FoldChange, x_low), x_high)
  }

  ### categorization ###
  if (!is.null(deg_list)) {
    # allow both character vector of gene IDs or logical/integer indexing
    if (is.character(deg_list)) {
      present <- intersect(deg_list, rownames(dea_res))
      missing <- setdiff(deg_list, rownames(dea_res))
      if (length(missing) > 0) {
        warning("deg_list contains ", length(missing), " gene(s) not found in rownames(dea_res). They will be ignored.")
      }
      idx <- present
    } else {
      # integer/logical indexing
      idx <- deg_list
    }

    up_down_genes <- dea_res[idx, , drop = FALSE]
    up_down_genes$category[up_down_genes$log2FoldChange > 0] <- "Up-regulated"
    up_down_genes$category[up_down_genes$log2FoldChange < 0] <- "Down-regulated"
    dea_res[idx, "category"] <- up_down_genes$category
  } else {
    dea_res$category[dea_res$padj < cutoff_padj & dea_res$log2FoldChange > cutoff_log2FC] <- "Up-regulated"
    dea_res$category[dea_res$padj < cutoff_padj & dea_res$log2FoldChange < -cutoff_log2FC] <- "Down-regulated"
  }

  # force legend/category order (even if a category is absent)
  legend_order <- c("Down-regulated", "Non-significant", "Up-regulated")
  dea_res$category <- factor(dea_res$category, levels = legend_order)

  ### summary of categories (used in legend labels) ###
  category_counts <- table(dea_res$category)
  n_up <- if ("Up-regulated" %in% names(category_counts)) category_counts["Up-regulated"] else 0
  n_down <- if ("Down-regulated" %in% names(category_counts)) category_counts["Down-regulated"] else 0
  n_other <- if ("Non-significant" %in% names(category_counts)) category_counts["Non-significant"] else 0

  legend_labels <- c(
    "Down-regulated" = paste0("Down-regulated (N = ", n_down, ")"),
    "Non-significant" = paste0("Non-significant (N = ", n_other, ")"),
    "Up-regulated" = paste0("Up-regulated (N = ", n_up, ")")
  )

  legend_breaks <- legend_order

  ### base plot ###
  p <- ggplot(
    dea_res,
    aes(
      x = x_plot,
      y = neg_log10_padj,
      color = category,
      shape = clipped
    )
  ) +
    geom_point(size = 2, alpha = 0.3, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Up-regulated" = col_up_genes,
        "Down-regulated" = col_down_genes,
        "Non-significant" = col_other_genes
      ),
      breaks = legend_breaks,
      labels = unname(legend_labels[legend_breaks]),
      drop = FALSE
    ) +
    scale_fill_manual(
      values = c(
        "Up-regulated" = col_up_genes,
        "Down-regulated" = col_down_genes,
        "Non-significant" = col_other_genes
      ),
      breaks = legend_breaks,
      labels = unname(legend_labels[legend_breaks]),
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 17)
    ) +
    guides(
      shape = "none",
      fill = "none",
      color = guide_legend(
        override.aes = list(
          alpha = 1,
          shape = 16,
          size = 3,
          stroke = 0
        )
      )
    ) +
    geom_hline(yintercept = -log10(cutoff_padj), linetype = "dashed") +
    geom_vline(xintercept = c(-cutoff_log2FC, cutoff_log2FC), linetype = "dashed") +
    labs(
      title = ifelse(
        report_cutoffs,
        yes = paste0(
          main, " (Cutoffs: |Log2FC| > ", cutoff_log2FC,
          " and Adjusted p-value < ", cutoff_padj, ")"
        ),
        no = ifelse(
          (is.null(descr) | descr == ""),
          yes = main,
          no = paste0(main, " (", descr, ")")
        )
      ),
      x = expression(Log[2] ~ "Fold Change"),
      y = expression(-Log[10] ~ "Adjusted p-value"),
      color = NULL,
      shape = NULL
    ) +
    theme_minimal() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.box = "vertical",
      plot.margin = margin(1, 1, 1, 1, unit = "lines")
    )

  # apply x-axis limits
  if (!is.null(xlim_range)) {
    p <- p + scale_x_continuous(limits = xlim_range)
  }

  ### gene labels (optional) ###
  if (!is.null(genes_to_label)) {
    genes_to_label <- unique(genes_to_label)
    present_lab <- intersect(genes_to_label, rownames(dea_res))
    missing_lab <- setdiff(genes_to_label, rownames(dea_res))
    if (length(missing_lab) > 0) {
      warning("genes_to_label contains ", length(missing_lab), " gene(s) not found in rownames(dea_res). They will be ignored.")
    }

    labeled_data <- dea_res[present_lab, , drop = FALSE]
    if (nrow(labeled_data) > 0) {
      labeled_data$gene <- rownames(labeled_data)

      labeled_in <- labeled_data[!labeled_data$clipped, , drop = FALSE]
      labeled_out <- labeled_data[labeled_data$clipped, , drop = FALSE]

      # black border for labelled points (within x-limits)
      if (nrow(labeled_in) > 0) {
        p <- p + geom_point(
          data = labeled_in,
          aes(
            x = x_plot,
            y = neg_log10_padj,
            fill = category
          ),
          shape = 21,
          colour = "black",
          stroke = 0.6,
          size = 2.8,
          alpha = 1,
          show.legend = FALSE,
          na.rm = TRUE
        )
      }

      # black border for labelled points (clipped to x-limits) - triangle
      if (nrow(labeled_out) > 0) {
        p <- p + geom_point(
          data = labeled_out,
          aes(
            x = x_plot,
            y = neg_log10_padj,
            fill = category
          ),
          shape = 24,
          colour = "black",
          stroke = 0.6,
          size = 3.1,
          alpha = 1,
          show.legend = FALSE,
          na.rm = TRUE
        )
      }

      # add labels (for all labeled genes)
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
        label.size = 0.25,
        label.r = unit(0.15, "lines"),
        fill = "white",
        colour = "black",
        arrow = arrow(length = unit(0.01, "npc"), type = "open")
      )
    }
  }

  # save to file if file_path is provided
  if (!is.null(file_path)) {
    ggsave(file_path, plot = p, width = width_in, height = height_in)
  } else {
    print(p)
  }
}
