# Load libraries ----------------------------------------------------------
library(paletteer)
library(corrplot)
library(writexl)

# Prepare folder for figures ----------------------------------------------
output_folder <- "nonsync/04_multicollinearity"
if (!dir.exists(output_folder)) dir.create(output_folder)

# Load clean data ---------------------------------------------------------
metadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds")
xdata <- readRDS("nonsync/01_clean_data/clean_cibersortx.rds")

cell_types <- colnames(xdata)

# sanitize cell types names
cell_types_sanitized <- gsub(" ", "_", cell_types)
cell_types_sanitized <- gsub("\\(|\\)", "", cell_types_sanitized)
# colnames(xdata) <- cell_types_sanitized

# put all data together
alldata <- cbind(metadata, xdata)
rm(metadata, xdata) # just to free some memory

# exclude biopsies on treatment
testdata <- subset(alldata, biopsy_time == "PRE-ICB") |> droplevels()
nrow(testdata) # now 170 rows

# Plot settings -----------------------------------------------------------
resol <- 300

# Inspect collinearity between predictors ---------------------------------

# define threshold for cutting tree
hclust_height_thr <- 0.5

# correlation matrix of cell types
xdata <- as.matrix(testdata[, cell_types])
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
    x = rank(testdata[, x_cell]), y = rank(testdata[, y_cell]),
    las = 1, xlab = paste(x_cell, "(Rank)"), ylab = paste(y_cell, "(Rank)"), asp = 1,
    main = paste0("Spearman correlation = ", strong_cor_df$spearman_correlation[i]),
    pch = 21, bg = grey(0.5, 0.3)
  )
  dev.off()
}
