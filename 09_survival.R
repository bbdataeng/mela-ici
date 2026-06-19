# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(paletteer) # v1.7.0
  library(survival) # v3.8-6
})


# Load data ---------------------------------------------------------------

# load random forests data
load("nonsync/08a_RF_metadataset.RData")
# remove unused objects
rm(list = setdiff(ls(), c("alldata", "rf_cv_binary_list", "to_do_list")))

# add missing metadata columns
xx <- readRDS("nonsync/04_clean_data/clean_metadata.rds")
alldata <- lapply(alldata, merge, y = xx, all.x = TRUE)
rm(xx)
alldata <- lapply(alldata, function(x) {
  names(x) <- gsub("\\(days\\)$", "", names(x))
  x
})

# select RF model to use
rf_2use <- "RF5_bin"
rf_2use_idx <- which(to_do_list[[1]]$rf_formula == rf_2use)

# Prepare output folders --------------------------------------------------

# main output folders
outdirs <- c(
  # output folder for complete dataset
  complete = "nonsync/09_survival/complete",
  # output folder for dataset without checkmate067
  nocheckmate067 = "nonsync/09_survival/nocheckmate067"
)

# create directories
for (i in seq_along(outdirs)) {
  if (!dir.exists(outdirs[[i]])) dir.create(outdirs[[i]], recursive = TRUE)
}

# Subset data based on survival availability ------------------------------

# subset data for PFS
pfs_data <- lapply(alldata, function(x) {
  x[!is.na(x$PFS) & !is.na(x$last_followup_status), ]
})
lapply(pfs_data, dim)

# subset data for OS
os_data <- lapply(alldata, function(x) {
  x[!is.na(x$OS) & !is.na(x$last_followup_status), ]
})
lapply(os_data, dim)

# check
identical(pfs_data[[1]]$accession, pfs_data[[2]]$accession)
identical(os_data[[1]]$accession, os_data[[2]]$accession)
# data are identical (they don't include CheckMate 067)
# can use only one
pfs_data <- pfs_data[[1]]
os_data <- os_data[[1]]


# Extract predicted classes -----------------------------------------------

# get predicted classes for PFS
pred_pfs <- lapply(rf_cv_binary_list, function(x) {
  predict(
    x[[rf_2use]],
    newdata = pfs_data
  )
})
stopifnot(length(pred_pfs[[1]]) == nrow(pfs_data))
stopifnot(length(pred_pfs[[2]]) == nrow(pfs_data))
do.call(table, pred_pfs)

# get predicted classes for OS
pred_os <- lapply(rf_cv_binary_list, function(x) {
  predict(
    x[[rf_2use]],
    newdata = os_data
  )
})
stopifnot(length(pred_os[[1]]) == nrow(os_data))
stopifnot(length(pred_os[[2]]) == nrow(os_data))
do.call(table, pred_os)

# compare with real classes
lapply(pred_pfs, table, true = pfs_data$response_2levels)
lapply(pred_os, table, true = os_data$response_2levels)


# Prepare data for survival analysis --------------------------------------

# data for overall survival analysis
survdata_os <- lapply(pred_os, function(x) {
  xx <- data.frame(
    time = os_data$OS,
    # event happened if patient is dead
    event = as.numeric(os_data$last_followup_status == "Dead"),
    real_class = os_data$response_2levels,
    pred_class = x
  )
  xx$surv <- Surv(
    time = xx$time,
    event = xx$event
  )
  xx
})

# data for progression-free survival analysis
survdata_pfs <- lapply(pred_pfs, function(x) {
  # data for progression-free survival
  xx <- data.frame(
    time = pfs_data$PFS,
    real_class = pfs_data$response_2levels,
    pred_class = x
  )
  # event happened unless PFS = OS and patient is alive
  xx$event <- 1
  xx$event[pfs_data$PFS == pfs_data$OS & pfs_data$last_followup_status == "Alive"] <- 0
  xx$surv <- Surv(
    time = xx$time,
    event = xx$event
  )
  xx
})


# Survival analysis -------------------------------------------------------

# fit models
fit_os_real <- lapply(survdata_os, function(x) survfit(surv ~ real_class, data = x))
fit_os_pred <- lapply(survdata_os, function(x) survfit(surv ~ pred_class, data = x))
fit_pfs_real <- lapply(survdata_pfs, function(x) survfit(surv ~ real_class, data = x))
fit_pfs_pred <- lapply(survdata_pfs, function(x) survfit(surv ~ pred_class, data = x))


# test
test_os_real <- lapply(survdata_os, function(x) survdiff(surv ~ real_class, data = x))
test_os_pred <- lapply(survdata_os, function(x) survdiff(surv ~ pred_class, data = x))
test_pfs_real <- lapply(survdata_pfs, function(x) survdiff(surv ~ real_class, data = x))
test_pfs_pred <- lapply(survdata_pfs, function(x) survdiff(surv ~ pred_class, data = x))


# Plot settings -----------------------------------------------------------

resol <- 300

colors_response2 <- c(
  NR = "#D5B3FF",
  R = "#FBD960"
)


# Plot survival curves ----------------------------------------------------

# plotting function
plot_survival <- function(fit, test, main, xlab, ylab, col,
                          lwd = 3, legend_pos = "bottomleft") {
  par(xpd = TRUE)
  plot(
    fit,
    mark.time = FALSE, bty = "o",
    col = col, lwd = lwd,
    xlab = xlab, ylab = ylab,
    main = main, xaxs = "S", yaxs = "i"
  )
  par(xpd = FALSE)
  xx <- data.frame(
    time = fit$time, cond = rep(1:2, times = fit$strata),
    lower = fit$lower, upper = fit$upper
  )
  xx <- xx[complete.cases(xx), ]
  for (i in 1:2) {
    xsub <- subset(xx, cond == i)
    polygon(
      x = c(xsub$time, rev(xsub$time)),
      y = c(xsub$lower, rev(xsub$upper)),
      border = NA, col = adjustcolor(col[i], 0.2)
    )
  }
  points(
    x = fit$time, y = fit$surv,
    pch = ifelse(fit$n.event == 0, 3, NA),
    cex = 1,
    col = col[rep(1:2, times = fit$strata)]
  )
  mtext(
    text = paste0(
      "Log-rank test:\n",
      test$pvalue |>
        insight::format_p(stars = TRUE, decimal_separator = ".", digits = 3)
    ),
    side = 1, adj = 1, line = -1, cex = 0.7
  )
  legend(
    x = legend_pos,
    bty = "n", cex = 0.8, legend = names(col),
    col = col, lwd = lwd, pch = 3, pt.lwd = 1
  )
}

for (x in names(outdirs)) {
  png(file.path(outdirs[x], "survival_plot.png"),
    width = 8 * resol, height = 6 * resol, res = resol
  )
  par(
    mfrow = c(2, 2), las = 1, font.main = 1,
    mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2, 0.7, 0), tcl = -0.3
  )
  plot_survival(
    fit = fit_os_real[[x]],
    test = test_os_real[[x]],
    main = "Overall survival\nstratified by true response",
    xlab = "Overall Survival (days)", ylab = "Probability",
    col = colors_response2, lwd = 3
  )
  plot_survival(
    fit = fit_os_pred[[x]],
    test = test_os_pred[[x]],
    main = "Overall survival\nstratified by predicted response",
    xlab = "Overall Survival (days)", ylab = "Probability",
    col = colors_response2, lwd = 3
  )
  plot_survival(
    fit = fit_pfs_real[[x]],
    test = test_pfs_real[[x]],
    main = "Progression-free survival\nstratified by true response",
    xlab = "Progression-Free Survival (days)", ylab = "Probability",
    col = colors_response2, lwd = 3, legend_pos = "topright"
  )
  plot_survival(
    fit = fit_pfs_pred[[x]],
    test = test_pfs_pred[[x]],
    main = "Progression-free survival\nstratified by predicted response",
    xlab = "Progression-Free Survival (days)", ylab = "Probability",
    col = colors_response2, lwd = 3, legend_pos = "topright"
  )
  dev.off()
}


# Summary statistics ------------------------------------------------------

times_days <- c(1, 2, 4) * 365.25 # get survival rates at 1, 2 and 4 years

for (x in names(outdirs)) {
  sink(file.path(outdirs[x], "OS_true_response.txt"))
  cat("===== Overall survival stratified by true response =====\n\n")
  print(fit_os_real[[x]])
  cat("\n\nSurvival Rates\n")
  print(summary(fit_os_real[[x]], times = times_days))
  sink()

  sink(file.path(outdirs[x], "OS_predicted_response.txt"))
  cat("===== Overall survival stratified by predicted response =====\n")
  cat(paste0("===== Predicted based on ", rf_2use, " =====\n\n"))
  print(fit_os_pred[[x]])
  cat("\n\nSurvival Rates\n")
  print(summary(fit_os_pred[[x]], times = times_days))
  sink()

  sink(file.path(outdirs[x], "PFS_true_response.txt"))
  cat("===== Progression-free survival stratified by true response =====\n\n")
  print(fit_pfs_real[[x]])
  cat("\n\nSurvival Rates\n")
  print(summary(fit_pfs_real[[x]], times = times_days))
  sink()

  sink(file.path(outdirs[x], "PFS_predicted_response.txt"))
  cat("===== Progression-free survival stratified by predicted response =====\n")
  cat(paste0("===== Predicted based on ", rf_2use, " =====\n\n"))
  print(fit_pfs_pred[[x]])
  cat("\n\nSurvival Rates\n")
  print(summary(fit_pfs_pred[[x]], times = times_days))
  sink()
}


# Compare true vs. predicted classes --------------------------------------

data_nr_os <- lapply(survdata_os, function(x) {
  xx <- x[c(which(x$real_class == "NR"), which(x$pred_class == "NR")), ]
  xx$type <- rep(c("true", "predicted"), times = c(sum(x$real_class == "NR"), sum(x$pred_class == "NR")))
  xx[, c("time", "event", "type", "surv")]
})
data_nr_pfs <- lapply(survdata_pfs, function(x) {
  xx <- x[c(which(x$real_class == "NR"), which(x$pred_class == "NR")), ]
  xx$type <- rep(c("true", "predicted"), times = c(sum(x$real_class == "NR"), sum(x$pred_class == "NR")))
  xx[, c("time", "event", "type", "surv")]
})
data_r_os <- lapply(survdata_os, function(x) {
  xx <- x[c(which(x$real_class == "R"), which(x$pred_class == "R")), ]
  xx$type <- rep(c("true", "predicted"), times = c(sum(x$real_class == "R"), sum(x$pred_class == "R")))
  xx[, c("time", "event", "type", "surv")]
})
data_r_pfs <- lapply(survdata_pfs, function(x) {
  xx <- x[c(which(x$real_class == "R"), which(x$pred_class == "R")), ]
  xx$type <- rep(c("true", "predicted"), times = c(sum(x$real_class == "R"), sum(x$pred_class == "R")))
  xx[, c("time", "event", "type", "surv")]
})

# fit models
fit_nr_os <- lapply(data_nr_os, function(x) survfit(surv ~ type, data = x))
fit_nr_pfs <- lapply(data_nr_pfs, function(x) survfit(surv ~ type, data = x))
fit_r_os <- lapply(data_r_os, function(x) survfit(surv ~ type, data = x))
fit_r_pfs <- lapply(data_r_pfs, function(x) survfit(surv ~ type, data = x))

# test
test_nr_os <- lapply(data_nr_os, function(x) survdiff(surv ~ type, data = x))
test_nr_pfs <- lapply(data_nr_pfs, function(x) survdiff(surv ~ type, data = x))
test_r_os <- lapply(data_r_os, function(x) survdiff(surv ~ type, data = x))
test_r_pfs <- lapply(data_r_pfs, function(x) survdiff(surv ~ type, data = x))

# colors
colors_type <- c(
  true = "brown",
  predicted = "darkgreen"
)
for (x in names(outdirs)) {
  png(file.path(outdirs[x], "survival_plot_type.png"),
    width = 8 * resol, height = 6 * resol, res = resol
  )
  par(
    mfrow = c(2, 2), las = 1, font.main = 1,
    mar = c(3.5, 3.5, 3.5, 0.5), mgp = c(2, 0.7, 0), tcl = -0.3
  )
  plot_survival(
    fit = fit_nr_os[[x]],
    test = test_nr_os[[x]],
    main = "Overall survival\nstratified by true and predicted non-responders",
    xlab = "Overall Survival (days)", ylab = "Probability",
    col = colors_type, lwd = 3
  )
  plot_survival(
    fit = fit_r_os[[x]],
    test = test_r_os[[x]],
    main = "Overall survival\nstratified by true and predicted responders",
    xlab = "Overall Survival (days)", ylab = "Probability",
    col = colors_type, lwd = 3
  )
  plot_survival(
    fit = fit_nr_pfs[[x]],
    test = test_pfs_real[[x]],
    main = "Progression-free survival\nstratified by true and predicted non-responders",
    xlab = "Progression-Free Survival (days)", ylab = "Probability",
    col = colors_type, lwd = 3, legend_pos = "topright"
  )
  plot_survival(
    fit = fit_nr_pfs[[x]],
    test = test_pfs_real[[x]],
    main = "Progression-free survival\nstratified by true and predicted responders",
    xlab = "Progression-Free Survival (days)", ylab = "Probability",
    col = colors_type, lwd = 3, legend_pos = "topright"
  )
  dev.off()
}
