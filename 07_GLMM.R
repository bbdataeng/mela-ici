# load("nonsync/07_GLMM.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(paletteer) # v1.7.0
  library(lme4) # v2.0-1
  library(car) # v3.1-5
  source("functions_GLMMs/diagnostic_fcns.r") # diagnostic functions for GLMMs
  source("functions_GLMMs/glmm_stability.r") # stability function for GLMMs
  source("functions_GLMMs/boot_glmm.r") # stability function for GLMMs
  source("functions_GLMMs/fit_function.R") # function to obtain fitted values
  source("plot_glmer_results.R") # function to plot results of binomial GLMM
})

# Prepare output folders --------------------------------------------------

# main output folders
outdirs <- c(
  # output folder for complete dataset
  complete = "nonsync/07_GLMM/complete",
  # output folder for dataset without checkmate067
  nocheckmate067 = "nonsync/07_GLMM/nocheckmate067"
)

# create directories
for (i in seq_along(outdirs)) {
  if (!dir.exists(outdirs[[i]])) dir.create(outdirs[[i]], recursive = TRUE)
}

# prepare instructions
to_do_list <- lapply(outdirs, function(x) {
  # define combinations of GLMMs to be created
  to_do <- expand.grid(
    response_type = "binary",
    age = c(TRUE, FALSE),
    gender = c(TRUE, FALSE)
  )
  # define formula names
  to_do$f <- paste0("F", 1:nrow(to_do))
  # prepare new folder names
  to_do$folder <- character(nrow(to_do))
  for (i in seq_len(nrow(to_do))) {
    to_do_list[[x]]$folder[i] <- file.path(
      x, # main output folder
      paste(to_do$f[i],
        ifelse(to_do$age[i] && to_do$gender[i], "withAge_withGender", ifelse(
          to_do$age[i], "withAge_woGender", ifelse(
            to_do$gender[i], "woAge_withGender", "woAge_woGender"
          )
        )),
        sep = "_"
      )
    )
    # create folders
    if (!dir.exists(to_do_list[[x]]$folder[i])) dir.create(to_do_list[[x]]$folder[i])
  }
  to_do # return object
})


# Load and prepare data ---------------------------------------------------

# load complete data
metadata_complete <- readRDS("nonsync/04_clean_data/clean_metadata.rds")
hed_complete <- readRDS("nonsync/04_clean_data/clean_hed.rds")

# add cibersortx absolute score to metadata
metadata_complete$cibersortx_Absolute_Score <- readRDS("nonsync/04_clean_data/clean_cibersortx.rds")$Absolute_Score

# transform binary response variable into 0s (NR) and 1s (R)
metadata_complete$response <- as.numeric(metadata_complete$response_2levels == "R")

# identify checkmate067 cohort
to_exclude <- which(
  metadata_complete$dataset == "Campbell-2023" &
    metadata_complete$enrichment_protocol == "targeted-mRNA-capture"
)

# keep only necessary metadata columns
metadata_complete <- metadata_complete[, c("accession", "response", "age", "gender", "dataset", "cibersortx_Absolute_Score")]

# remove mean HED column
hed_complete <- hed_complete[, setdiff(names(hed_complete), "HED_mean")]

# merge metadata and HED and create lists of data
alldata <- list(
  complete = cbind(metadata_complete, hed_complete),
  nocheckmate067 = cbind(metadata_complete, hed_complete)[-to_exclude, ]
)
rm(metadata_complete, hed_complete, to_exclude)

# load PCA of cibersortx
pca <- list(
  readRDS("nonsync/06_PCA_multicollinearity/complete/PCA/results_pca.rds"),
  readRDS("nonsync/06_PCA_multicollinearity/nocheckmate067/PCA/results_pca.rds")
)
names(pca) <- names(to_do_list)

# how many PCs to include in the analysis?
# criterion: first n PCs that sum up for 80% of the total variance
n_pcs <- sapply(pca, function(x) {
  cumulative_proportion <- summary(x)[[6]]["Cumulative Proportion", ]
  which(cumulative_proportion >= 0.8) |>
    head(1) |>
    as.numeric()
})
n_pcs

# add PCs to the data
for (i in seq_along(alldata)) {
  alldata[[i]] <- cbind(alldata[[i]], as.data.frame(pca[[i]]$x)[, seq_len(n_pcs[i])])
}


# Prepare test data -------------------------------------------------------

df_data_list <- lapply(seq_along(to_do_list), function(x) {
  dfs <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    vars_to_include <- names(alldata[[x]])
    # exclude age when appropriate
    if (!to_do_list[[x]]$age[i]) vars_to_include <- setdiff(vars_to_include, "age")
    # exclude gender when appropriate
    if (!to_do_list[[x]]$gender[i]) vars_to_include <- setdiff(vars_to_include, "gender")
    # exclude incomplete cases
    xx <- alldata[[x]][, vars_to_include] |>
      na.exclude() |>
      droplevels()
    # z-scale age (if present) and HED data
    vars_to_scale <- grepv("^HED", vars_to_include)
    if (to_do_list[[x]]$age[i]) vars_to_scale <- c("age", vars_to_scale)
    for (xcol in vars_to_scale) {
      xx[, xcol] <- xx[, xcol] |>
        scale() |>
        as.vector()
    }
    # rename z-scaled variables
    names(xx) <- ifelse(
      vars_to_include %in% vars_to_scale,
      paste0("z.", vars_to_include), vars_to_include
    )
    # save data in the respective folder
    write.csv(xx,
      file = file.path(to_do_list[[x]]$folder[i], "data.csv"),
      row.names = FALSE
    )
    # return object
    xx
  })
  names(dfs) <- to_do_list[[x]]$f
  dfs
})
names(df_data_list) <- names(to_do_list)


# Define formulas ---------------------------------------------------------

for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    # define fixed effect predictors
    predictors <- df_data_list[[x]][[i]] |>
      names() |>
      setdiff(y = c("accession", "response", "dataset")) |>
      paste(collapse = " + ")
    # build formula
    to_do_list[[x]]$formula[i] <- paste0(
      "response ~ ", # response variable
      predictors, # fixed effect predictors
      " + (1 | dataset)" # random effect predictors
    )

    # save formula, number of observations, number of datasets in each folder
    n_data <- nrow(df_data_list[[x]][[i]])
    sink(file.path(to_do_list[[x]]$folder[i], "formula_sampleSizes.txt"))
    cat("===== Formula =====", sep = "\n")
    cat(to_do_list[[x]]$formula[i], "\n")
    cat("\n===== Excluded data (NAs present) =====", sep = "\n")
    cat(paste0(nrow(alldata[[x]]) - n_data, " observations"), sep = "\n")
    cat("\n===== Data used for model fitting =====", sep = "\n")
    cat(paste0(n_data, " observations"), sep = "\n")
    cat("\nResponse levels")
    cat(
      "\nNR (=0): ", sum(df_data_list[[x]][[i]]$response == 0), # number of non-respondants
      "\t\tR (=1): ", sum(df_data_list[[x]][[i]]$response == 1) # number of respondants
    )
    cat("\n\nNumber of 'dataset' levels: ", nlevels(df_data_list[[x]][[i]]$dataset))
    print(table(df_data_list[[x]][[i]]$dataset))
    sink()
  }
}


# Fit models --------------------------------------------------------------

# fit models
models <- lapply(seq_along(to_do_list), function(x) {
  mods <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    glmer(
      as.formula(to_do_list[[x]]$formula[i]), # model formula
      data = df_data_list[[x]][[i]], # data
      family = binomial # binomial model
    )
  })
  names(mods) <- to_do_list[[x]]$f
  mods
}) # all models converged but with boundary singular fit message
names(models) <- names(to_do_list)


# export model objects and summaries
for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    # save model objects
    saveRDS(
      object = models[[x]][[i]],
      file = file.path(to_do_list[[x]]$folder[i], "model.rds")
    )
    # write model summary
    print(summary(models[[x]][[i]]), correlation = TRUE) |>
      capture.output(file = file.path(to_do_list[[x]]$folder[i], "model_summary.txt"))
  }
}


# Model diagnostics -------------------------------------------------------

# calculate model stability
stability <- lapply(seq_along(to_do_list), function(x) {
  stab <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    glmm.model.stab(model.res = models[[x]][[i]])
  })
  names(stab) <- to_do_list[[x]]$f
  stab
})
names(stability) <- names(to_do_list)

# export diagnostics
for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    ### collinearity and overdispersion ##
    sink(file.path(to_do_list[[x]]$folder[i], "model_diagnostics.txt"))
    cat("===== VIFs (Variance Inflation Factors) =====
")
    print(vif(models[[x]][[i]]))
    cat("
===== Overdispersion test =====
")
    print(overdisp.test(models[[x]][[i]]))
    sink()

    ### distribution of BLUPs ###
    png(file.path(to_do_list[[x]]$folder[i], "BLUPs_distribution.png"),
      width = 3 * 300, height = 3 * 300, res = 300
    )
    ranef.diagn.plot(models[[x]][[i]])
    dev.off()

    ### model stability results ###
    # table
    stability[[x]][[i]]$summary[, -1] |>
      as.matrix() |>
      round(3) |>
      capture.output(file = file.path(to_do_list[[x]]$folder[i], "model_stability.txt"))

    # plot
    png(file.path(to_do_list[[x]]$folder[i], "model_stability.png"),
      width = 6 * 300, height = 4 * 300, res = 300
    )
    m.stab.plot(stability[[x]][[i]]$summary[, -1])
    dev.off()
  }
}

# Parametric bootstrap ----------------------------------------------------

# run parametric bootstrap
bootstrap <- lapply(seq_along(to_do_list), function(x) {
  boot <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    cat("Running parametric bootstrap for F", i, "...", sep = "")
    xx <- boot.glmm.pred(
      model.res = models[[x]][[i]],
      nboots = 1000, para = TRUE, n.cores = 12, level = 0.95, use = "PC1"
    )
    cat(" done!", sep = "\n")
    xx
  })
  names(boot) <- to_do_list[[x]]$f
  boot
})
names(bootstrap) <- names(to_do_list)

# export bootstrap results
for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    # table
    bootstrap[[x]][[i]]$ci.estimates |>
      round(digits = 3) |> # round to 3 digits
      capture.output(file = file.path(to_do_list[[x]]$folder[i], "bootstrap_CIs.txt"))
    # quick plot
    png(file.path(to_do_list[[x]]$folder[i], "bootstrap_CIs.png"),
      width = 6 * 300, height = 4 * 300, res = 300
    )
    m.stab.plot(bootstrap[[x]][[i]]$ci.estimates)
    dev.off()
    # save bootstrap object
    saveRDS(bootstrap[[x]][[i]], file.path(to_do_list[[x]]$folder[i], "bootstrap_results.rds"))
  }
}

# Model results tables ----------------------------------------------------

# function to get table of model results
get_results_table <- function(model, anova_res, boot_object) {
  restab <- summary(model)$coefficients |> as.data.frame()
  restab$OR <- exp(restab$Estimate)
  restab$lower_CI <- boot_object$ci.estimates$X2.5.
  restab$upper_CI <- boot_object$ci.estimates$X97.5.
  restab$Chisq <- c(NA, anova_res$Chisq)
  restab <- as.matrix(restab) |>
    round(3) |>
    as.data.frame() # round to 3 digits
  restab$p.val <- c( # get p values:
    restab["(Intercept)", "Pr(>|z|)"], # intercept: p-value obtained from z-test
    anova_res$`Pr(>Chisq)`
  ) # predictors: p-value obtained from ANOVA
  restab$p.val <- insight::format_p( # format p-values
    restab$p.val,
    stars = TRUE, name = NULL, digits = 3
  )
  restab <- restab[, c(
    "Estimate", "OR", "lower_CI", "upper_CI", "Std. Error",
    "Chisq", "p.val"
  )]
  return(restab)
}

# get tables of model results
res_tabs <- lapply(seq_along(to_do_list), function(x) {
  xres <- lapply(seq_len(nrow(to_do_list[[x]])), function(i) {
    get_results_table(
      model = models[[x]][[i]], # fitted model
      anova_res = Anova( # LRT for p-values of fixed effect predictors
        mod = models[[x]][[i]], type = "II", test.statistic = "Chisq"
      ),
      boot_object = bootstrap[[x]][[i]] # parametric bootstrap
    )
  })
  names(xres) <- to_do_list[[x]]$f
  xres
})
names(res_tabs) <- names(to_do_list)

# export tables of model results
for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    capture.output(
      res_tabs[[x]][[i]],
      file = file.path(to_do_list[[x]]$folder[i], "model_results_table.txt")
    )
  }
}

# Plots showing PCs -------------------------------------------------------

# plot resolution (pixel per inch)
resol <- 300

for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    for (var_to_plot in paste0("PC", 1:n_pcs[x])) {
      xx_path <- file.path(
        to_do_list[[x]]$folder[i], paste0("fitted_", var_to_plot, ".png")
      )
      xx_n_features <- sum(abs( # number of features with loading >= 0.3
        pca[[x]]$rotation[, var_to_plot]
      ) >= 0.3)
      plot_glmer_fitted_PCA(
        model = models[[x]][[i]], # binomial model fitted with lme4::glmer()
        model_formula = as.formula(to_do_list[[x]]$formula[i]), # formula used to fit the model
        results_table = res_tabs[[x]][[i]], # results_table, output of get_results_table()
        data = df_data_list[[x]][[i]], # data used to fit the model
        response_levels = c("NR", "R"), # levels of response, corresponding to 0 and 1
        bootstrap_object = bootstrap[[x]][[i]], # bootstrap object returned by boot.glmm.pred()
        pca_results = pca[[x]], # PCA results returned by prcomp()
        link = "logit", # link function
        fitted_resolution = 100, # resolution of calculated fitted values
        PC_to_plot = var_to_plot, # which PC should be plotted
        n_top_feautures = xx_n_features, # number of top features to plot
        file_path = xx_path, # optional: path to save the plot
        res_ppi = 300, # resolution (pixels per inch)
        relative_heights = c(2, 1), # relative heights of the two plot sections
        width_in = 6, height_in = 5 # width and height of the plot in inches
      )
    }
  }
}


# Plots showing CIBERSORTx absolute score ---------------------------------

for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    plot_glmer_fitted_continuous(
      model = models[[x]][[i]], # binomial model fitted with lme4::glmer()
      model_formula = as.formula(to_do_list[[x]]$formula[i]), # formula used to fit the model
      results_table = res_tabs[[x]][[i]], # results_table, output of get_results_table()
      data = df_data_list[[x]][[i]], # data used to fit the model
      response_levels = c("NR", "R"), # levels of response, corresponding to 0 and 1
      bootstrap_object = bootstrap[[x]][[i]], # bootstrap object returned by boot.glmm.pred()
      link = "logit", # link function
      covariate_to_plot = "cibersortx_Absolute_Score", # which covariate should be plotted
      fitted_resolution = 100, # resolution of calculated fitted values
      file_path = file.path(to_do_list[[x]]$folder[i], paste0("fitted_cibersortx_absScore.png")), # path to save the plot
      res_ppi = 300, # resolution (pixels per inch)
      width_in = 8, height_in = 4 # width and height of the plot in inches
    )
  }
}


# Plots showing HED data --------------------------------------------------

# plot separately
for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    for (var_to_plot in c("z.HED_locusA", "z.HED_locusB", "z.HED_locusC")) {
      plot_glmer_fitted_continuous(
        model = models[[x]][[i]], # binomial model fitted with lme4::glmer()
        model_formula = as.formula(to_do_list[[x]]$formula[i]), # formula used to fit the model
        results_table = res_tabs[[x]][[i]], # results_table, output of get_results_table()
        data = df_data_list[[x]][[i]], # data used to fit the model
        response_levels = c("NR", "R"), # levels of response, corresponding to 0 and 1
        bootstrap_object = bootstrap[[x]][[i]], # bootstrap object returned by boot.glmm.pred()
        link = "logit", # link function
        covariate_to_plot = var_to_plot, # which covariate should be plotted
        fitted_resolution = 100, # resolution of calculated fitted values
        file_path = file.path(to_do_list[[x]]$folder[i], paste0("fitted_", var_to_plot, ".png")), # path to save the plot
        res_ppi = 300, # resolution (pixels per inch)
        width_in = 7, height_in = 4 # width and height of the plot in inches
      )
    }
  }
}

# plot together
for (x in seq_along(to_do_list)) {
  for (i in seq_len(nrow(to_do_list[[x]]))) {
    plot_glmer_fitted_multiple_continuous(
      model = models[[x]][[i]], # binomial model fitted with lme4::glmer()
      model_formula = as.formula(to_do_list[[x]]$formula[i]), # formula used to fit the model
      results_table = res_tabs[[x]][[i]], # results_table, output of get_results_table()
      data = df_data_list[[x]][[i]], # data used to fit the model
      response_levels = c("NR", "R"), # levels of response, corresponding to 0 and 1
      bootstrap_object = bootstrap[[x]][[i]], # bootstrap object returned by boot.glmm.pred()
      link = "logit", # link function
      covariates_to_plot = c("z.HED_locusA", "z.HED_locusB", "z.HED_locusC"), # which covariate should be plotted
      fitted_resolution = 100, # resolution of calculated fitted values
      file_path = file.path(to_do_list[[x]]$folder[i], paste0("fitted_z.HED_lociABC.png")), # path to save the plot
      res_ppi = 300, # resolution (pixels per inch)
      width_in = 10, height_in = 4 # width and height of the plot in inches
    )
  }
}


# Plots showing age and gender --------------------------------------------

# age
for (x in seq_along(to_do_list)) {
  for (i in which(to_do_list[[x]]$age)) {
    plot_glmer_fitted_continuous(
      model = models[[x]][[i]], # binomial model fitted with lme4::glmer()
      model_formula = as.formula(to_do_list[[x]]$formula[i]), # formula used to fit the model
      results_table = res_tabs[[x]][[i]], # results_table, output of get_results_table()
      data = df_data_list[[x]][[i]], # data used to fit the model
      response_levels = c("NR", "R"), # levels of response, corresponding to 0 and 1
      bootstrap_object = bootstrap[[x]][[i]], # bootstrap object returned by boot.glmm.pred()
      link = "logit", # link function
      covariate_to_plot = "z.age", # which covariate should be plotted
      fitted_resolution = 100, # resolution of calculated fitted values
      file_path = file.path(to_do_list[[x]]$folder[i], paste0("fitted_z.age.png")), # path to save the plot
      res_ppi = 300, # resolution (pixels per inch)
      width_in = 7, height_in = 4 # width and height of the plot in inches
    )
  }
}

# gender
for (x in seq_along(to_do_list)) {
  for (i in which(to_do_list[[x]]$gender)) {
    plot_glmer_fitted_categorical(
      model = models[[x]][[i]], # binomial model fitted with lme4::glmer()
      model_formula = as.formula(to_do_list[[x]]$formula[i]), # formula used to fit the model
      results_table = res_tabs[[x]][[i]], # results_table, output of get_results_table()
      data = df_data_list[[x]][[i]], # data used to fit the model
      bootstrap_object = bootstrap[[x]][[i]], # bootstrap object returned by boot.glmm.pred()
      link = "logit", # link function
      variable_to_plot = "gender", # which variable should be plotted
      file_path = file.path(to_do_list[[x]]$folder[i], paste0("fitted_gender.png")), # path to save the plot
      res_ppi = 300, # resolution (pixels per inch)
      width_in = 7, height_in = 4 # width and height of the plot in inches
    )
  }
}


# Compare estimates -------------------------------------------------------

for (x in seq_along(to_do_list)) {
  # extract model terms names
  xterms <- res_tabs[[x]][[1]] |> rownames()

  # create summary dataframe of model results
  estimates <- lapply(X = seq_len(nrow(to_do_list[[x]])), FUN = function(i) {
    data.frame(
      term = xterms,
      model = to_do_list[[x]]$f[i],
      res_tabs[[x]][[i]][xterms, ],
      row.names = NULL
    )
  }) |> do.call(what = rbind)

  # define range of estimate values
  xrange <- estimates[, c("Estimate", "lower_CI", "upper_CI")] |>
    as.matrix() |>
    range(na.rm = TRUE)
  # prepare data for plotting
  estimates$model <- as.factor(estimates$model)
  levels(estimates$model) <- gsub("^F[0-9]_", "", basename(to_do_list[[x]]$folder))
  estimates$y <- seq_along(xterms)
  estimates$y <- estimates$y + (
    as.numeric(estimates$model) - mean(as.numeric(estimates$model))
  ) * 0.18
  # prepare colors
  xcols <- paletteer_d("RColorBrewer::Dark2", nlevels(estimates$model))
  names(xcols) <- levels(estimates$model)
  # create plot
  resol <- 300
  png(file.path(outdirs[x], "model_estimates.png"),
    width = 8 * resol, height = 6 * resol, res = resol
  )
  par(mar = c(3, 11, 0.5, 0.5), las = 1, xpd = FALSE)
  plot(NULL,
    ylim = c(length(xterms), 0) + 0.5, xlim = c(-max(abs(xrange)), max(abs(xrange))),
    xlab = "", ylab = "", yaxt = "n", bty = "o", yaxs = "i"
  )
  xx <- seq(from = 1, to = length(xterms), by = 2) - 0.5
  rect( # add rectangles for improved visibility
    xleft = rep(par("usr")[1], length(xx)),
    xright = rep(par("usr")[2], length(xx)),
    ybottom = xx, ytop = xx + 1, border = NA, col = grey(0.5, 0.2)
  )
  abline(v = 0, lty = 1, col = "grey40") # add vertical line at 0
  mtext(side = 1, text = "Estimate", line = par("mar")[1] - 1)
  mtext(side = 2, text = "Term", line = par("mar")[2] - 1, las = 0)
  axis(side = 2, at = seq_along(xterms), labels = xterms, tcl = 0)
  segments( # segments spanning over 95% CIs
    x0 = estimates$lower_CI, x1 = estimates$upper_CI,
    y0 = estimates$y, y1 = estimates$y,
    col = xcols[estimates$model], lwd = 2
  )
  points( # points depicting estimates
    x = estimates$Estimate, y = estimates$y,
    pch = 18, cex = 1, col = xcols[estimates$model]
  )
  # add legend
  xx <- levels(estimates$model)
  xx <- paste0(xx, " (F", 1:4, ")")
  legend(
    x = "bottomright", legend = xx,
    lty = 1, pch = 18, lwd = 2, pt.cex = 1, cex = 0.6, col = xcols,
    title = "Model", title.font = 2
  )
  dev.off() # close graphic device
}


# Save image --------------------------------------------------------------

save.image("nonsync/07_GLMM.RData")
