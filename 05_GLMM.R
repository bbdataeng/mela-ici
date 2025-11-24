# load("nonsync/05_GLMM.RData") # run to restore working space

# Load libraries ----------------------------------------------------------
library(paletteer)
library(lme4)
library(car)
source("functions_GLMMs/diagnostic_fcns.r") # diagnostic functions for GLMMs
source("functions_GLMMs/glmm_stability.r") # stability function for GLMMs
source("functions_GLMMs/boot_glmm.r") # stability function for GLMMs
source("functions_GLMMs/fit_function.R") # function to obtain fitted values
source("plot_glmer_results.R") # function to plot results of binomial GLMM

# Prepare new directories -------------------------------------------------

# main output folder
output_folder <- "nonsync/05_GLMM"
if (!dir.exists(output_folder)) dir.create(output_folder)

# define combinations of GLMMs to be created
to_do <- expand.grid(
  response_type = "binary",
  age = c(TRUE, FALSE),
  gender = c(TRUE, FALSE)
)

# define formula names
to_do$f <- paste0("F", 1:nrow(to_do))

# create new folders
to_do$folder <- character(nrow(to_do))
for (i in seq_len(nrow(to_do))) {
  to_do$folder[i] <- file.path(
    output_folder, # main output folder
    paste(to_do$f[i],
      ifelse(to_do$age[i] && to_do$gender[i], "withAge_withGender", ifelse(
        to_do$age[i], "withAge_woGender", ifelse(
          to_do$gender[i], "woAge_withGender", "woAge_woGender"
        )
      )),
      sep = "_"
    )
  )
  if (!dir.exists(to_do$folder[i])) dir.create(to_do$folder[i])
}

# plot resolution (pixel per inch)
resol <- 300


# Load and prepare data ---------------------------------------------------

# load metadata
metadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds") |> as.data.frame()
# load HED data
hed_data <- readRDS("nonsync/01_clean_data/clean_hed.rds") |> as.data.frame()
# load PCA of cibersortx
pca_cibersortx <- readRDS("nonsync/04_multicollinearity/PCA/results_pca.rds")
pcs_data_cibersortx <- pca_cibersortx$x |> as.data.frame()

# how many PCs to include in the analysis?
# criterion: first n PCs that sum up for 80% of the total variance
cumulative_proportion <- summary(pca_cibersortx)[[6]]["Cumulative Proportion", ]
n_pcs <- which(cumulative_proportion >= 0.8) |>
  head(1) |>
  as.numeric()
n_pcs

# which columns have NAs?
metadata |>
  apply(2, function(x) any(is.na(x))) |>
  which() # age, gender, enrichment_protocol
pcs_data_cibersortx |>
  apply(2, function(x) any(is.na(x))) |>
  which() # no missing values
hed_data |>
  apply(2, function(x) any(is.na(x))) |>
  which() # no missing values

# check response levels
str(metadata$response_6levels) # ordered factor (6 levels)
str(metadata$response_2levels) # ordered factor (2 levels)

# transform binary response variable into 0s (NR) and 1s (R)
metadata$response <- as.numeric(metadata$response_2levels == "R")
str(metadata$response)

# keep only necessary metadata columns
metadata <- metadata[, c("accession", "response", "age", "gender", "dataset")]

# remove mean HED column
hed_data <- hed_data[, setdiff(names(hed_data), "HED_mean")]

# merge metadata, HED data, and first PCs
alldata <- cbind(metadata, hed_data, pcs_data_cibersortx[, 1:n_pcs])
rm(metadata, hed_data, pcs_data_cibersortx)



# Inspect data ------------------------------------------------------------

# inspect data tabulations with binary response
table(alldata$response, useNA = "ifany") |> plot() # very good
table(alldata$gender, alldata$response, useNA = "ifany") # not too bad but many NAs in NR

# inspect distribution of continuous variables
hist(alldata$age) # ok
par(mfrow = c(2, 2))
for (xx in grepv("^HED", names(alldata))) hist(alldata[, xx], main = xx) # ok
par(mfrow = c(2, 3))
for (xx in grepv("^PC", names(alldata))) hist(alldata[, xx], main = xx)
# mostly nicely distributed

# age and HED data should be z-transformed for better model convergence and interpretation


# Prepare test data -------------------------------------------------------

# define new dataframe for each row of to_do
to_do$df_object <- paste("df", "data", to_do$f, sep = "_")

# create new dataframes
for (i in seq_len(nrow(to_do))) {
  vars_to_include <- names(alldata)
  # exclude age when appropriate
  if (!to_do$age[i]) vars_to_include <- setdiff(vars_to_include, "age")
  # exclude gender when appropriate
  if (!to_do$gender[i]) vars_to_include <- setdiff(vars_to_include, "gender")
  # exclude incomplete cases
  xx <- alldata[, vars_to_include] |>
    na.exclude() |>
    droplevels()
  # z-scale age (if present) and HED data
  vars_to_scale <- grepv("^HED", vars_to_include)
  if (to_do$age[i]) vars_to_scale <- c("age", vars_to_scale)
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
  # assign new object
  assign(x = to_do$df_object[i], value = xx)
}
rm(xx, vars_to_include, vars_to_scale, xcol)
ls(pattern = "^df") # new dataframes created

# save data in each folder
for (i in seq_len(nrow(to_do))) {
  write.csv(x = get(to_do$df_object[i]), file = file.path(
    to_do$folder[i], "data.csv"
  ), row.names = FALSE)
}


# Define formulas ---------------------------------------------------------

# define formulas
for (i in seq_len(nrow(to_do))) {
  # define fixed effect predictors
  predictors <- get(to_do$df_object[i]) |>
    names() |>
    setdiff(y = c("accession", "response", "dataset")) |>
    paste(collapse = " + ")
  # build formula
  to_do$formula[i] <- paste0(
    "response ~ ", # response variable
    predictors, # fixedeffect predictors
    " + (1 | dataset)" # random effect predictors
  )
}

# save formula, number of observations, number of datasets in each folder
for (i in seq_len(nrow(to_do))) {
  n_data <- nrow(get(to_do$df_object[i]))
  sink(file.path(to_do$folder[i], "formula_sampleSizes.txt"))
  cat("===== Formula =====", sep = "\n")
  cat(to_do$formula[i], "\n")
  cat("\n===== Excluded data (NAs present) =====", sep = "\n")
  cat(paste0(nrow(alldata) - n_data, " observations"), sep = "\n")
  cat("\n===== Data used for model fitting =====", sep = "\n")
  cat(paste0(n_data, " observations"), sep = "\n")
  cat("\nResponse levels")
  cat(
    "\nNR (=0): ", sum(get(to_do$df_object[i])$response == 0), # number of non-respondants
    "\t\tR (=1): ", sum(get(to_do$df_object[i])$response == 1) # number of respondants
  )
  cat("\n\nNumber of 'dataset' levels: ", nlevels(get(to_do$df_object[i])$dataset))
  print(table(get(to_do$df_object[i])$dataset))
  sink()
}



# Fit models --------------------------------------------------------------

# prepare names of model objects
to_do$model <- paste("model", to_do$f, sep = "_")

# fit models
for (i in seq_len(nrow(to_do))) {
  glmer(
    as.formula(to_do$formula[i]), # model formula
    data = get(to_do$df_object[i]), # data
    family = binomial # binomial model
  ) |> assign(x = to_do$model[i]) # assign new object
}
# some models converged but with boundary singular fit message
ls(pattern = "^model") # new model objects

# export model objects and summaries
for (i in seq_len(nrow(to_do))) {
  # save model objects
  saveRDS(
    object = get(to_do$model[i]),
    file = file.path(to_do$folder[i], "model.rds")
  )
  # write model summary
  print(summary(get(to_do$model[i])), correlation = TRUE) |>
    capture.output(file = file.path(to_do$folder[i], "model_summary.txt"))
}


# Model diagnostics -------------------------------------------------------

# define names of new objects for model stability results
to_do$stability <- paste("stability", to_do$f, sep = "_")

# calculate model stability
for (i in seq_len(nrow(to_do))) {
  glmm.model.stab(
    model.res = get(to_do$model[i])
  ) |> assign(x = to_do$stability[i]) # assign new object
}

for (i in seq_len(nrow(to_do))) {
  ### collinearity and overdispersion ##
  sink(file.path(to_do$folder[i], "model_diagnostics.txt"))
  cat("===== VIFs (Variance Inflation Factors) =====
")
  print(vif(get(to_do$model[i])))
  cat("
===== Overdispersion test =====
")
  print(overdisp.test(get(to_do$model[i])))
  sink()

  ### distribution of BLUPs ###
  png(file.path(to_do$folder[i], "BLUPs_distribution.png"),
    width = 3 * 300, height = 3 * 300, res = 300
  )
  ranef.diagn.plot(get(to_do$model[i]))
  dev.off()

  ### model stability results ###
  # table
  get(to_do$stability[i])$summary[, -1] |>
    as.matrix() |>
    round(3) |>
    capture.output(file = file.path(to_do$folder[i], "model_stability.txt"))

  # plot
  png(file.path(to_do$folder[i], "model_stability.png"),
    width = 6 * 300, height = 4 * 300, res = 300
  )
  m.stab.plot(get(to_do$stability[i])$summary[, -1])
  dev.off()
}

# Parametric bootstrap ----------------------------------------------------

# define names of new objects for bootstrap results
to_do$bootstrap <- paste("boostrap", to_do$f, sep = "_")

# run parametric bootstrap
for (i in seq_len(nrow(to_do))) {
  cat("Running parametric bootstrap for F", i, "...", sep = "")
  boot.glmm.pred(
    model.res = get(to_do$model[i]),
    nboots = 1000, para = TRUE, n.cores = 6, level = 0.95, use = "PC1"
  ) |> assign(x = to_do$bootstrap[i]) # assign new object
  cat(" done!", sep = "\n")
}

# export bootstrap results
for (i in seq_len(nrow(to_do))) {
  # table
  get(to_do$bootstrap[i])$ci.estimates |>
    round(digits = 3) |> # round to 3 digits
    capture.output(file = file.path(to_do$folder[i], "bootstrap_CIs.txt"))
  # quick plot
  png(file.path(to_do$folder[i], "bootstrap_CIs.png"),
    width = 6 * 300, height = 4 * 300, res = 300
  )
  m.stab.plot(get(to_do$bootstrap[i])$ci.estimates)
  dev.off()
  # save bootstrap object
  saveRDS(get(to_do$bootstrap[i]), file.path(to_do$folder[i], "bootstrap_results.rds"))
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

# define names of new objects for model results tables
to_do$restab <- paste("restab", to_do$f, sep = "_")

# get results table and save it
for (i in seq_len(nrow(to_do))) {
  # get table of model results
  get_results_table(
    model = get(to_do$model[i]), # fitted model
    anova_res = Anova( # LRT for p-values of fixed effect predictors
      mod = get(to_do$model[i]), type = "II", test.statistic = "Chisq"
    ),
    boot_object = get(to_do$bootstrap[i]) # parametric bootstrap
  ) |> assign(x = to_do$restab[i]) # assign new object

  # save it
  capture.output(
    get(to_do$restab[i]),
    file = file.path(to_do$folder[i], "model_results_table.txt")
  )
}



# Plots showing PCs -------------------------------------------------------

for (i in seq_len(nrow(to_do))) {
  for (var_to_plot in paste0("PC", 1:n_pcs)) {
    xx_path <- file.path(
      to_do$folder[i], paste0("fitted_", var_to_plot, ".png")
    )
    plot_glmer_fitted_PCA(
      model = get(to_do$model[i]), # binomial model fitted with lme4::glmer()
      model_formula = as.formula(to_do$formula[i]), # formula used to fit the model
      results_table = get(to_do$restab[i]), # results_table, output of get_results_table()
      data = get(to_do$df_object[i]), # data used to fit the model
      response_levels = c("NR", "R"), # levels of response, corresponding to 0 and 1
      bootstrap_object = get(to_do$bootstrap[i]), # bootstrap object returned by boot.glmm.pred()
      pca_results = pca_cibersortx, # PCA results returned by prcomp()
      link = "logit", # link function
      fitted_resolution = 100, # resolution of calculated fitted values
      PC_to_plot = var_to_plot, # which PC should be plotted
      n_top_feautures = 10, # optionally, an integer determining the number of top features to plot
      file_path = xx_path, # optional: path to save the plot
      res_ppi = 300, # resolution (pixels per inch)
      relative_heights = c(2, 1), # relative heights of the two plot sections
      width_in = 6, height_in = 5 # width and height of the plot in inches
    )
  }
}

# Plots showing HED data --------------------------------------------------

for (i in seq_len(nrow(to_do))) {
  for (var_to_plot in c("z.HED_locusA", "z.HED_locusB", "z.HED_locusC")) {
    plot_glmer_fitted_continuous(
      model = get(to_do$model[i]), # binomial model fitted with lme4::glmer()
      model_formula = as.formula(to_do$formula[i]), # formula used to fit the model
      results_table = get(to_do$restab[i]), # results_table, output of get_results_table()
      data = get(to_do$df_object[i]), # data used to fit the model
      response_levels = c("NR", "R"), # levels of response, corresponding to 0 and 1
      bootstrap_object = get(to_do$bootstrap[i]), # bootstrap object returned by boot.glmm.pred()
      link = "logit", # link function
      covariate_to_plot = var_to_plot, # which covariate should be plotted
      fitted_resolution = 100, # resolution of calculated fitted values
      file_path = file.path(to_do$folder[i], paste0("fitted_", var_to_plot, ".png")), # path to save the plot
      res_ppi = 300, # resolution (pixels per inch)
      width_in = 7, height_in = 4 # width and height of the plot in inches
    )
  }
}


# Plots showing age and gender --------------------------------------------

# age
for (i in which(to_do$age)) {
  plot_glmer_fitted_continuous(
    model = get(to_do$model[i]), # binomial model fitted with lme4::glmer()
    model_formula = as.formula(to_do$formula[i]), # formula used to fit the model
    results_table = get(to_do$restab[i]), # results_table, output of get_results_table()
    data = get(to_do$df_object[i]), # data used to fit the model
    response_levels = c("NR", "R"), # levels of response, corresponding to 0 and 1
    bootstrap_object = get(to_do$bootstrap[i]), # bootstrap object returned by boot.glmm.pred()
    link = "logit", # link function
    covariate_to_plot = "z.age", # which covariate should be plotted
    fitted_resolution = 100, # resolution of calculated fitted values
    file_path = file.path(to_do$folder[i], paste0("fitted_z.age.png")), # path to save the plot
    res_ppi = 300, # resolution (pixels per inch)
    width_in = 7, height_in = 4 # width and height of the plot in inches
  )
}

# gender
for (i in which(to_do$gender)) {
  plot_glmer_fitted_categorical(
    model = get(to_do$model[i]), # binomial model fitted with lme4::glmer()
    model_formula = as.formula(to_do$formula[i]), # formula used to fit the model
    results_table = get(to_do$restab[i]), # results_table, output of get_results_table()
    data = get(to_do$df_object[i]), # data used to fit the model
    bootstrap_object = get(to_do$bootstrap[i]), # bootstrap object returned by boot.glmm.pred()
    link = "logit", # link function
    variable_to_plot = "gender", # which variable should be plotted
    file_path = file.path(to_do$folder[i], paste0("fitted_gender.png")), # path to save the plot
    res_ppi = 300, # resolution (pixels per inch)
    width_in = 7, height_in = 4 # width and height of the plot in inches
  )
}


# Compare estimates -------------------------------------------------------

# extract model terms names
xterms <- get(to_do$restab[1]) |> rownames()

# create summary dataframe of model results
estimates <- lapply(X = seq_len(nrow(to_do)), FUN = function(i) {
  data.frame(
    term = xterms,
    model = to_do$f[i],
    get(to_do$restab[i])[xterms, ],
    row.names = NULL
  )
}) |> do.call(what = rbind)

# define range of estimate values
xrange <- estimates[, c("Estimate", "lower_CI", "upper_CI")] |>
  as.matrix() |>
  range(na.rm = TRUE)
# prepare data for plotting
estimates$model <- as.factor(estimates$model)
levels(estimates$model) <- gsub("^F[0-9]_", "", basename(to_do$folder))
estimates$y <- seq_along(xterms)
estimates$y <- estimates$y + (
  as.numeric(estimates$model) - mean(as.numeric(estimates$model))
) * 0.18
# prepare colors
xcols <- paletteer_d("RColorBrewer::Dark2", nlevels(estimates$model))
names(xcols) <- levels(estimates$model)
# create plot
resol <- 300
png(file.path(output_folder, "model_estimates.png"),
  width = 8 * resol, height = 6 * resol, res = resol
)
par(mar = c(3, 7, 0.5, 0.5), las = 1, xpd = FALSE)
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


# Save image --------------------------------------------------------------

save.image("nonsync/05_GLMM.RData")
