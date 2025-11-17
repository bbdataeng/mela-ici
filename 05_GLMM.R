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

# folders for different model formulas
folder_all_predictors <- file.path(output_folder, "all_predictors")
folder_w0_age <- file.path(output_folder, "w0_age")
folder_w0_gender <- file.path(output_folder, "w0_gender")
folder_w0_age_gender <- file.path(output_folder, "w0_age_gender")
for (xx in ls(pattern = "^folder")) { # create folders
  if (!dir.exists(get(xx))) dir.create(get(xx))
}

# Load and prepare data ---------------------------------------------------
metadata <- readRDS("nonsync/01_clean_data/clean_metadata.rds") |> as.data.frame()

# load HED data
hed_data <- readRDS("nonsync/01_clean_data/clean_hed.rds") |> as.data.frame()

# load PCs of cibersortx data
pca_cibersortx <- readRDS("nonsync/04_multicollinearity/PCA/results_pca.rds")
pcs_data_cibersortx <- pca_cibersortx$x |> as.data.frame()

# how many PCs to include in the analysis?
# criterion: first n PCs that sum up for 80% of the total variance
cumulative_proportion <- summary(pca_cibersortx)[[6]]["Cumulative Proportion", ]
n_pcs <- which(cumulative_proportion >= 0.8) |>
  head(1) |>
  as.numeric()
n_pcs

# merge metadata, HED data, and first PCs
testdata <- cbind(metadata, hed_data, pcs_data_cibersortx[, 1:n_pcs])

# which columns have NAs?
testdata |>
  apply(2, function(x) any(is.na(x))) |>
  which() # age, gender, enrichment_protocol

# check response levels
str(testdata$response) # ordered factor (6 levels)
str(testdata$response_group) # ordered factor (2 levels)
testdata$response_group <- factor( # response_group non-ordered
  testdata$response_group,
  ordered = FALSE
)

# remove column enrichment protocol (highly collinear with dataset, NAs present)
testdata <- testdata[, names(testdata) != "enrichment_protocol"]

# remove columns not included in the analysis
testdata <- testdata[, setdiff(names(testdata), c(
  "accession", "patient_id", "sample_name", "response", "treatment",
  "biopsy_time", "HED_mean"
))]


# Get subsets without NAs for specific variables --------------------------

# data without column "gender"
data_w0_gender <- testdata[, setdiff(names(testdata), "gender")]

# data without column "age"
data_w0_age <- testdata[, setdiff(names(testdata), "age")]

# data without columns "age" and "gender"
data_w0_age_gender <- testdata[, setdiff(names(testdata), c("age", "gender"))]

# keep only complete case in each dataset
data_all <- testdata[complete.cases(testdata), ] |> droplevels()
data_w0_gender <- data_w0_gender[complete.cases(data_w0_gender), ] |> droplevels()
data_w0_age <- data_w0_age[complete.cases(data_w0_age), ] |> droplevels()
data_w0_age_gender <- data_w0_age_gender[complete.cases(data_w0_age_gender), ] |> droplevels()


# check which datasets are left
data_all$dataset |>
  table() |>
  addmargins() # 106 cases, 3 datasets
data_w0_gender$dataset |>
  table() |>
  addmargins() # 116 cases, 4 datasets
data_w0_age$dataset |>
  table() |>
  addmargins() # 106 cases, 3 datasets
data_w0_age_gender$dataset |>
  table() |>
  addmargins() # 165 cases, 5 datasets


# Prepare test data for LMM -----------------------------------------------

# inspect data tabulations with 2-level response
table(testdata$response_group, useNA = "ifany") |> plot() # very good
table(testdata$gender, testdata$response_group, useNA = "ifany") # not too bad but many NAs in NR

# inspect distribution of continuous variables
hist(testdata$age) # ok
par(mfrow = c(2, 2))
for (xx in grepv("^HED", names(testdata))) hist(testdata[, xx], main = xx) # ok
par(mfrow = c(2, 3))
for (xx in grepv("^PC", names(testdata))) hist(testdata[, xx], main = xx)
# nicely distributed

# z-scale numeric predictors (age)
data_all$z.age <- data_all$age |>
  scale() |>
  as.vector()
data_w0_gender$z.age <- data_w0_gender$age |>
  scale() |>
  as.vector()
for (xcol in grepv("^HED", names(testdata))) {
  data_all[, paste0("z.", xcol)] <- as.vector(scale(data_all[, xcol]))
  data_w0_age[, paste0("z.", xcol)] <- as.vector(scale(data_w0_age[, xcol]))
  data_w0_gender[, paste0("z.", xcol)] <- as.vector(scale(data_w0_gender[, xcol]))
  data_w0_age_gender[, paste0("z.", xcol)] <- as.vector(scale(data_w0_age_gender[, xcol]))
}


# Fit binomial GLMMs ------------------------------------------------------

# define model formula and write them as txt
model_formula_all <- paste(
  "response_group ~ z.age + gender + z.HED_locusA + z.HED_locusB + z.HED_locusC",
  paste(grepv("^PC", names(testdata)), collapse = " + "),
  "(1 | dataset)",
  sep = " + "
) |> as.formula()
model_formula_w0_age <- paste(
  "response_group ~ gender + z.HED_locusA + z.HED_locusB + z.HED_locusC",
  paste(grepv("^PC", names(testdata)), collapse = " + "),
  "(1 | dataset)",
  sep = " + "
) |> as.formula()
model_formula_w0_gender <- paste(
  "response_group ~ z.age + z.HED_locusA + z.HED_locusB + z.HED_locusC",
  paste(grepv("^PC", names(testdata)), collapse = " + "),
  "(1 | dataset)",
  sep = " + "
) |> as.formula()
model_formula_w0_age_gender <- paste(
  "response_group ~ z.HED_locusA + z.HED_locusB + z.HED_locusC",
  paste(grepv("^PC", names(testdata)), collapse = " + "),
  "(1 | dataset)",
  sep = " + "
) |> as.formula()

# fit models
full_bin_glmm_all <- glmer(model_formula_all, data = data_all, family = binomial)
# model converged but with boundary singular fit message
summary(full_bin_glmm_all)$varcor
full_bin_glmm_w0_age <- glmer(model_formula_w0_age, data = data_w0_age, family = binomial)
# model converged but with boundary singular fit message
summary(full_bin_glmm_w0_age)$varcor
full_bin_glmm_w0_gender <- glmer(model_formula_w0_gender, data = data_w0_gender, family = binomial)
# model converged but with boundary singular fit message
summary(full_bin_glmm_w0_gender)$varcor
full_bin_glmm_w0_age_gender <- glmer(model_formula_w0_age_gender, data = data_w0_age_gender, family = binomial)
# model converged
summary(full_bin_glmm_w0_age_gender)$varcor

# export models
saveRDS(full_bin_glmm_all, file.path(folder_all_predictors, "model.rds"))
saveRDS(full_bin_glmm_w0_age, file.path(folder_w0_age, "model.rds"))
saveRDS(full_bin_glmm_w0_gender, file.path(folder_w0_gender, "model.rds"))
saveRDS(full_bin_glmm_w0_age_gender, file.path(folder_w0_age_gender, "model.rds"))

# write model summaries as txt
print(summary(full_bin_glmm_all), correlation = TRUE) |>
  capture.output(file = file.path(folder_all_predictors, "model_summary.txt"))
print(summary(full_bin_glmm_w0_age), correlation = TRUE) |>
  capture.output(file = file.path(folder_w0_age, "model_summary.txt"))
print(summary(full_bin_glmm_w0_gender), correlation = TRUE) |>
  capture.output(file = file.path(folder_w0_gender, "model_summary.txt"))
print(summary(full_bin_glmm_w0_age_gender), correlation = TRUE) |>
  capture.output(file = file.path(folder_w0_age_gender, "model_summary.txt"))



# Model diagnostics -------------------------------------------------------

### collinearity and overdispersion ###

sink(file.path(folder_all_predictors, "model_diagnostics.txt"))
cat("===== VIFs (Variance Inflation Factors) =====
")
print(vif(full_bin_glmm_all))
cat("
===== Overdispersion test =====
")
print(overdisp.test(full_bin_glmm_all))
sink()

sink(file.path(folder_w0_age, "model_diagnostics.txt"))
cat("===== VIFs (Variance Inflation Factors) =====
")
print(vif(full_bin_glmm_w0_age))
cat("
===== Overdispersion test =====
")
print(overdisp.test(full_bin_glmm_w0_age))
sink()

sink(file.path(folder_w0_gender, "model_diagnostics.txt"))
cat("===== VIFs (Variance Inflation Factors) =====
")
print(vif(full_bin_glmm_w0_gender))
cat("
===== Overdispersion test =====
")
print(overdisp.test(full_bin_glmm_w0_gender))
sink()

sink(file.path(folder_w0_age_gender, "model_diagnostics.txt"))
cat("===== VIFs (Variance Inflation Factors) =====
")
print(vif(full_bin_glmm_w0_age_gender))
cat("
===== Overdispersion test =====
")
print(overdisp.test(full_bin_glmm_w0_age_gender))
sink()

###  normal distribution of BLUPs ###
ranef.diagn.plot(full_bin_glmm_all)
ranef.diagn.plot(full_bin_glmm_w0_age)
ranef.diagn.plot(full_bin_glmm_w0_gender)
ranef.diagn.plot(full_bin_glmm_w0_age_gender)

# model stability
stability_full_bin_glmm_all <- glmm.model.stab(model.res = full_bin_glmm_all)
table(stability_full_bin_glmm_all$detailed$lme4.warnings) # 2 singular fits
stability_full_bin_glmm_all$summary[, -1] |>
  as.matrix() |>
  round(3) |>
  capture.output(file = file.path(folder_all_predictors, "model_stability.txt"))
png(file.path(folder_all_predictors, "model_stability.png"),
  width = 6 * 300, height = 4 * 300, res = 300
)
m.stab.plot(stability_full_bin_glmm_all$summary[, -1])
dev.off()

stability_full_bin_glmm_w0_age <- glmm.model.stab(model.res = full_bin_glmm_w0_age)
table(stability_full_bin_glmm_w0_age$detailed$lme4.warnings) # 2 singular fits
stability_full_bin_glmm_w0_age$summary[, -1] |>
  as.matrix() |>
  round(3) |>
  capture.output(file = file.path(folder_w0_age, "model_stability.txt"))
png(file.path(folder_w0_age, "model_stability.png"),
  width = 6 * 300, height = 4 * 300, res = 300
)
m.stab.plot(stability_full_bin_glmm_w0_age$summary[, -1])
dev.off()

stability_full_bin_glmm_w0_gender <- glmm.model.stab(model.res = full_bin_glmm_w0_gender)
table(stability_full_bin_glmm_w0_gender$detailed$lme4.warnings) # 4 singular fits
stability_full_bin_glmm_w0_gender$summary[, -1] |>
  as.matrix() |>
  round(3) |>
  capture.output(file = file.path(folder_w0_gender, "model_stability.txt"))
png(file.path(folder_w0_gender, "model_stability.png"),
  width = 6 * 300, height = 4 * 300, res = 300
)
m.stab.plot(stability_full_bin_glmm_w0_gender$summary[, -1])
dev.off()

stability_full_bin_glmm_w0_age_gender <- glmm.model.stab(model.res = full_bin_glmm_w0_age_gender)
table(stability_full_bin_glmm_w0_age_gender$detailed$lme4.warnings) # 1 singular fit
stability_full_bin_glmm_w0_age_gender$summary[, -1] |>
  as.matrix() |>
  round(3) |>
  capture.output(file = file.path(folder_w0_age_gender, "model_stability.txt"))
png(file.path(folder_w0_age_gender, "model_stability.png"),
  width = 6 * 300, height = 4 * 300, res = 300
)
m.stab.plot(stability_full_bin_glmm_w0_age_gender$summary[, -1])
dev.off()


# Parametric bootstrap ----------------------------------------------------

# model with all predictors
boot_full_bin_glmm_all <- boot.glmm.pred(
  model.res = full_bin_glmm_all, nboots = 1000,
  para = TRUE, n.cores = 6, level = 0.95, use = "PC1"
)
round(boot_full_bin_glmm_all$ci.estimates, 3) |>
  capture.output(file = file.path(folder_all_predictors, "bootstrap_CIs.txt"))
png(file.path(folder_all_predictors, "bootstrap_CIs.png"),
  width = 6 * 300, height = 4 * 300, res = 300
)
m.stab.plot(boot_full_bin_glmm_all$ci.estimates)
dev.off()
saveRDS(boot_full_bin_glmm_all, file.path(folder_all_predictors, "bootstrap_results.rds"))

# model without age
boot_full_bin_glmm_w0_age <- boot.glmm.pred(
  model.res = full_bin_glmm_w0_age, nboots = 1000,
  para = TRUE, n.cores = 6, level = 0.95, use = "PC1"
)
round(boot_full_bin_glmm_w0_age$ci.estimates, 3) |>
  capture.output(file = file.path(folder_w0_age, "bootstrap_CIs.txt"))
png(file.path(folder_w0_age, "bootstrap_CIs.png"),
  width = 6 * 300, height = 4 * 300, res = 300
)
m.stab.plot(boot_full_bin_glmm_w0_age$ci.estimates)
dev.off()
saveRDS(boot_full_bin_glmm_w0_age, file.path(folder_w0_age, "bootstrap_results.rds"))

# model without gender
boot_full_bin_glmm_w0_gender <- boot.glmm.pred(
  model.res = full_bin_glmm_w0_gender, nboots = 1000,
  para = TRUE, n.cores = 6, level = 0.95, use = "PC1"
)
round(boot_full_bin_glmm_w0_gender$ci.estimates, 3) |>
  capture.output(file = file.path(folder_w0_gender, "bootstrap_CIs.txt"))
png(file.path(folder_w0_gender, "bootstrap_CIs.png"),
  width = 6 * 300, height = 4 * 300, res = 300
)
m.stab.plot(boot_full_bin_glmm_w0_gender$ci.estimates)
dev.off()
saveRDS(boot_full_bin_glmm_w0_gender, file.path(folder_w0_gender, "bootstrap_results.rds"))

# model without age and gender
boot_full_bin_glmm_w0_age_gender <- boot.glmm.pred(
  model.res = full_bin_glmm_w0_age_gender, nboots = 1000,
  para = TRUE, n.cores = 6, level = 0.95, use = "PC1"
)
round(boot_full_bin_glmm_w0_age_gender$ci.estimates, 3) |>
  capture.output(file = file.path(folder_w0_age_gender, "bootstrap_CIs.txt"))
png(file.path(folder_w0_age_gender, "bootstrap_CIs.png"),
  width = 6 * 300, height = 4 * 300, res = 300
)
m.stab.plot(boot_full_bin_glmm_w0_age_gender$ci.estimates)
dev.off()
saveRDS(boot_full_bin_glmm_w0_age_gender, file.path(folder_w0_age_gender, "bootstrap_results.rds"))



# Model results tables ----------------------------------------------------

anova_all <- Anova(full_bin_glmm_all, type = "II", test.statistic = "Chisq")
anova_w0_age <- Anova(full_bin_glmm_w0_age, type = "II", test.statistic = "Chisq")
anova_w0_gender <- Anova(full_bin_glmm_w0_gender, type = "II", test.statistic = "Chisq")
anova_w0_age_gender <- Anova(full_bin_glmm_w0_age_gender, type = "II", test.statistic = "Chisq")

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

restab_all <- get_results_table(
  full_bin_glmm_all, anova_all, boot_full_bin_glmm_all
)
restab_w0_age <- get_results_table(
  full_bin_glmm_w0_age, anova_w0_age, boot_full_bin_glmm_w0_age
)
restab_w0_gender <- get_results_table(
  full_bin_glmm_w0_gender, anova_w0_gender, boot_full_bin_glmm_w0_gender
)
restab_w0_age_gender <- get_results_table(
  full_bin_glmm_w0_age_gender, anova_w0_age_gender, boot_full_bin_glmm_w0_age_gender
)

capture.output(
  restab_all,
  file = file.path(folder_all_predictors, "model_results_table.txt")
)
capture.output(
  restab_w0_age,
  file = file.path(folder_w0_age, "model_results_table.txt")
)
capture.output(
  restab_w0_gender,
  file = file.path(folder_w0_gender, "model_results_table.txt")
)
capture.output(
  restab_w0_age_gender,
  file = file.path(folder_w0_age_gender, "model_results_table.txt")
)




# Plots showing PCs -------------------------------------------------------

for (var_to_plot in paste0("PC", 1:n_pcs)) {
  xx_path <- file.path(
    folder_all_predictors, paste0("fitted_", var_to_plot, ".png")
  )
  plot_glmer_fitted_PCA(
    model = full_bin_glmm_all, # binomial model fitted with lme4::glmer()
    model_formula = model_formula_all, # formula used to fit the model
    results_table = restab_all, # results_table, output of get_results_table()
    data = data_all, # data used to fit the model
    bootstrap_object = boot_full_bin_glmm_all, # bootstrap object returned by boot.glmm.pred()
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
for (var_to_plot in paste0("PC", 1:n_pcs)) {
  xx_path <- file.path(
    folder_w0_age_gender, paste0("fitted_", var_to_plot, ".png")
  )
  plot_glmer_fitted_PCA(
    model = full_bin_glmm_w0_age_gender, # binomial model fitted with lme4::glmer()
    model_formula = model_formula_w0_age_gender, # formula used to fit the model
    results_table = restab_w0_age_gender, # results_table, output of get_results_table()
    data = data_w0_age_gender, # data used to fit the model
    bootstrap_object = boot_full_bin_glmm_w0_age_gender, # bootstrap object returned by boot.glmm.pred()
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
for (var_to_plot in paste0("PC", 1:n_pcs)) {
  xx_path <- file.path(
    folder_w0_age, paste0("fitted_", var_to_plot, ".png")
  )
  plot_glmer_fitted_PCA(
    model = full_bin_glmm_w0_age, # binomial model fitted with lme4::glmer()
    model_formula = model_formula_w0_age, # formula used to fit the model
    results_table = restab_w0_age, # results_table, output of get_results_table()
    data = data_w0_age, # data used to fit the model
    bootstrap_object = boot_full_bin_glmm_w0_age, # bootstrap object returned by boot.glmm.pred()
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
for (var_to_plot in paste0("PC", 1:n_pcs)) {
  xx_path <- file.path(
    folder_w0_gender, paste0("fitted_", var_to_plot, ".png")
  )
  plot_glmer_fitted_PCA(
    model = full_bin_glmm_w0_gender, # binomial model fitted with lme4::glmer()
    model_formula = model_formula_w0_gender, # formula used to fit the model
    results_table = restab_w0_gender, # results_table, output of get_results_table()
    data = data_w0_gender, # data used to fit the model
    bootstrap_object = boot_full_bin_glmm_w0_gender, # bootstrap object returned by boot.glmm.pred()
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

# Plots showing HED data --------------------------------------------------
for (xvar in c("z.HED_locusA", "z.HED_locusB", "z.HED_locusC")) {
  plot_glmer_fitted_continuous(
    model = full_bin_glmm_all, # binomial model fitted with lme4::glmer()
    model_formula = model_formula_all, # formula used to fit the model
    results_table = restab_all, # results_table, output of get_results_table()
    data = data_all, # data used to fit the model
    bootstrap_object = boot_full_bin_glmm_all, # bootstrap object returned by boot.glmm.pred()
    link = "logit", # link function
    covariate_to_plot = xvar, # which covariate should be plotted
    fitted_resolution = 100, # resolution of calculated fitted values
    file_path = file.path(folder_all_predictors, paste0("fitted_", xvar, ".png")), # optional: path to save the plot
    res_ppi = 300, # resolution (pixels per inch)
    width_in = 7, height_in = 4 # width and height of the plot in inches
  )
  plot_glmer_fitted_continuous(
    model = full_bin_glmm_w0_gender, # binomial model fitted with lme4::glmer()
    model_formula = model_formula_w0_gender, # formula used to fit the model
    results_table = restab_w0_gender, # results_table, output of get_results_table()
    data = data_w0_gender, # data used to fit the model
    bootstrap_object = boot_full_bin_glmm_w0_gender, # bootstrap object returned by boot.glmm.pred()
    link = "logit", # link function
    covariate_to_plot = xvar, # which covariate should be plotted
    fitted_resolution = 100, # resolution of calculated fitted values
    file_path = file.path(folder_w0_gender, paste0("fitted_", xvar, ".png")), # optional: path to save the plot
    res_ppi = 300, # resolution (pixels per inch)
    width_in = 7, height_in = 4 # width and height of the plot in inches
  )
  plot_glmer_fitted_continuous(
    model = full_bin_glmm_w0_age, # binomial model fitted with lme4::glmer()
    model_formula = model_formula_w0_age, # formula used to fit the model
    results_table = restab_w0_age, # results_table, output of get_results_table()
    data = data_w0_age, # data used to fit the model
    bootstrap_object = boot_full_bin_glmm_w0_age, # bootstrap object returned by boot.glmm.pred()
    link = "logit", # link function
    covariate_to_plot = xvar, # which covariate should be plotted
    fitted_resolution = 100, # resolution of calculated fitted values
    file_path = file.path(folder_w0_age, paste0("fitted_", xvar, ".png")), # optional: path to save the plot
    res_ppi = 300, # resolution (pixels per inch)
    width_in = 7, height_in = 4 # width and height of the plot in inches
  )
  plot_glmer_fitted_continuous(
    model = full_bin_glmm_w0_age_gender, # binomial model fitted with lme4::glmer()
    model_formula = model_formula_w0_age_gender, # formula used to fit the model
    results_table = restab_w0_age_gender, # results_table, output of get_results_table()
    data = data_w0_age_gender, # data used to fit the model
    bootstrap_object = boot_full_bin_glmm_w0_age_gender, # bootstrap object returned by boot.glmm.pred()
    link = "logit", # link function
    covariate_to_plot = xvar, # which covariate should be plotted
    fitted_resolution = 100, # resolution of calculated fitted values
    file_path = file.path(folder_w0_age_gender, paste0("fitted_", xvar, ".png")), # optional: path to save the plot
    res_ppi = 300, # resolution (pixels per inch)
    width_in = 7, height_in = 4 # width and height of the plot in inches
  )
}

# Plots showing age and gender --------------------------------------------

# age
plot_glmer_fitted_continuous(
  model = full_bin_glmm_all, # binomial model fitted with lme4::glmer()
  model_formula = model_formula_all, # formula used to fit the model
  results_table = restab_all, # results_table, output of get_results_table()
  data = data_all, # data used to fit the model
  bootstrap_object = boot_full_bin_glmm_all, # bootstrap object returned by boot.glmm.pred()
  link = "logit", # link function
  covariate_to_plot = "z.age", # which covariate should be plotted
  fitted_resolution = 100, # resolution of calculated fitted values
  file_path = file.path(folder_all_predictors, "fitted_z.age.png"), # optional: path to save the plot
  res_ppi = 300, # resolution (pixels per inch)
  width_in = 7, height_in = 4 # width and height of the plot in inches
)
plot_glmer_fitted_continuous(
  model = full_bin_glmm_w0_gender, # binomial model fitted with lme4::glmer()
  model_formula = model_formula_w0_gender, # formula used to fit the model
  results_table = restab_w0_gender, # results_table, output of get_results_table()
  data = data_w0_gender, # data used to fit the model
  bootstrap_object = boot_full_bin_glmm_w0_gender, # bootstrap object returned by boot.glmm.pred()
  link = "logit", # link function
  covariate_to_plot = "z.age", # which covariate should be plotted
  fitted_resolution = 100, # resolution of calculated fitted values
  file_path = file.path(folder_w0_gender, "fitted_z.age.png"), # optional: path to save the plot
  res_ppi = 300, # resolution (pixels per inch)
  width_in = 7, height_in = 4 # width and height of the plot in inches
)

# gender
plot_glmer_fitted_categorical(
  model = full_bin_glmm_all, # binomial model fitted with lme4::glmer()
  model_formula = model_formula_all, # formula used to fit the model
  results_table = restab_all, # results_table, output of get_results_table()
  data = data_all, # data used to fit the model
  bootstrap_object = boot_full_bin_glmm_all, # bootstrap object returned by boot.glmm.pred()
  link = "logit", # link function
  variable_to_plot = "gender", # which variable should be plotted
  file_path = file.path(folder_all_predictors, "fitted_gender.png"), # optional: path to save the plot
  res_ppi = 300, # resolution (pixels per inch)
  width_in = 7, height_in = 4 # width and height of the plot in inches
)
plot_glmer_fitted_categorical(
  model = full_bin_glmm_w0_age, # binomial model fitted with lme4::glmer()
  model_formula = model_formula_w0_age, # formula used to fit the model
  results_table = restab_w0_age, # results_table, output of get_results_table()
  data = data_w0_age, # data used to fit the model
  bootstrap_object = boot_full_bin_glmm_w0_age, # bootstrap object returned by boot.glmm.pred()
  link = "logit", # link function
  variable_to_plot = "gender", # which variable should be plotted
  file_path = file.path(folder_w0_age, "fitted_gender.png"), # optional: path to save the plot
  res_ppi = 300, # resolution (pixels per inch)
  width_in = 7, height_in = 4 # width and height of the plot in inches
)


# Compare estimates -------------------------------------------------------

# extract model terms names
xterms <- rownames(restab_all)
# create summary dataframe of model results
estimates <- rbind(
  data.frame( # model with all predictors
    term = xterms,
    model = "all_predictors",
    restab_all,
    row.names = NULL
  ),
  data.frame( # model without age
    term = xterms,
    model = "wo_age",
    restab_w0_age[xterms, ],
    row.names = NULL
  ),
  data.frame( # model without gender
    term = xterms,
    model = "wo_gender",
    restab_w0_gender[xterms, ],
    row.names = NULL
  ),
  data.frame( # model without age and gender
    term = xterms,
    model = "wo_age_and_gender",
    restab_w0_age_gender[xterms, ],
    row.names = NULL
  )
)
# define range of estimate values
xrange <- estimates[, c("Estimate", "lower_CI", "upper_CI")] |>
  as.matrix() |>
  range(na.rm = TRUE)
# prepare data for plotting
estimates$model <- factor(estimates$model, levels = c(
  "all_predictors", "wo_age", "wo_gender", "wo_age_and_gender"
))
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
