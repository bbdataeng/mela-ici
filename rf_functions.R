# Stratified split --------------------------------------------------------
# function for split stratified over the levels of a variable y
stratified_split <- function(y, p_train = 0.8) {
  stopifnot(is.factor(y))
  idx_train <- integer(0)
  for (lvl in levels(y)) {
    idx <- which(y == lvl)
    if (length(idx) == 0) next
    n_tr <- max(1, floor(length(idx) * p_train))
    idx_train <- c(idx_train, sample(idx, n_tr))
  }
  sort(unique(idx_train))
}

# Confusion matrix --------------------------------------------------------
get_confusion_matrix <- function(rf_object, testdata, show_sum = TRUE) {
  if (class(rf_object) == "ranger") {
    # Predict probabilities
    pred <- predict(rf_object, data = testdata, type = "response")
    # For probability = TRUE, this is a matrix of probs, with olumns = levels of response
    prob_mat <- pred$predictions
    # Turn probabilities into hard class calls
    pred_class <- colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]
    # Make it a factor with same levels as true response
    pred_class <- factor(pred_class, levels = levels(testdata$response))
    # Create confusion matrix
    cm_test <- table(Observed = testdata$response, Predicted = pred_class)
    # return object
    if (show_sum) {
      return(addmargins(cm_test))
    } else {
      return(cm_test)
    }
  } else if (class(rf_object) == "ordfor") {
    # predict test data
    pred <- predict(rf_object, newdata = testdata)
    # class predictions
    pred_class <- pred$ypred
    # Make it a factor with same levels as true response
    pred_class <- factor(pred_class, levels = levels(testdata$response), ordered = TRUE)
    # Create confusion matrix
    cm_test <- table(Observed = testdata$response, Predicted = pred_class)
    # return object
    if (show_sum) {
      return(addmargins(cm_test))
    } else {
      return(cm_test)
    }
  } else {
    stop("Unvalid class of 'rf_object'")
  }
}


# Plot ROC ----------------------------------------------------------------
# only works for binary RF created with ranger
plot_roc_auc <- function(rf_object, testdata, positive_level = "R", ...) {
  # predicted probabilities
  pred <- predict(
    rf_object,
    data = testdata, type = "response", probability = TRUE
  )
  # probability of positive level
  prob_pos <- pred$predictions[, positive_level]
  # get ROC
  xroc <- roc(
    response = testdata$response, predictor = prob_pos, quiet = TRUE
  )
  # plot ROC
  plot(xroc, print.auc = TRUE, auc.polygon = TRUE, ...)
}


# Accuracy metrics --------------------------------------------------------
get_accuracy_metrics <- function(rf_object, testdata, confusion_matrix, positive_level = "R") {
  # remove the "Sum" row and column, if present
  confusion_matrix <- confusion_matrix[
    rownames(confusion_matrix) != "Sum",
    colnames(confusion_matrix) != "Sum",
    drop = FALSE
  ]
  if (class(rf_object) == "ranger") { ## binary random forest ##
    ### calculate AUC ##
    pred <- predict( # predicted probabilities
      rf_object,
      data = testdata, type = "response", probability = TRUE
    )
    prob_pos <- pred$predictions[, positive_level] # probability of positive level
    # get AUC
    xauc <- roc(
      response = testdata$response, predictor = prob_pos, quiet = TRUE
    ) |>
      auc() |>
      as.numeric()
    ### calculate prediction performance metrics from Confusion Matrix
    # extract cells
    negative_level <- colnames(confusion_matrix) |> setdiff(positive_level)
    stopifnot(length(negative_level) == 1) # make sure that it is a 2x2 design
    TN <- confusion_matrix[negative_level, negative_level] # true negative
    FP <- confusion_matrix[negative_level, positive_level] # false positive
    FN <- confusion_matrix[positive_level, negative_level] # false negative
    TP <- confusion_matrix[positive_level, positive_level] # true positive
    N <- sum(confusion_matrix) # total observations
    # accuracy
    accuracy <- (TP + TN) / N
    # sensitivity (true positive rate)
    sensitivity <- TP / (TP + FN)
    # specificity (true negative rate)
    specificity <- TN / (TN + FP)
    # precision (positive predictive value)
    precision <- TP / (TP + FP)
    # F1 Score
    f1_score <- 2 * TP / (2 * TP + FP + FN)
    # Balanced Accuracy
    balanced_accuracy <- (sensitivity + specificity) / 2
    # Matthews Correlation Coefficient (MCC)
    den_mcc <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    if (den_mcc == 0) {
      mcc <- NA_real_
    } else {
      mcc <- (TP * TN - FP * FN) / den_mcc
    }
    # put everything in a named vector
    metrics <- c(
      AUC = xauc,
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      f1_score = f1_score,
      balanced_accuracy = balanced_accuracy,
      mcc = mcc
    )
  } else if (class(rf_object) == "ordfor") { ## ordinal random forest ##
    K <- nrow(confusion_matrix) # number of ordinal classes
    N <- sum(confusion_matrix) # total observations

    # overall accuracy
    overall_accuracy <- sum(diag(confusion_matrix)) / N

    # per-class sensitivity (recall aka true positive rate)
    # Sensitivity_k = TP_k / (true class k total) = diag / rowSums
    true_totals <- rowSums(confusion_matrix)
    pred_totals <- colSums(confusion_matrix)
    TP <- diag(confusion_matrix) # true positives
    sensitivity_per_class <- TP / true_totals


    # per-class precision
    # Precision_k = TP_k / (predicted class k total) = diag / colSums
    precision_per_class <- TP / pred_totals

    # per-class F1, macro-F1, weighted-F1
    ## F1_k = 2 * precision_k * sensitivity_k / (precision_k + sensitivity_k)
    f1_per_class <- 2 * precision_per_class * sensitivity_per_class /
      (precision_per_class + sensitivity_per_class)
    ## avoid NaN if precision + sensitivity = 0 for some class
    f1_per_class[is.nan(f1_per_class)] <- NA_real_
    ## macro-F1: unweighted average over classes
    macro_F1 <- mean(f1_per_class, na.rm = TRUE)
    ## weighted-F1: weighted by support (true class counts = rowSums)
    weighted_F1 <- weighted.mean(f1_per_class, w = true_totals, na.rm = TRUE)

    # Mean Absolute Error (MAE) in ordinal distance
    # treat classes as 1..K in the order of rows/cols
    class_indices <- seq_len(K)
    # matrix of absolute distances |i - j|
    dist_mat <- abs(outer(class_indices, class_indices, FUN = "-"))
    ## MAE = sum( count_ij * |i-j| ) / N
    MAE <- sum(confusion_matrix * dist_mat) / N

    # adjacent / non-adjacent error rates
    ## adjacent error: |i-j| = 1 (but i != j)
    ## non-adjacent error: |i-j| >= 2
    # we compute rates among all misclassified cases
    # logical matrices for adjacency
    adjacent_mask <- dist_mat == 1
    nonadjacent_mask <- dist_mat >= 2
    total_errors <- N - sum(diag(confusion_matrix))
    adjacent_errors <- sum(confusion_matrix * adjacent_mask)
    nonadjacent_errors <- sum(confusion_matrix * nonadjacent_mask)

    if (total_errors == 0) {
      adjacent_error_rate <- NA_real_
      nonadjacent_error_rate <- NA_real_
    } else {
      adjacent_error_rate <- adjacent_errors / total_errors
      nonadjacent_error_rate <- nonadjacent_errors / total_errors
    }

    # Quadratic Weighted Kappa (QWK)
    ## standard definition using:
    ##  - observed confusion matrix
    ##  - expected counts under independence
    ##  - quadratic weights w_ij = (i-j)^2 / (K-1)^2
    row_marginals <- true_totals
    col_marginals <- pred_totals
    # expected counts under independence
    expected <- outer(row_marginals, col_marginals, FUN = "*") / N
    ## quadratic weight matrix
    W <- outer(class_indices, class_indices,
      FUN = function(i, j) ((i - j)^2) / ((K - 1)^2)
    )
    ## observed and expected weighted sums
    obs_weighted <- sum(W * confusion_matrix)
    exp_weighted <- sum(W * expected)
    if (exp_weighted == 0) {
      qwk <- NA_real_
    } else {
      qwk <- 1 - obs_weighted / exp_weighted
    }

    # collect metrics in a convenient structure
    metrics <- list(
      overall_accuracy = overall_accuracy,
      sensitivity_per_class = sensitivity_per_class,
      precision_per_class = precision_per_class,
      f1_per_class = f1_per_class,
      macro_F1 = macro_F1,
      weighted_F1 = weighted_F1,
      MAE = MAE,
      adjacent_error_rate = adjacent_error_rate,
      nonadjacent_error_rate = nonadjacent_error_rate,
      QWK = qwk
    )
  } else {
    stop("Unvalid class of 'rf_object'")
  }
  return(metrics)
}


# Variable importance -----------------------------------------------------
get_importance <- function(rf_object) {
  if (class(rf_object) == "ranger") { ## binary random forest ##
    vimp <- rf_object$variable.importance |> sort(decreasing = TRUE)
  } else if (class(rf_object) == "ordfor") { ## ordinal random forest ##
    vimp <- rf_object$varimp |> sort(decreasing = TRUE)
  } else {
    stop("Unvalid class of 'rf_object'")
  }
  return(vimp)
}
