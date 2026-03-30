dir.create("output", showWarnings = FALSE)
dir.create("output/plots", showWarnings = FALSE)

options(stringsAsFactors = FALSE)

auc_rank <- function(y, p) {
  y <- as.integer(y)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calc_metrics <- function(y, p, threshold = 0.5) {
  pred <- ifelse(p >= threshold, 1L, 0L)
  tp <- sum(pred == 1 & y == 1)
  tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)

  accuracy <- (tp + tn) / length(y)
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  sensitivity <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  f1 <- if (
    is.na(precision) || is.na(sensitivity) || (precision + sensitivity) == 0
  ) {
    NA_real_
  } else {
    2 * precision * sensitivity / (precision + sensitivity)
  }

  list(
    metrics = data.frame(
      accuracy = accuracy,
      error_rate = 1 - accuracy,
      precision = precision,
      sensitivity = sensitivity,
      specificity = specificity,
      f1 = f1,
      auc = auc_rank(y, p)
    ),
    confusion = matrix(
      c(tn, fp, fn, tp),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        Actual = c("0", "1"),
        Predicted = c("0", "1")
      )
    ),
    pred = pred
  )
}

roc_df <- function(y, p) {
  thresholds <- sort(unique(c(0, p, 1)))
  out <- lapply(thresholds, function(thr) {
    pred <- ifelse(p >= thr, 1L, 0L)
    tp <- sum(pred == 1 & y == 1)
    tn <- sum(pred == 0 & y == 0)
    fp <- sum(pred == 1 & y == 0)
    fn <- sum(pred == 0 & y == 1)
    tpr <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    fpr <- if ((fp + tn) == 0) 0 else fp / (fp + tn)
    data.frame(threshold = thr, tpr = tpr, fpr = fpr)
  })
  unique(do.call(rbind, out))
}

train <- read.csv("data/crime-training-data_modified.csv")
eval_data <- read.csv("data/crime-evaluation-data_modified.csv")

train$target <- as.integer(train$target)

summary_train <- train

train$chas <- factor(train$chas)
eval_data$chas <- factor(eval_data$chas, levels = levels(train$chas))

train_transformed <- within(train, {
  nox_sq <- nox^2
  rm_sq <- rm^2
  dis_sq <- dis^2
  medv_sq <- medv^2
  log_lstat <- log(lstat)
  log_tax <- log(tax)
})

eval_transformed <- within(eval_data, {
  nox_sq <- nox^2
  rm_sq <- rm^2
  dis_sq <- dis^2
  medv_sq <- medv^2
  log_lstat <- log(lstat)
  log_tax <- log(tax)
})

summary_stats <- data.frame(
  variable = names(summary_train)[names(summary_train) != "target"],
  mean = sapply(summary_train[names(summary_train) != "target"], mean),
  sd = sapply(summary_train[names(summary_train) != "target"], sd),
  median = sapply(summary_train[names(summary_train) != "target"], median),
  min = sapply(summary_train[names(summary_train) != "target"], min),
  max = sapply(summary_train[names(summary_train) != "target"], max),
  row.names = NULL
)

numeric_train <- train
numeric_train$chas <- as.integer(as.character(numeric_train$chas))
correlations <- sort(cor(numeric_train)[, "target"], decreasing = TRUE)
correlation_df <- data.frame(
  variable = names(correlations),
  correlation_with_target = as.numeric(correlations),
  row.names = NULL
)

spearman_correlations <- sapply(
  numeric_train[names(numeric_train) != "target"],
  function(x) suppressWarnings(cor(x, numeric_train$target, method = "spearman"))
)
spearman_df <- data.frame(
  variable = names(spearman_correlations),
  spearman_with_target = as.numeric(spearman_correlations),
  row.names = NULL
)
correlation_compare <- merge(correlation_df, spearman_df, by = "variable", all = TRUE)

model_1 <- glm(target ~ ., data = train, family = binomial())
model_2 <- MASS::stepAIC(model_1, trace = FALSE)
model_3_full <- glm(
  target ~ zn + indus + chas + nox + nox_sq + rm + rm_sq + age + dis + dis_sq +
    rad + log_tax + ptratio + lstat + log_lstat + medv + medv_sq,
  data = train_transformed,
  family = binomial()
)
model_3 <- MASS::stepAIC(model_3_full, trace = FALSE)

models <- list(
  "Model 1: Full main-effects" = model_1,
  "Model 2: Stepwise main-effects" = model_2,
  "Model 3: Stepwise transformed" = model_3
)

model_metrics <- list()
model_coefficients <- list()
model_confusions <- list()
roc_data <- list()

for (name in names(models)) {
  model <- models[[name]]
  preds <- predict(model, type = "response")
  model_eval <- calc_metrics(train$target, preds)

  model_metrics[[name]] <- cbind(
    model = name,
    data.frame(
      aic = AIC(model),
      log_likelihood = as.numeric(logLik(model))
    ),
    model_eval$metrics
  )
  model_coefficients[[name]] <- cbind(
    model = name,
    term = rownames(summary(model)$coefficients),
    data.frame(summary(model)$coefficients, row.names = NULL)
  )
  conf <- as.data.frame(as.table(model_eval$confusion))
  names(conf) <- c("actual", "predicted", "count")
  conf$model <- name
  model_confusions[[name]] <- conf
  roc_piece <- roc_df(train$target, preds)
  roc_piece$model <- name
  roc_data[[name]] <- roc_piece
}

model_metrics_df <- do.call(rbind, model_metrics)
row.names(model_metrics_df) <- NULL
coefficients_df <- do.call(rbind, model_coefficients)
row.names(coefficients_df) <- NULL
confusions_df <- do.call(rbind, model_confusions)
row.names(confusions_df) <- NULL
roc_df_all <- do.call(rbind, roc_data)
row.names(roc_df_all) <- NULL

best_model <- model_3
best_model_name <- "Model 3: Stepwise transformed"
best_predictions <- predict(best_model, newdata = eval_transformed, type = "response")

evaluation_predictions <- data.frame(
  row_id = seq_len(nrow(eval_data)),
  probability = round(best_predictions, 6),
  classification = ifelse(best_predictions >= 0.5, 1L, 0L)
)

write.csv(summary_stats, "output/summary_statistics.csv", row.names = FALSE)
write.csv(correlation_df, "output/correlation_with_target.csv", row.names = FALSE)
write.csv(correlation_compare, "output/correlation_compare.csv", row.names = FALSE)
write.csv(model_metrics_df, "output/model_metrics.csv", row.names = FALSE)
write.csv(coefficients_df, "output/model_coefficients.csv", row.names = FALSE)
write.csv(confusions_df, "output/model_confusions.csv", row.names = FALSE)
write.csv(roc_df_all, "output/roc_curve_points.csv", row.names = FALSE)
write.csv(evaluation_predictions, "output/evaluation_predictions.csv", row.names = FALSE)

png("output/plots/target_distribution.png", width = 900, height = 600)
barplot(
  table(train$target),
  col = c("#7AA6C2", "#D95F02"),
  main = "Crime Risk Class Distribution",
  xlab = "Target Class",
  ylab = "Count"
)
dev.off()

png("output/plots/key_boxplots.png", width = 1200, height = 800)
par(mfrow = c(2, 2))
boxplot(nox ~ target, data = train, col = c("#7AA6C2", "#D95F02"),
        main = "NOX by Target", xlab = "Target", ylab = "NOX")
boxplot(dis ~ target, data = train, col = c("#7AA6C2", "#D95F02"),
        main = "DIS by Target", xlab = "Target", ylab = "DIS")
boxplot(rad ~ target, data = train, col = c("#7AA6C2", "#D95F02"),
        main = "RAD by Target", xlab = "Target", ylab = "RAD")
boxplot(medv ~ target, data = train, col = c("#7AA6C2", "#D95F02"),
        main = "MEDV by Target", xlab = "Target", ylab = "MEDV")
dev.off()

png("output/plots/roc_comparison.png", width = 900, height = 700)
plot(c(0, 1), c(0, 1), type = "l", lty = 2, col = "gray50",
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = "ROC Curves for Candidate Models")
palette <- c("#1B9E77", "#D95F02", "#7570B3")
for (i in seq_along(models)) {
  name <- names(models)[i]
  piece <- roc_df_all[roc_df_all$model == name, ]
  piece <- piece[order(piece$fpr, piece$tpr), ]
  lines(piece$fpr, piece$tpr, col = palette[i], lwd = 2)
}
legend(
  "bottomright",
  legend = sprintf(
    "%s (AUC = %.3f)",
    model_metrics_df$model,
    model_metrics_df$auc
  ),
  col = palette,
  lwd = 2,
  bty = "n"
)
dev.off()

saveRDS(
  list(
    train = train,
    eval_data = eval_data,
    summary_stats = summary_stats,
    correlation_df = correlation_df,
    correlation_compare = correlation_compare,
    model_metrics = model_metrics_df,
    coefficients = coefficients_df,
    confusions = confusions_df,
    evaluation_predictions = evaluation_predictions,
    best_model_name = best_model_name
  ),
  "output/report_objects.rds"
)

cat("Analysis complete. Outputs written to output/.\n")
