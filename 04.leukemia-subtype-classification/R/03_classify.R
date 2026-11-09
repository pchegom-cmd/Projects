# Multiclass subtype classifier evaluated with repeated, nested cross-validation.
# Gene selection and scaling are fitted inside each training fold to avoid leakage.
suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
})

set.seed(42)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

d <- readRDS("data/all_b_lineage.rds")
X <- t(d$expr)
y <- factor(d$meta$subtype)

stratified_folds <- function(y, k) {
  folds <- integer(length(y))
  for (level in levels(droplevels(y))) {
    idx <- which(y == level)
    folds[idx] <- sample(rep_len(sample(k), length(idx)))
  }
  folds
}

fit_model <- function(X, y, n_genes = 1000) {
  genes <- order(apply(X, 2, mad), decreasing = TRUE)[seq_len(n_genes)]
  centre <- colMeans(X[, genes])
  spread <- apply(X[, genes], 2, sd)
  z <- scale(X[, genes], centre, spread)
  fit <- cv.glmnet(z, y, family = "multinomial", type.multinomial = "grouped",
                   alpha = 0.5, foldid = stratified_folds(y, 3))
  list(fit = fit, genes = genes, centre = centre, spread = spread)
}

predict_model <- function(model, X) {
  z <- scale(X[, model$genes, drop = FALSE], model$centre, model$spread)
  as.character(predict(model$fit, z, s = "lambda.1se", type = "class"))
}

balanced_accuracy <- function(truth, pred) {
  mean(sapply(levels(truth), function(l) mean(pred[truth == l] == l)))
}

run_cv <- function(X, y, folds = 5) {
  id <- stratified_folds(y, folds)
  pred <- character(length(y))
  for (f in seq_len(folds)) {
    test <- id == f
    pred[test] <- predict_model(fit_model(X[!test, ], y[!test]), X[test, , drop = FALSE])
  }
  pred
}

repeats <- 5
predictions <- lapply(seq_len(repeats), function(r) run_cv(X, y))
scores <- sapply(predictions, function(p) balanced_accuracy(y, p))

# Negative control: with shuffled labels, performance should fall to chance level.
shuffled <- sample(y)
control <- balanced_accuracy(shuffled, run_cv(X, shuffled))

summary_table <- data.frame(
  metric = c("balanced_accuracy_mean", "balanced_accuracy_sd", "shuffled_label_control", "chance_level"),
  value  = round(c(mean(scores), sd(scores), control, 1 / nlevels(y)), 3)
)
print(summary_table)
write.csv(summary_table, "results/tables/classifier_performance.csv", row.names = FALSE)

confusion <- table(truth = rep(y, repeats), predicted = factor(unlist(predictions), levels = levels(y)))
print(confusion)
write.csv(confusion, "results/tables/confusion_matrix.csv")

df <- as.data.frame(confusion)
p <- ggplot(df, aes(predicted, truth, fill = Freq)) + geom_tile() +
  geom_text(aes(label = Freq)) + scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = paste("Confusion matrix over", repeats, "repeats of 5-fold CV")) +
  theme_minimal() + theme(legend.position = "none")
ggsave("results/figures/04_confusion_matrix.png", p, width = 6, height = 4.5, dpi = 150)

# Final model trained on all samples: for inspecting genes, not for estimating performance.
final <- fit_model(X, y)
coefs <- coef(final$fit, s = "lambda.1se")[[1]]
selected <- setdiff(rownames(coefs)[as.vector(coefs != 0)], "(Intercept)")
write.csv(data.frame(probe_set = selected), "results/tables/selected_probe_sets.csv", row.names = FALSE)
cat("Final model keeps", length(selected), "probe sets\n")
