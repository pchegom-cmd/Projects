# Quality control, UMAP and unsupervised clustering with stability assessment.
suppressPackageStartupMessages({
  library(ggplot2)
  library(uwot)
  library(cluster)
  library(mclust)
})

set.seed(42)
dir.create("results/figures_2000", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables_2000", recursive = TRUE, showWarnings = FALSE)

d <- readRDS("data/all_b_lineage.rds")
expr <- d$expr
meta <- d$meta

# --- Quality control ---
long <- data.frame(sample = rep(colnames(expr), each = nrow(expr)), value = as.vector(expr))
p <- ggplot(long, aes(sample, value)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Sample", y = "log2 expression", title = "Expression distribution per sample") +
  theme_minimal() + theme(axis.text.x = element_blank())
ggsave("results/figures/01_qc_boxplot.png", p, width = 9, height = 4, dpi = 150)

# --- Feature selection: the most variable genes, without using subtype labels ---
n_genes <- 2000
top <- order(apply(expr, 1, mad), decreasing = TRUE)[seq_len(n_genes)]
x <- t(scale(t(expr[top, ])))   # samples in columns; each gene z-scored
samples <- t(x)                 # samples in rows, for clustering functions

# --- PCA ---
pca <- prcomp(samples)
var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
meta$PC1 <- pca$x[, 1]; meta$PC2 <- pca$x[, 2]
p <- ggplot(meta, aes(PC1, PC2, colour = subtype)) + geom_point(size = 2) +
  labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"), title = "PCA") +
  theme_minimal()
ggsave("results/figures/02_pca.png", p, width = 6, height = 4.5, dpi = 150)

# --- Consensus clustering by subsampling ---
consensus_stats <- function(samples, k, B = 100, frac = 0.8) {
  n <- nrow(samples)
  together <- matrix(0, n, n)
  sampled <- matrix(0, n, n)
  for (b in seq_len(B)) {
    idx <- sample(n, floor(frac * n))
    cl <- kmeans(samples[idx, ], centers = k, nstart = 5, iter.max = 50)$cluster
    sampled[idx, idx] <- sampled[idx, idx] + 1
    same <- outer(cl, cl, "==") * 1
    together[idx, idx] <- together[idx, idx] + same
  }
  cons <- together / pmax(sampled, 1)
  values <- cons[upper.tri(cons)]
  # PAC: share of ambiguous pairs. Lower means more stable clusters.
  mean(values > 0.1 & values < 0.9)
}

distances <- dist(samples)
stats <- do.call(rbind, lapply(2:6, function(k) {
  full <- kmeans(samples, centers = k, nstart = 25, iter.max = 50)$cluster
  data.frame(
    k = k,
    PAC = consensus_stats(samples, k),
    silhouette = mean(silhouette(full, distances)[, "sil_width"]),
    ARI_vs_subtype = adjustedRandIndex(full, meta$subtype)
  )
}))
print(stats)
write.csv(stats, "results/tables/cluster_stability.csv", row.names = FALSE)

best_k <- stats$k[which.min(stats$PAC)]
cat("Most stable k (lowest PAC):", best_k, "\n")
meta$cluster <- factor(kmeans(samples, centers = best_k, nstart = 25, iter.max = 50)$cluster)
write.csv(table(cluster = meta$cluster, subtype = meta$subtype), "results/tables/cluster_vs_subtype.csv")

# --- UMAP, for visualisation only ---
um <- umap(samples, n_neighbors = 15, min_dist = 0.1, n_threads = 1)
meta$UMAP1 <- um[, 1]; meta$UMAP2 <- um[, 2]
p <- ggplot(meta, aes(UMAP1, UMAP2, colour = subtype, shape = cluster)) + geom_point(size = 2.5) +
  labs(title = "UMAP of B-lineage ALL samples") + theme_minimal()
ggsave("results/figures/03_umap.png", p, width = 6.5, height = 4.5, dpi = 150)
