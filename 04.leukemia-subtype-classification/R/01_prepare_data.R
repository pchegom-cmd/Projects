# Prepare the B-lineage subset of the public Bioconductor `ALL` dataset.
suppressPackageStartupMessages({
  library(Biobase)
  library(ALL)
})

data(ALL)
dir.create("data", showWarnings = FALSE)

# Keep B-lineage samples only: T-ALL would dominate the main axis of variation.
eset <- ALL[, grepl("^B", as.character(ALL$BT))]

# Keep subtypes with at least 5 samples: smaller groups cannot be cross-validated.
subtype <- as.character(eset$mol.biol)
counts <- table(subtype)
print(counts)
eset <- eset[, subtype %in% names(counts[counts >= 5])]

meta <- data.frame(
  sample  = sampleNames(eset),
  subtype = as.character(eset$mol.biol),
  stage   = as.character(eset$BT),
  age     = eset$age,
  sex     = as.character(eset$sex),
  stringsAsFactors = FALSE
)

saveRDS(list(expr = exprs(eset), meta = meta), "data/all_b_lineage.rds")
cat("Saved", ncol(eset), "samples and", nrow(eset), "probe sets\n")
