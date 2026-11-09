# Molecular subtype classification of B-lineage acute lymphoblastic leukaemia

A small, reproducible learning project in R: unsupervised subtype discovery and supervised classification from public gene-expression data.

## Data

The Bioconductor [`ALL`](https://bioconductor.org/packages/ALL/) data package: pre-processed Affymetrix HG-U95Av2 microarray profiles from adult ALL patients (Chiaretti et al., *Blood*, 2004). The analysis uses B-lineage samples from molecular groups with at least 5 samples.

## Workflow

1. **Quality control** of per-sample expression distributions.
2. **Feature selection** of the most variable probe sets by median absolute deviation (2,000 for clustering, 1,000 for classification).
3. **Unsupervised analysis**: PCA, UMAP and k-means with consensus-based stability (PAC), silhouette, and adjusted Rand index against known subtypes.
4. **Classification**: elastic-net multinomial regression, assessed with repeated stratified 5-fold cross-validation. Gene selection and scaling are fitted only on each training fold to avoid data leakage. A shuffled-label run acts as a negative control.

## Reproduce

```bash
Rscript install_packages.R
Rscript run_all.R
```

Outputs are written to `results/`. Software versions are recorded in `results/sessionInfo.txt`.

## Results

- Balanced accuracy: **0.79 ± 0.03** (mean ± SD over 5 repeats of 5-fold cross-validation)
- Shuffled-label control: **0.24**, close to the chance level of 0.25
- Most stable number of clusters: **k = 3** (PAC 0.05), but with low silhouette (0.11) and almost no agreement with known subtypes (ARI 0.02)
- Recall by subtype: ALL1/AF4 98%, BCR/ABL 89%, NEG 85%, E2A/PBX1 44%. E2A/PBX1 has only 5 patients and was mostly misclassified as NEG.
- Main observation: subtype labels can be predicted well with supervised learning, but clustering on the most variable genes finds a stable structure that does not correspond to these subtypes. High-variance genes may therefore reflect other sources of variation, so stable clusters should not be assumed to be biological subtypes.

![UMAP](results/figures/03_umap.png)
![Confusion matrix](results/figures/04_confusion_matrix.png)

## Limitations

This is a small, single-cohort microarray dataset with very small groups, so performance estimates are uncertain and have not been validated in an independent cohort. The "NEG" group is defined by the absence of tested abnormalities and is biologically heterogeneous. This is a learning project and is not intended for clinical use.
