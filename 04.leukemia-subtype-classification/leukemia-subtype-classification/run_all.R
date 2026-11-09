# Run from the repository root: Rscript run_all.R
source("R/01_prepare_data.R", echo = FALSE)
source("R/02_explore_and_cluster.R", echo = FALSE)
source("R/03_classify.R", echo = FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
