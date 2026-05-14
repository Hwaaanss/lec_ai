#!/usr/bin/env Rscript

cran_packages <- c(
  "data.table",
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "ggplot2",
  "pheatmap",
  "yaml",
  "jsonlite",
  "e1071",
  "glmnet",
  "xgboost",
  "pROC",
  "PRROC",
  "umap"
)

bioc_packages <- c(
  "BiocManager",
  "limma",
  "GSVA",
  "msigdbr",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "annotate",
  "hgu133a.db"
)

ensure_installed <- function(pkgs, installer) {
  installed <- rownames(installed.packages())
  missing <- setdiff(pkgs, installed)
  if (!length(missing)) {
    message("All packages already installed for installer: ", deparse(substitute(installer)))
    return(invisible(NULL))
  }
  installer(missing)
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

ensure_installed(cran_packages, function(pkgs) {
  install.packages(pkgs)
})

if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager")
}

ensure_installed(bioc_packages, function(pkgs) {
  BiocManager::install(pkgs, ask = FALSE, update = FALSE)
})

message("R package installation check complete.")
