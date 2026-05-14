#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)][1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

config <- read_config()
paths <- config$paths

ensure_dir(paths$summaries_dir)

read_if_exists <- function(path) {
  if (file.exists(path)) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    NULL
  }
}

cohort_tbl <- read_if_exists(file.path(paths$tables_dir, "cohort_summary.tsv"))
probe_tbl <- read_if_exists(file.path(paths$tables_dir, "probe_mapping_summary.tsv"))
assoc_tbl <- read_if_exists(file.path(paths$tables_dir, "eda_clinical_association_tests.tsv"))
go_tbl <- read_if_exists(file.path(paths$tables_dir, "go_bp_ora_up.tsv"))
kegg_tbl <- read_if_exists(file.path(paths$tables_dir, "kegg_ora_up.tsv"))
model_tbl <- read_if_exists(file.path(paths$tables_dir, "model_performance_nested_cv.tsv"))
transfer_tbl <- read_if_exists(file.path(paths$tables_dir, "model_performance_source_transfer.tsv"))
signature_tbl <- read_if_exists(file.path(paths$tables_dir, "top_predictive_signatures.tsv"))
pasnet_tuning_tbl <- read_if_exists(file.path(paths$tables_dir, "pasnet_tuning_summary.tsv"))

lines <- c(
  "# GSE25066 Analysis Summary",
  "",
  "## Cohort",
  ""
)

if (!is.null(cohort_tbl)) {
  lines <- c(lines, capture.output(print(cohort_tbl)), "")
}

if (!is.null(probe_tbl)) {
  lines <- c(lines, "## Probe-to-Gene Mapping", "", capture.output(print(probe_tbl)), "")
}

if (!is.null(assoc_tbl)) {
  lines <- c(lines, "## EDA Highlights", "", capture.output(print(head(arrange(assoc_tbl, p_value), 10))), "")
}

if (!is.null(go_tbl) || !is.null(kegg_tbl)) {
  lines <- c(lines, "## Biological Interpretation", "")
  if (!is.null(go_tbl)) {
    lines <- c(lines, "Top GO BP terms up in pCR:", "", capture.output(print(head(go_tbl, 10))), "")
  }
  if (!is.null(kegg_tbl)) {
    lines <- c(lines, "Top KEGG terms up in pCR:", "", capture.output(print(head(kegg_tbl, 10))), "")
  }
}

if (!is.null(model_tbl)) {
  lines <- c(lines, "## Modeling", "", "Nested CV performance:", "", capture.output(print(model_tbl)), "")
}

if (!is.null(transfer_tbl)) {
  lines <- c(lines, "Source-transfer performance:", "", capture.output(print(transfer_tbl)), "")
}

if (!is.null(pasnet_tuning_tbl)) {
  lines <- c(lines, "PASNet tuning summary:", "", capture.output(print(pasnet_tuning_tbl)), "")
}

if (!is.null(signature_tbl)) {
  lines <- c(lines, "Top predictive signatures:", "", capture.output(print(signature_tbl)), "")
}

writeLines(lines, file.path(paths$summaries_dir, "analysis_summary.md"))

message("[04] Summary markdown written")
