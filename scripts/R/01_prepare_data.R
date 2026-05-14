#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)][1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

config <- read_config()

paths <- config$paths
cohort_cfg <- config$cohort
feature_cfg <- config$features

ensure_dir(paths$interim_dir)
ensure_dir(paths$processed_dir)
ensure_dir(paths$tables_dir)
ensure_dir(file.path(paths$results_dir, "logs"))

message("[01] Reading series matrix metadata")
metadata_raw <- read_series_metadata(paths$series_matrix)
metadata_clean <- standardize_metadata(metadata_raw, feature_cfg$metadata_keep_fields)

message("[01] Reading expression matrix")
expr_probe <- read_expression_matrix(paths$series_matrix)

message("[01] Annotating GPL96 probes")
probe_annotation <- annotate_probes_hgu133a(rownames(expr_probe))
collapsed <- collapse_probes_to_genes(expr_probe, probe_annotation)
expr_gene <- collapsed$expression_gene
probe_selection <- collapsed$probe_selection

primary_meta <- subset_primary_cohort(
  metadata_clean,
  filter_field = cohort_cfg$her2_filter_field,
  filter_value = cohort_cfg$her2_negative_value,
  endpoint_field = cohort_cfg$primary_endpoint_field,
  positive = cohort_cfg$primary_positive_label,
  negative = cohort_cfg$primary_negative_label
)

secondary_meta <- subset_secondary_cohort(
  metadata_clean,
  filter_field = cohort_cfg$her2_filter_field,
  filter_value = cohort_cfg$her2_negative_value,
  endpoint_field = cohort_cfg$secondary_endpoint_field,
  positive_labels = cohort_cfg$secondary_positive_labels,
  negative_labels = cohort_cfg$secondary_negative_labels
)

primary_ids <- primary_meta$geo_accession
secondary_ids <- secondary_meta$geo_accession

expr_probe_primary <- expr_probe[, primary_ids, drop = FALSE]
expr_gene_primary <- expr_gene[, primary_ids, drop = FALSE]
expr_gene_secondary <- expr_gene[, secondary_ids, drop = FALSE]

cohort_summary <- bind_rows(
  metadata_clean %>% summarize(
    cohort = "all_samples",
    n = n(),
    pcr = sum(pathologic_response_pcr_rd == "pCR", na.rm = TRUE),
    rd = sum(pathologic_response_pcr_rd == "RD", na.rm = TRUE)
  ),
  primary_meta %>% summarize(
    cohort = "her2neg_primary",
    n = n(),
    pcr = sum(pathologic_response_pcr_rd == "pCR", na.rm = TRUE),
    rd = sum(pathologic_response_pcr_rd == "RD", na.rm = TRUE)
  ),
  secondary_meta %>% summarize(
    cohort = "her2neg_secondary",
    n = n(),
    rcb_0_i = sum(response_rcb_binary == "RCB-0/I", na.rm = TRUE),
    rcb_ii_iii = sum(response_rcb_binary == "RCB-II/III", na.rm = TRUE)
  )
)

probe_mapping_summary <- tibble::tibble(
  probe_count = nrow(expr_probe),
  annotated_probe_count = n_distinct(probe_annotation$PROBEID),
  unique_gene_count = nrow(expr_gene),
  collapsed_probe_count = nrow(probe_selection),
  dropped_probe_count = nrow(expr_probe) - n_distinct(probe_annotation$PROBEID)
)

write_table_safe(metadata_clean, file.path(paths$processed_dir, "metadata_clean.tsv"))
write_table_safe(primary_meta, file.path(paths$processed_dir, "metadata_primary.tsv"))
write_table_safe(secondary_meta, file.path(paths$processed_dir, "metadata_secondary.tsv"))
write_table_safe(probe_annotation, file.path(paths$processed_dir, "probe_annotation.tsv"))
write_table_safe(probe_selection, file.path(paths$processed_dir, "probe_selection.tsv"))
write_table_safe(cohort_summary, file.path(paths$tables_dir, "cohort_summary.tsv"))
write_table_safe(probe_mapping_summary, file.path(paths$tables_dir, "probe_mapping_summary.tsv"))

save_rds_safe(expr_probe, file.path(paths$interim_dir, "expression_probe_all.rds"))
save_rds_safe(expr_probe_primary, file.path(paths$interim_dir, "expression_probe_primary.rds"))
save_rds_safe(expr_gene, file.path(paths$processed_dir, "expression_gene_all.rds"))
save_rds_safe(expr_gene_primary, file.path(paths$processed_dir, "expression_gene_primary.rds"))
save_rds_safe(expr_gene_secondary, file.path(paths$processed_dir, "expression_gene_secondary.rds"))

message("[01] Completed preprocessing and cohort extraction")
