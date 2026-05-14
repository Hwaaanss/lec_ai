#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)][1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

config <- read_config()
paths <- config$paths
feature_cfg <- config$features

ensure_dir(file.path(paths$plots_dir, "eda"))

metadata <- readr::read_tsv(file.path(paths$processed_dir, "metadata_primary.tsv"), show_col_types = FALSE)
expr_probe <- readRDS(file.path(paths$interim_dir, "expression_probe_primary.rds"))
expr_gene <- readRDS(file.path(paths$processed_dir, "expression_gene_primary.rds"))

metadata <- metadata %>%
  mutate(
    response_primary = factor(response_primary, levels = c("RD", "pCR")),
    source = replace_na(source, "NA"),
    pam50_class = replace_na(pam50_class, "NA"),
    grade = replace_na(grade, "NA"),
    clinical_ajcc_stage = replace_na(clinical_ajcc_stage, "NA"),
    clinical_nodal_status = replace_na(clinical_nodal_status, "NA"),
    `er_status_ihc_esr1_for indeterminate` = replace_na(`er_status_ihc_esr1_for indeterminate`, "NA")
  )

clinical_fields <- c(
  "source",
  "response_primary",
  "er_status_ihc_esr1_for indeterminate",
  "grade",
  "clinical_ajcc_stage",
  "clinical_nodal_status",
  "pam50_class"
)

dist_tbl <- bind_rows(lapply(clinical_fields, function(field) {
  metadata %>%
    count(.data[[field]], response_primary, name = "n") %>%
    transmute(variable = field, level = .data[[field]], response_primary, n)
}))
write_table_safe(dist_tbl, file.path(paths$tables_dir, "eda_clinical_distributions.tsv"))

missing_tbl <- tibble::tibble(
  variable = names(metadata),
  missing_n = vapply(metadata, function(x) sum(is.na(x) | x == "NA"), numeric(1)),
  missing_pct = missing_n / nrow(metadata)
)
write_table_safe(missing_tbl, file.path(paths$tables_dir, "eda_missing_summary.tsv"))

p_class <- metadata %>%
  count(response_primary, name = "n") %>%
  ggplot(aes(x = response_primary, y = n, fill = response_primary)) +
  geom_col(width = 0.7) +
  labs(title = "Primary Response Class Distribution", x = NULL, y = "Samples") +
  theme_report()
ggsave(file.path(paths$plots_dir, "eda", "class_distribution_barplot.png"), p_class, width = 8, height = 6, dpi = 300)

p_missing <- missing_tbl %>%
  ggplot(aes(x = reorder(variable, missing_pct), y = missing_pct)) +
  geom_col(fill = "#4C78A8") +
  coord_flip() +
  labs(title = "Missingness Summary", x = NULL, y = "Missing fraction") +
  theme_report()
ggsave(file.path(paths$plots_dir, "eda", "missing_summary.png"), p_missing, width = 9, height = 7, dpi = 300)

density_probe_idx <- seq_len(min(1000, nrow(expr_probe)))
density_long <- as.data.frame(expr_probe[density_probe_idx, , drop = FALSE]) %>%
  tibble::rownames_to_column("probe_id") %>%
  tidyr::pivot_longer(-probe_id, names_to = "sample", values_to = "expression") %>%
  left_join(metadata %>% dplyr::select(geo_accession, response_primary), by = c("sample" = "geo_accession"))

p_density <- ggplot(density_long, aes(x = expression, group = sample, color = response_primary)) +
  geom_density(alpha = 0.15) +
  labs(title = "Probe-Level Expression Density", x = "Expression", y = "Density") +
  theme_report()
ggsave(file.path(paths$plots_dir, "eda", "expression_density.png"), p_density, width = 10, height = 7, dpi = 300)

box_df <- tibble::tibble(
  sample = colnames(expr_gene),
  median_expression = apply(expr_gene, 2, median, na.rm = TRUE)
) %>%
  left_join(metadata %>% dplyr::select(geo_accession, response_primary, source), by = c("sample" = "geo_accession"))

p_box <- ggplot(box_df, aes(x = response_primary, y = median_expression, fill = response_primary)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(title = "Sample Median Expression by Response", x = NULL, y = "Median expression") +
  theme_report()
ggsave(file.path(paths$plots_dir, "eda", "expression_boxplot.png"), p_box, width = 8, height = 6, dpi = 300)

pca <- prcomp(t(expr_gene), center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  tibble::rownames_to_column("geo_accession") %>%
  left_join(metadata, by = "geo_accession")

for (color_field in c("response_primary", "source", "er_status_ihc_esr1_for indeterminate", "pam50_class")) {
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_field]])) +
    geom_point(size = 2.3, alpha = 0.85) +
    labs(title = paste("PCA colored by", color_field), x = "PC1", y = "PC2", color = color_field) +
    theme_report()
  ggsave(file.path(paths$plots_dir, "eda", paste0("pca_", color_field, ".png")), p, width = 9, height = 7, dpi = 300)
}

mi_top <- extract_top_mi_genes(
  expr_gene = expr_gene[, metadata$geo_accession, drop = FALSE],
  response = metadata$response_primary,
  n = 500,
  n_bins = 5
)
mi_score_tbl <- tibble::tibble(
  gene_symbol = names(mi_top$scores),
  mutual_information = as.numeric(mi_top$scores)
)
write_table_safe(mi_score_tbl, file.path(paths$tables_dir, "eda_mutual_information_response_primary.tsv"))

mi_pca <- prcomp(t(expr_gene[mi_top$genes, metadata$geo_accession, drop = FALSE]), center = TRUE, scale. = TRUE)
mi_pca_df <- as.data.frame(mi_pca$x[, 1:2]) %>%
  tibble::rownames_to_column("geo_accession") %>%
  left_join(metadata, by = "geo_accession")

p_mi_pca <- ggplot(mi_pca_df, aes(x = PC1, y = PC2, color = response_primary)) +
  geom_point(size = 2.3, alpha = 0.85) +
  labs(
    title = "PCA colored by response_primary (Mutual Information-selected genes)",
    x = "PC1",
    y = "PC2",
    color = "response_primary"
  ) +
  theme_report()
ggsave(
  file.path(paths$plots_dir, "eda", "pca_response_primary_mutual_information.png"),
  p_mi_pca,
  width = 9,
  height = 7,
  dpi = 300
)

top_genes <- extract_top_variable_genes(expr_gene, feature_cfg$top_variable_genes_heatmap)
heat_mat <- scale_rows(expr_gene[top_genes, , drop = FALSE])
annotation_col <- metadata %>%
  dplyr::select(geo_accession, response_primary, source, pam50_class) %>%
  tibble::column_to_rownames("geo_accession")
png(file.path(paths$plots_dir, "eda", "top_variable_gene_heatmap.png"), width = 2200, height = 1600, res = 220)
pheatmap::pheatmap(
  heat_mat,
  annotation_col = annotation_col[colnames(heat_mat), , drop = FALSE],
  show_colnames = FALSE,
  fontsize_row = 7,
  main = "Top Variable Genes"
)
close_plot_device()

assoc_results <- bind_rows(
  tibble::tibble(
    variable = "age_years",
    test = "wilcox",
    p_value = wilcox.test(age_years ~ response_primary, data = metadata)$p.value
  ),
  lapply(
    c("source", "grade", "clinical_t_stage", "clinical_nodal_status", "clinical_ajcc_stage", "er_status_ihc_esr1_for indeterminate", "pam50_class"),
    function(field) {
      test_res <- categorical_association_test(metadata[[field]], metadata$response_primary)
      tibble::tibble(variable = field, test = test_res$test, p_value = test_res$p_value)
    }
  ) %>% bind_rows()
) %>%
  arrange(p_value)
write_table_safe(assoc_results, file.path(paths$tables_dir, "eda_clinical_association_tests.tsv"))

assoc_plot_df <- metadata %>%
  count(pam50_class, response_primary, name = "n")

p_assoc <- ggplot(assoc_plot_df, aes(x = pam50_class, y = n, fill = response_primary)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Clinical Response Association: PAM50", x = "PAM50 class", y = "Within-class response fraction") +
  theme_report() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(paths$plots_dir, "eda", "clinical_response_association.png"), p_assoc, width = 10, height = 7, dpi = 300)

message("[02] EDA outputs generated")
