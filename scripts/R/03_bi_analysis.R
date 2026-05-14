#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)][1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

config <- read_config()
paths <- config$paths
cohort_cfg <- config$cohort
feature_cfg <- config$features

ensure_dir(file.path(paths$plots_dir, "bi"))
ensure_dir(paths$models_dir)

metadata <- readr::read_tsv(file.path(paths$processed_dir, "metadata_primary.tsv"), show_col_types = FALSE)
expr_gene <- readRDS(file.path(paths$processed_dir, "expression_gene_primary.rds"))

metadata <- metadata %>%
  mutate(
    response_primary = factor(response_primary, levels = c("RD", "pCR")),
    response_binary = ifelse(response_primary == "pCR", 1L, 0L)
  )

hallmark_sets <- make_collection_gene_sets("H")
kegg_sets <- make_collection_gene_sets("C2", "CP:KEGG_LEGACY")
gobp_sets <- make_collection_gene_sets("C5", "GO:BP")
immune_sets <- load_custom_immune_signatures(paths$immune_signatures_tsv)
kegg_term2gene <- tibble::enframe(kegg_sets, name = "gs_name", value = "gene_symbol") %>%
  tidyr::unnest(gene_symbol)

message("[03] Computing ssGSEA scores")
hallmark_scores <- compute_ssgsea_scores(expr_gene, hallmark_sets, min_size = feature_cfg$min_genes_per_set)
kegg_scores <- compute_ssgsea_scores(expr_gene, kegg_sets, min_size = feature_cfg$min_genes_per_set)
gobp_scores <- compute_ssgsea_scores(expr_gene, gobp_sets, min_size = feature_cfg$min_genes_per_set)
immune_scores <- compute_ssgsea_scores(expr_gene, immune_sets, min_size = 2)

score_df <- bind_rows(
  as.data.frame(hallmark_scores) %>% tibble::rownames_to_column("feature") %>% mutate(collection = "Hallmark"),
  as.data.frame(kegg_scores) %>% tibble::rownames_to_column("feature") %>% mutate(collection = "KEGG"),
  as.data.frame(gobp_scores) %>% tibble::rownames_to_column("feature") %>% mutate(collection = "GO_BP"),
  as.data.frame(immune_scores) %>% tibble::rownames_to_column("feature") %>% mutate(collection = "ImmuneSignature")
) %>%
  tidyr::pivot_longer(-c(feature, collection), names_to = "geo_accession", values_to = "score")

write_table_safe(score_df, file.path(paths$processed_dir, "pathway_scores_long.tsv"))

immune_sig_tbl <- readr::read_tsv(paths$immune_signatures_tsv, show_col_types = FALSE)
immune_axis_map <- immune_sig_tbl %>%
  distinct(signature_name, axis)

immune_score_wide <- as.data.frame(immune_scores) %>%
  tibble::rownames_to_column("feature") %>%
  tidyr::pivot_longer(-feature, names_to = "geo_accession", values_to = "score") %>%
  left_join(immune_axis_map, by = c("feature" = "signature_name"))

balance_df <- immune_score_wide %>%
  group_by(geo_accession, axis) %>%
  summarize(axis_mean = mean(score, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = axis, values_from = axis_mean) %>%
  mutate(
    immune_balance_index = activation - suppression,
    immune_balance_ratio = activation / pmax(abs(suppression), 1e-6)
  )

write_table_safe(balance_df, file.path(paths$processed_dir, "immune_balance_features.tsv"))

message("[03] Differential expression with limma")
design <- model.matrix(~ response_primary, data = metadata)
fit <- limma::lmFit(expr_gene, design)
fit <- limma::eBayes(fit)
de_tbl <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
de_tbl <- de_tbl %>%
  tibble::rownames_to_column("SYMBOL") %>%
  left_join(readr::read_tsv(file.path(paths$processed_dir, "probe_selection.tsv"), show_col_types = FALSE) %>% dplyr::select(SYMBOL, ENTREZID), by = "SYMBOL")

write_table_safe(de_tbl, file.path(paths$tables_dir, "de_primary_limma.tsv"))

sig_up <- de_tbl %>%
  filter(adj.P.Val <= feature_cfg$de_adj_p_threshold, logFC >= feature_cfg$de_logfc_threshold) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique() %>%
  as.character()

sig_down <- de_tbl %>%
  filter(adj.P.Val <= feature_cfg$de_adj_p_threshold, logFC <= -feature_cfg$de_logfc_threshold) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique() %>%
  as.character()

universe_ids <- de_tbl %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique() %>%
  as.character()

sig_up_symbols <- de_tbl %>%
  filter(adj.P.Val <= feature_cfg$de_adj_p_threshold, logFC >= feature_cfg$de_logfc_threshold) %>%
  filter(!is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

sig_down_symbols <- de_tbl %>%
  filter(adj.P.Val <= feature_cfg$de_adj_p_threshold, logFC <= -feature_cfg$de_logfc_threshold) %>%
  filter(!is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

universe_symbols <- de_tbl %>%
  filter(!is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

go_up <- if (length(sig_up) > 0) enrichGO(gene = sig_up, universe = universe_ids, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
go_down <- if (length(sig_down) > 0) enrichGO(gene = sig_down, universe = universe_ids, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
kegg_up <- if (length(sig_up_symbols) > 0) enricher(gene = sig_up_symbols, universe = universe_symbols, TERM2GENE = kegg_term2gene)
kegg_down <- if (length(sig_down_symbols) > 0) enricher(gene = sig_down_symbols, universe = universe_symbols, TERM2GENE = kegg_term2gene)

if (!is.null(go_up)) write_table_safe(as.data.frame(go_up), file.path(paths$tables_dir, "go_bp_ora_up.tsv"))
if (!is.null(go_down)) write_table_safe(as.data.frame(go_down), file.path(paths$tables_dir, "go_bp_ora_down.tsv"))
if (!is.null(kegg_up)) write_table_safe(as.data.frame(kegg_up), file.path(paths$tables_dir, "kegg_ora_up.tsv"))
if (!is.null(kegg_down)) write_table_safe(as.data.frame(kegg_down), file.path(paths$tables_dir, "kegg_ora_down.tsv"))

safe_enrichplot(go_up, file.path(paths$plots_dir, "bi", "go_bp_dotplot_up.png"), "GO BP ORA: up in pCR")
safe_enrichplot(kegg_up, file.path(paths$plots_dir, "bi", "kegg_dotplot_up.png"), "KEGG ORA: up in pCR")

ranked_stats <- de_tbl$logFC
names(ranked_stats) <- de_tbl$ENTREZID
ranked_stats <- ranked_stats[!is.na(names(ranked_stats))]
ranked_stats <- tapply(ranked_stats, names(ranked_stats), max)
ranked_stats <- sort(unlist(ranked_stats), decreasing = TRUE)
ranked_stats <- make_ranked_stats_unique(ranked_stats)

ranked_stats_symbol <- de_tbl$logFC
names(ranked_stats_symbol) <- de_tbl$SYMBOL
ranked_stats_symbol <- ranked_stats_symbol[!is.na(names(ranked_stats_symbol))]
ranked_stats_symbol <- tapply(ranked_stats_symbol, names(ranked_stats_symbol), max)
ranked_stats_symbol <- sort(unlist(ranked_stats_symbol), decreasing = TRUE)
ranked_stats_symbol <- make_ranked_stats_unique(ranked_stats_symbol)

gsea_go <- tryCatch(
  suppressWarnings(gseGO(geneList = ranked_stats, OrgDb = org.Hs.eg.db, ont = "BP", verbose = FALSE, eps = 0)),
  error = function(e) NULL
)
gsea_kegg <- tryCatch(
  suppressWarnings(GSEA(geneList = ranked_stats_symbol, TERM2GENE = kegg_term2gene, verbose = FALSE, eps = 0)),
  error = function(e) NULL
)

if (!is.null(gsea_go)) write_table_safe(as.data.frame(gsea_go), file.path(paths$tables_dir, "gsea_go_bp.tsv"))
if (!is.null(gsea_kegg)) write_table_safe(as.data.frame(gsea_kegg), file.path(paths$tables_dir, "gsea_kegg.tsv"))

volcano_tbl <- de_tbl %>%
  mutate(
    neg_log10_adj_p = -log10(pmax(adj.P.Val, 1e-300)),
    significance = case_when(
      adj.P.Val <= feature_cfg$de_adj_p_threshold & logFC >= feature_cfg$de_logfc_threshold ~ "Up in pCR",
      adj.P.Val <= feature_cfg$de_adj_p_threshold & logFC <= -feature_cfg$de_logfc_threshold ~ "Up in RD",
      TRUE ~ "NS"
    )
  )

p_volcano <- ggplot(volcano_tbl, aes(x = logFC, y = neg_log10_adj_p, color = significance)) +
  geom_point(alpha = 0.7, size = 1.4) +
  scale_color_manual(values = c("Up in pCR" = "#D64D4D", "Up in RD" = "#3E7CB1", "NS" = "grey70")) +
  labs(title = "Volcano Plot", x = "logFC (pCR vs RD)", y = "-log10 adjusted p-value") +
  theme_report()
ggsave(file.path(paths$plots_dir, "bi", "volcano_plot.png"), p_volcano, width = 9, height = 7, dpi = 300)

top_de_genes <- volcano_tbl %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 40) %>%
  pull(SYMBOL)

heat_de <- scale_rows(expr_gene[top_de_genes, , drop = FALSE])
annotation_col <- metadata %>%
  dplyr::select(geo_accession, response_primary, source, pam50_class) %>%
  tibble::column_to_rownames("geo_accession")
png(file.path(paths$plots_dir, "bi", "top_de_heatmap.png"), width = 2200, height = 1600, res = 220)
pheatmap::pheatmap(
  heat_de,
  annotation_col = annotation_col[colnames(heat_de), , drop = FALSE],
  show_colnames = FALSE,
  main = "Top Differentially Expressed Genes"
)
close_plot_device()

if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
  top_term <- as.data.frame(gsea_go)$ID[[1]]
  png(file.path(paths$plots_dir, "bi", "gsea_go_curve.png"), width = 2000, height = 1400, res = 220)
  print(clusterProfiler::gseaplot2(gsea_go, geneSetID = top_term, title = top_term))
  close_plot_device()
}

balance_plot_df <- balance_df %>%
  left_join(metadata %>% dplyr::select(geo_accession, response_primary), by = "geo_accession")

p_balance <- ggplot(balance_plot_df, aes(x = response_primary, y = immune_balance_index, fill = response_primary)) +
  geom_boxplot() +
  geom_jitter(width = 0.12, alpha = 0.55) +
  labs(title = "Immune Activation-Suppression Balance", x = NULL, y = "Balance index") +
  theme_report()
ggsave(file.path(paths$plots_dir, "bi", "activation_suppression_balance_boxplot.png"), p_balance, width = 8, height = 6, dpi = 300)

pathway_top <- score_df %>%
  left_join(metadata %>% dplyr::select(geo_accession, response_primary), by = "geo_accession") %>%
  group_by(feature) %>%
  summarize(stat = abs(wilcox.test(score ~ response_primary)$statistic), .groups = "drop") %>%
  arrange(desc(stat)) %>%
  slice_head(n = 40) %>%
  pull(feature)

pathway_heat <- score_df %>%
  filter(feature %in% pathway_top) %>%
  dplyr::select(feature, geo_accession, score) %>%
  tidyr::pivot_wider(names_from = geo_accession, values_from = score) %>%
  as.data.frame()
rownames(pathway_heat) <- pathway_heat$feature
pathway_heat$feature <- NULL
pathway_heat <- as.matrix(pathway_heat)

pathway_heat_scaled <- scale_rows(pathway_heat)
heatmap_res <- 220
cell_size_pt <- 5.5
device_width_in <- ((ncol(pathway_heat_scaled) * cell_size_pt) / 72) + 5.5
device_height_in <- ((nrow(pathway_heat_scaled) * cell_size_pt) / 72) + 2.8

png(
  file.path(paths$plots_dir, "bi", "pathway_score_heatmap.png"),
  width = ceiling(device_width_in * heatmap_res),
  height = ceiling(device_height_in * heatmap_res),
  res = heatmap_res
)
pheatmap::pheatmap(
  pathway_heat_scaled,
  annotation_col = annotation_col[colnames(pathway_heat), , drop = FALSE],
  show_colnames = FALSE,
  cellwidth = cell_size_pt,
  cellheight = cell_size_pt,
  fontsize_row = 9,
  main = "Top Pathway Scores"
)
close_plot_device()

hallmark_wide <- as.data.frame(hallmark_scores) %>%
  tibble::rownames_to_column("feature") %>%
  tidyr::pivot_longer(-feature, names_to = "geo_accession", values_to = "score") %>%
  tidyr::pivot_wider(names_from = feature, values_from = score)

immune_wide <- as.data.frame(immune_scores) %>%
  tibble::rownames_to_column("feature") %>%
  tidyr::pivot_longer(-feature, names_to = "geo_accession", values_to = "score") %>%
  tidyr::pivot_wider(names_from = feature, values_from = score)

pathway_wide <- hallmark_wide %>%
  left_join(immune_wide, by = "geo_accession") %>%
  left_join(balance_df, by = "geo_accession")

write_table_safe(pathway_wide, file.path(paths$processed_dir, "pathway_scores_wide.tsv"))
write_table_safe(
  tibble::tibble(
    collection = c(rep("Hallmark", nrow(hallmark_scores)), rep("ImmuneSignature", nrow(immune_scores)), rep("Balance", 4)),
    feature = c(rownames(hallmark_scores), rownames(immune_scores), c("activation", "suppression", "immune_balance_index", "immune_balance_ratio"))
  ),
  file.path(paths$tables_dir, "model_bi_feature_manifest.tsv")
)

pasnet_expression <- as.data.frame(t(expr_gene)) %>%
  tibble::rownames_to_column("geo_accession")
write_table_safe(pasnet_expression, file.path(paths$processed_dir, "pasnet_expression_primary.tsv"))

pasnet_gene_sets <- bind_rows(
  tibble::enframe(hallmark_sets, name = "pathway", value = "gene_symbol") %>%
    tidyr::unnest(gene_symbol) %>%
    mutate(collection = "Hallmark"),
  tibble::enframe(immune_sets, name = "pathway", value = "gene_symbol") %>%
    tidyr::unnest(gene_symbol) %>%
    mutate(collection = "ImmuneSignature")
) %>%
  dplyr::select(collection, pathway, gene_symbol) %>%
  distinct()
write_table_safe(pasnet_gene_sets, file.path(paths$processed_dir, "pasnet_gene_sets.tsv"))

clinical_input <- metadata %>%
  dplyr::select(geo_accession, all_of(feature_cfg$clinical_fields), response_primary, response_binary, source)
clinical_input <- coerce_categorical_na(clinical_input, setdiff(names(clinical_input), c("geo_accession", "age_years", "response_binary")))

bi_input <- pathway_wide %>%
  left_join(metadata %>% dplyr::select(geo_accession, response_primary, response_binary, source), by = "geo_accession")

combined_input <- clinical_input %>%
  left_join(pathway_wide, by = "geo_accession")

write_table_safe(clinical_input, file.path(paths$processed_dir, "model_input_clinical_primary.tsv"))
write_table_safe(bi_input, file.path(paths$processed_dir, "model_input_bi_primary.tsv"))
write_table_safe(combined_input, file.path(paths$processed_dir, "model_input_combined_primary.tsv"))

message("[03] BI analysis and model input export completed")
