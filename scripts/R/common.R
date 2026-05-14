required_packages <- c(
  "yaml",
  "data.table",
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "ggplot2",
  "pheatmap",
  "limma",
  "GSVA",
  "msigdbr",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "pROC"
)

load_or_stop <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing, collapse = ", "),
      ". Run scripts/R/install_r_packages.R first."
    )
  }
  invisible(TRUE)
}

load_or_stop(required_packages)

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(limma)
  library(GSVA)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(pROC)
})

read_config <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- if (length(args) > 0) args[[1]] else "config/analysis_config.yaml"
  yaml::read_yaml(config_path)
}

project_path <- function(...) {
  file.path(...)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

na_clean <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "null", "<NA>")] <- NA_character_
  x
}

safe_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out
}

read_series_metadata <- function(series_path) {
  con <- gzfile(series_path, open = "rt")
  on.exit(close(con), add = TRUE)

  sample_ids <- NULL
  sample_titles <- NULL
  char_lines <- list()
  other_sample_lines <- list()

  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0L) {
      break
    }
    if (!nzchar(trimws(line))) {
      next
    }
    if (startsWith(line, "!series_matrix_table_begin")) {
      break
    }
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) == 0L || !nzchar(fields[[1]])) {
      next
    }
    tag <- fields[[1]]
    values <- gsub('^"|"$', "", fields[-1])
    if (identical(tag, "!Sample_geo_accession")) {
      sample_ids <- values
    } else if (identical(tag, "!Sample_title")) {
      sample_titles <- values
    } else if (identical(tag, "!Sample_characteristics_ch1")) {
      char_lines[[length(char_lines) + 1L]] <- values
    } else if (startsWith(tag, "!Sample_")) {
      other_sample_lines[[tag]] <- values
    }
  }

  if (is.null(sample_ids)) {
    stop("Failed to parse sample accessions from series matrix.")
  }

  meta <- data.frame(
    geo_accession = sample_ids,
    title = sample_titles,
    stringsAsFactors = FALSE
  )

  if (length(other_sample_lines) > 0) {
    for (nm in names(other_sample_lines)) {
      clean_nm <- sub("^!Sample_", "", nm)
      meta[[clean_nm]] <- other_sample_lines[[nm]]
    }
  }

  if (length(char_lines) > 0) {
    for (values in char_lines) {
      pieces <- strsplit(values, ":", fixed = TRUE)
      key <- trimws(tolower(vapply(pieces, `[`, character(1), 1)))
      value <- trimws(vapply(pieces, function(x) {
        if (length(x) < 2) {
          ""
        } else {
          paste(x[-1], collapse = ":")
        }
      }, character(1)))
      key_name <- key[[1]]
      meta[[key_name]] <- value
    }
  }

  meta[] <- lapply(meta, na_clean)
  meta
}

read_expression_matrix <- function(series_path) {
  expr <- data.table::fread(
    cmd = sprintf(
      "gzip -dc %s | awk 'BEGIN{flag=0} /^!series_matrix_table_begin$/{flag=1; next} /^!series_matrix_table_end$/{flag=0} flag'",
      shQuote(series_path)
    ),
    header = TRUE,
    data.table = FALSE
  )
  rownames(expr) <- expr[[1]]
  expr[[1]] <- NULL
  as.matrix(expr)
}

require_annotation_package <- function() {
  if (!requireNamespace("hgu133a.db", quietly = TRUE)) {
    stop(
      "Package 'hgu133a.db' is required for GPL96 probe annotation. ",
      "Install it with scripts/R/install_r_packages.R."
    )
  }
}

annotate_probes_hgu133a <- function(probe_ids) {
  require_annotation_package()
  anno_db <- get("hgu133a.db", envir = asNamespace("hgu133a.db"))
  probe_ids <- unique(probe_ids)
  quiet_map_ids <- function(column) {
    suppressMessages(
      AnnotationDbi::mapIds(
        anno_db,
        keys = probe_ids,
        keytype = "PROBEID",
        column = column,
        multiVals = "first"
      )
    )
  }
  ann <- data.frame(
    PROBEID = probe_ids,
    SYMBOL = quiet_map_ids("SYMBOL"),
    GENENAME = quiet_map_ids("GENENAME"),
    ENTREZID = quiet_map_ids("ENTREZID"),
    stringsAsFactors = FALSE
  )
  ann$PROBEID <- as.character(ann$PROBEID)
  ann$SYMBOL <- as.character(ann$SYMBOL)
  ann$GENENAME <- as.character(ann$GENENAME)
  ann$ENTREZID <- as.character(ann$ENTREZID)
  ann <- ann %>%
    mutate(
      SYMBOL = na_clean(SYMBOL),
      GENENAME = na_clean(GENENAME),
      ENTREZID = na_clean(ENTREZID)
    ) %>%
    filter(!is.na(SYMBOL))
  ann
}

collapse_probes_to_genes <- function(expr_probe, annotation_tbl) {
  annotation_tbl <- annotation_tbl %>%
    filter(PROBEID %in% rownames(expr_probe)) %>%
    distinct(PROBEID, SYMBOL, .keep_all = TRUE)

  if (nrow(annotation_tbl) == 0) {
    stop("No probe annotations matched the expression matrix.")
  }

  probe_medians <- apply(expr_probe[annotation_tbl$PROBEID, , drop = FALSE], 1, median, na.rm = TRUE)
  annotation_tbl$probe_median <- probe_medians[annotation_tbl$PROBEID]

  best <- annotation_tbl %>%
    arrange(desc(probe_median)) %>%
    group_by(SYMBOL) %>%
    dplyr::slice(1) %>%
    ungroup()

  gene_mat <- expr_probe[best$PROBEID, , drop = FALSE]
  rownames(gene_mat) <- best$SYMBOL

  list(
    expression_gene = gene_mat,
    probe_selection = best
  )
}

standardize_metadata <- function(meta, keep_fields) {
  if (!"sample id" %in% names(meta)) {
    stop("Expected 'sample id' in parsed metadata.")
  }

  meta$sample_id <- meta[["sample id"]]

  meta <- meta %>%
    mutate(
      age_years = safe_numeric(age_years),
      grade = na_clean(grade),
      grade = ifelse(grade == "4=Indeterminate", "Indeterminate", grade),
      response_primary = na_clean(pathologic_response_pcr_rd),
      response_rcb = na_clean(pathologic_response_rcb_class),
      response_rcb_binary = case_when(
        response_rcb == "RCB-0/I" ~ "RCB-0/I",
        response_rcb %in% c("RCB-II", "RCB-III") ~ "RCB-II/III",
        TRUE ~ NA_character_
      ),
      response_binary = case_when(
        response_primary == "pCR" ~ 1L,
        response_primary == "RD" ~ 0L,
        TRUE ~ NA_integer_
      )
    )

  keep <- unique(c(keep_fields, "response_primary", "response_rcb", "response_rcb_binary", "response_binary"))
  meta %>%
    mutate(across(everything(), na_clean)) %>%
    mutate(age_years = suppressWarnings(as.numeric(age_years))) %>%
    dplyr::select(any_of(keep))
}

subset_primary_cohort <- function(meta, filter_field, filter_value, endpoint_field, positive, negative) {
  meta %>%
    filter(.data[[filter_field]] == filter_value) %>%
    filter(.data[[endpoint_field]] %in% c(positive, negative))
}

subset_secondary_cohort <- function(meta, filter_field, filter_value, endpoint_field, positive_labels, negative_labels) {
  meta %>%
    filter(.data[[filter_field]] == filter_value) %>%
    filter(.data[[endpoint_field]] %in% c(positive_labels, negative_labels))
}

write_table_safe <- function(df, path) {
  ensure_dir(dirname(path))
  data.table::fwrite(df, path, sep = "\t", quote = FALSE, na = "NA")
}

save_rds_safe <- function(obj, path) {
  ensure_dir(dirname(path))
  saveRDS(obj, path)
}

theme_report <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

make_collection_gene_sets <- function(collection, subcollection = NULL, species = "Homo sapiens") {
  tbl <- tryCatch(
    {
      if (is.null(subcollection)) {
        msigdbr::msigdbr(species = species, collection = collection)
      } else {
        msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection)
      }
    },
    error = function(e) {
      if (is.null(subcollection)) {
        msigdbr::msigdbr(species = species, category = collection)
      } else {
        msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection)
      }
    }
  )
  tbl <- tbl %>%
    dplyr::select(gs_name, gene_symbol) %>%
    distinct()
  split(tbl$gene_symbol, tbl$gs_name)
}

load_custom_immune_signatures <- function(path) {
  sig_tbl <- readr::read_tsv(path, show_col_types = FALSE)
  split(sig_tbl$gene_symbol, sig_tbl$signature_name)
}

compute_ssgsea_scores <- function(expr_gene, gene_sets, min_size = 5) {
  gene_sets <- gene_sets[vapply(gene_sets, function(gs) sum(unique(gs) %in% rownames(expr_gene)) >= min_size, logical(1))]
  if (length(gene_sets) == 0) {
    stop("No gene sets passed the minimum size filter.")
  }

  expr_gene <- as.matrix(expr_gene)

  res <- tryCatch(
    {
      param <- GSVA::ssgseaParam(exprData = expr_gene, geneSets = gene_sets, normalize = TRUE)
      GSVA::gsva(param, verbose = FALSE)
    },
    error = function(e) {
      GSVA::gsva(
        expr = expr_gene,
        gset.idx.list = gene_sets,
        method = "ssgsea",
        kcdf = "Gaussian",
        abs.ranking = FALSE,
        verbose = FALSE
      )
    }
  )

  as.matrix(res)
}

scale_rows <- function(mat) {
  t(scale(t(mat)))
}

compute_pr_auc <- function(truth, prob) {
  ord <- order(prob, decreasing = TRUE)
  truth <- truth[ord]
  prob <- prob[ord]

  tp <- cumsum(truth == 1)
  fp <- cumsum(truth == 0)
  precision <- tp / pmax(tp + fp, 1)
  recall <- tp / sum(truth == 1)

  precision <- c(1, precision)
  recall <- c(0, recall)
  sum((recall[-1] - recall[-length(recall)]) * precision[-1])
}

categorical_association_test <- function(x, y, simulate_b = 10000) {
  tab <- table(x, y, useNA = "ifany")
  if (!all(dim(tab) > 1)) {
    return(list(test = "NA", p_value = NA_real_))
  }

  chi_fit <- suppressWarnings(chisq.test(tab))
  if (any(chi_fit$expected < 5)) {
    if (nrow(tab) == 2 && ncol(tab) == 2) {
      return(list(test = "fisher", p_value = fisher.test(tab)$p.value))
    }
    sim_fit <- chisq.test(tab, simulate.p.value = TRUE, B = simulate_b)
    return(list(test = "chisq_simulated", p_value = sim_fit$p.value))
  }

  list(test = "chisq", p_value = chi_fit$p.value)
}

make_ranked_stats_unique <- function(stats_named, jitter_scale = 1e-9) {
  stats_named <- stats_named[!is.na(stats_named)]
  if (length(stats_named) == 0) {
    return(stats_named)
  }

  ord <- order(-stats_named, names(stats_named))
  stats_named <- stats_named[ord]
  offsets <- seq_along(stats_named) * jitter_scale
  adjusted <- stats_named - offsets
  names(adjusted) <- names(stats_named)
  adjusted
}

close_plot_device <- function() {
  invisible(dev.off())
}

extract_top_variable_genes <- function(expr_gene, n = 50) {
  vars <- apply(expr_gene, 1, var, na.rm = TRUE)
  names(sort(vars, decreasing = TRUE))[seq_len(min(n, length(vars)))]
}

discretize_equal_frequency <- function(x, n_bins = 5) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    return(rep(NA_character_, length(x)))
  }

  unique_vals <- unique(x[!is.na(x)])
  if (length(unique_vals) < 2) {
    return(rep(NA_character_, length(x)))
  }

  n_bins <- max(2L, min(as.integer(n_bins), length(unique_vals)))
  probs <- seq(0, 1, length.out = n_bins + 1L)
  breaks <- unique(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 8))

  if (length(breaks) < 2) {
    return(rep(NA_character_, length(x)))
  }

  as.character(cut(x, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE))
}

compute_empirical_mutual_information <- function(x, y, n_bins = 5) {
  y <- na_clean(y)
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) < 10) {
    return(NA_real_)
  }

  x_disc <- discretize_equal_frequency(x[keep], n_bins = n_bins)
  keep_disc <- !is.na(x_disc)
  if (sum(keep_disc) < 10) {
    return(NA_real_)
  }

  x_disc <- x_disc[keep_disc]
  y_keep <- y[keep][keep_disc]

  if (length(unique(x_disc)) < 2 || length(unique(y_keep)) < 2) {
    return(0)
  }

  joint <- table(x_disc, y_keep)
  pxy <- joint / sum(joint)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  idx <- which(pxy > 0, arr.ind = TRUE)

  sum(vapply(seq_len(nrow(idx)), function(i) {
    row_i <- idx[i, 1]
    col_i <- idx[i, 2]
    pxy_val <- pxy[row_i, col_i]
    pxy_val * log2(pxy_val / (px[row_i] * py[col_i]))
  }, numeric(1)))
}

extract_top_mi_genes <- function(expr_gene, response, n = 500, n_bins = 5) {
  expr_gene <- as.matrix(expr_gene)
  scores <- apply(expr_gene, 1, compute_empirical_mutual_information, y = response, n_bins = n_bins)
  scores <- scores[!is.na(scores)]
  scores <- sort(scores, decreasing = TRUE)
  top_n <- min(as.integer(n), length(scores))
  list(
    scores = scores,
    genes = names(scores)[seq_len(top_n)]
  )
}

safe_enrichplot <- function(enrich_obj, path, title_text) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(invisible(FALSE))
  }
  p <- dotplot(enrich_obj, showCategory = 15) + ggtitle(title_text)
  ggsave(path, p, width = 10, height = 7, dpi = 300)
  invisible(TRUE)
}

coerce_categorical_na <- function(df, columns) {
  for (col in columns) {
    if (col %in% names(df)) {
      vals <- as.character(df[[col]])
      vals[is.na(vals)] <- "NA"
      df[[col]] <- vals
    }
  }
  df
}
