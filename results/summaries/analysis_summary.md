# GSE25066 Analysis Summary

## Cohort

# A tibble: 3 × 6
  cohort                n   pcr    rd rcb_0_i rcb_ii_iii
  <chr>             <dbl> <dbl> <dbl>   <dbl>      <dbl>
1 all_samples         508    57   249      NA         NA
2 her2neg_primary     294    55   239      NA         NA
3 her2neg_secondary   289    NA    NA      83        206

## Probe-to-Gene Mapping

# A tibble: 1 × 5
  probe_count annotated_probe_count unique_gene_count collapsed_probe_count
        <dbl>                 <dbl>             <dbl>                 <dbl>
1       22283                 21056             13039                 13039
# ℹ 1 more variable: dropped_probe_count <dbl>

## EDA Highlights

# A tibble: 8 × 3
  variable                             test                  p_value
  <chr>                                <chr>                   <dbl>
1 er_status_ihc_esr1_for indeterminate chisq           0.00000000179
2 grade                                chisq_simulated 0.0001000    
3 pam50_class                          chisq_simulated 0.0001000    
4 age_years                            wilcox          0.101        
5 clinical_t_stage                     chisq_simulated 0.297        
6 clinical_nodal_status                chisq           0.680        
7 clinical_ajcc_stage                  chisq_simulated 0.697        
8 source                               chisq           0.925        

## Biological Interpretation

Top GO BP terms up in pCR:

# A tibble: 10 × 12
   ID    Description GeneRatio BgRatio RichFactor FoldEnrichment zScore   pvalue
   <chr> <chr>       <chr>     <chr>        <dbl>          <dbl>  <dbl>    <dbl>
 1 GO:0… chromosome… 56/430    337/11…      0.166           4.63  13.0  1.71e-22
 2 GO:0… nuclear di… 56/430    338/11…      0.166           4.62  13.0  1.98e-22
 3 GO:0… nuclear ch… 48/430    249/11…      0.193           5.37  13.5  3.20e-22
 4 GO:0… mitotic nu… 46/430    234/11…      0.197           5.48  13.3  1.07e-21
 5 GO:0… mitotic si… 38/430    161/11…      0.236           6.58  13.7  5.15e-21
 6 GO:0… organelle … 56/430    374/11…      0.150           4.17  12.0  2.96e-20
 7 GO:0… sister chr… 40/430    195/11…      0.205           5.72  12.8  1.10e-19
 8 GO:1… microtubul… 30/430    143/11…      0.210           5.85  11.2  2.54e-15
 9 GO:0… mitotic ce… 46/430    355/11…      0.130           3.61   9.64 2.42e-14
10 GO:0… cell cycle… 51/430    438/11…      0.116           3.25   9.24 6.58e-14
# ℹ 4 more variables: p.adjust <dbl>, qvalue <dbl>, geneID <chr>, Count <dbl>

Top KEGG terms up in pCR:

# A tibble: 4 × 12
  ID     Description GeneRatio BgRatio RichFactor FoldEnrichment zScore   pvalue
  <chr>  <chr>       <chr>     <chr>        <dbl>          <dbl>  <dbl>    <dbl>
1 KEGG_… KEGG_CELL_… 21/142    115/43…     0.183            5.55   9.12 6.38e-11
2 KEGG_… KEGG_P53_S… 9/142     62/4317     0.145            4.41   4.99 1.61e- 4
3 KEGG_… KEGG_OOCYT… 11/142    95/4317     0.116            3.52   4.58 2.41e- 4
4 KEGG_… KEGG_CHEMO… 14/142    171/43…     0.0819           2.49   3.66 1.31e- 3
# ℹ 4 more variables: p.adjust <dbl>, qvalue <dbl>, geneID <chr>, Count <dbl>

## Modeling

Nested CV performance:

# A tibble: 195 × 9
   feature_set   model    `repeat`  fold auroc auprc balanced_accuracy threshold
   <chr>         <chr>       <dbl> <dbl> <dbl> <dbl>             <dbl>     <dbl>
 1 clinical_only ridge_l…        1     1 0.714 0.301             0.652     0.524
 2 clinical_only ridge_l…        1     2 0.752 0.336             0.683     0.522
 3 clinical_only ridge_l…        1     3 0.748 0.421             0.616     0.496
 4 clinical_only ridge_l…        1     4 0.722 0.446             0.632     0.515
 5 clinical_only ridge_l…        1     5 0.877 0.495             0.833     0.509
 6 clinical_only ridge_l…        2     1 0.803 0.543             0.707     0.506
 7 clinical_only ridge_l…        2     2 0.775 0.367             0.732     0.498
 8 clinical_only ridge_l…        2     3 0.718 0.388             0.707     0.502
 9 clinical_only ridge_l…        2     4 0.816 0.515             0.726     0.502
10 clinical_only ridge_l…        2     5 0.727 0.282             0.673     0.517
# ℹ 185 more rows
# ℹ 1 more variable: best_params <chr>

Source-transfer performance:

# A tibble: 26 × 9
   feature_set   model    train_source test_source auroc auprc balanced_accuracy
   <chr>         <chr>    <chr>        <chr>       <dbl> <dbl>             <dbl>
 1 clinical_only ridge_l… MDACC        ISPY        0.699 0.274             0.660
 2 clinical_only ridge_l… ISPY         MDACC       0.678 0.326             0.631
 3 clinical_only elastic… MDACC        ISPY        0.702 0.281             0.660
 4 clinical_only elastic… ISPY         MDACC       0.690 0.336             0.653
 5 clinical_only linear_… MDACC        ISPY        0.709 0.285             0.645
 6 clinical_only linear_… ISPY         MDACC       0.696 0.346             0.627
 7 clinical_only xgboost  MDACC        ISPY        0.687 0.385             0.645
 8 clinical_only xgboost  ISPY         MDACC       0.737 0.354             0.670
 9 bi_only       ridge_l… MDACC        ISPY        0.670 0.258             0.708
10 bi_only       ridge_l… ISPY         MDACC       0.742 0.335             0.696
# ℹ 16 more rows
# ℹ 2 more variables: threshold <dbl>, best_params <chr>

PASNet tuning summary:

# A tibble: 1 × 6
  model  feature_set  auroc_mean baseline_best_auroc target_auroc note          
  <chr>  <chr>             <dbl>               <dbl>        <dbl> <chr>         
1 pasnet gene_pathway      0.753               0.778        0.758 PASNet did no…

Top predictive signatures:

# A tibble: 195 × 4
   feature                                mean_abs_importance feature_set  model
   <chr>                                                <dbl> <chr>        <chr>
 1 er_status_ihc_esr1_for indeterminate_P              0.205  clinical_on… ridg…
 2 er_status_ihc_esr1_for indeterminate_N              0.177  clinical_on… ridg…
 3 grade_3                                             0.151  clinical_on… ridg…
 4 grade_2                                             0.143  clinical_on… ridg…
 5 age_years                                           0.0827 clinical_on… ridg…
 6 clinical_t_stage_T4                                 0.0527 clinical_on… ridg…
 7 clinical_nodal_status_N0                            0.0386 clinical_on… ridg…
 8 grade_1                                             0.0239 clinical_on… ridg…
 9 clinical_ajcc_stage_IIB                             0.0222 clinical_on… ridg…
10 clinical_t_stage_T1                                 0.0213 clinical_on… ridg…
# ℹ 185 more rows

