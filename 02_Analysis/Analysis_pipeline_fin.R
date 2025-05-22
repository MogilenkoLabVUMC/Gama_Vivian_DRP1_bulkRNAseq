###############################################################################
##  Analysis pipeline – human bulk RNA-seq (refactored 2025-05-19)          ##
###############################################################################
##  Key improvements ----------------------------------------------------------
##  • robust path handling (here::here) and helper-sourcing
##  • automatic package checks / installs
##  • all constants collected in one config section
##  • species-aware GSEA (HS vs MM) + SynGO support
##  • fixes: helper sourcing, get_db_plot_params, empty first PDF page,
##           unreadable volcano labels, hclust NaN, heat-map guards, etc.
###############################################################################

# -------------------------------------------------------------------- #
# 0.  Configuration – edit in ONE place                                #
# -------------------------------------------------------------------- #
config <- list(
  ## raw data
  counts_file   = "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/sorted_counts_matrix.txt",
  metadata_file = "03_Results/01_Preprocessing/04_FeatureCounts/count_matrices_fc/metadata.csv",
  ## output roots
  out_root      = "03_Results/02_Analysis",
  helper_root   = "01_Scripts/GSEA_module",   # <<-- external module mount
  ## analysis parameters
  p_cutoff      = 0.05,
  fc_cutoff     = 2,
  calcium_genes = c(
      "NEURONATIN","CACNG3","CACNA1C","CACNA1S","ATP2A1",
      "RYR1","MYLK3","CASR","VDR","STIM1","STIM2",
      "ORAI1","CALB1","CALR"),
  ## SynGO data (relative to repo root)
  syngo_dir     = "00_Data/SynGO_bulk_20231201",
  syngo_ns      = "CC"                       # GO cellular-component ontology
)

# -------------------------------------------------------------------- #
# 1.  Packages & helper sourcing                                       #
# -------------------------------------------------------------------- #
required_pkgs <- c(
  "edgeR","limma","dplyr","ggplot2","pheatmap","RColorBrewer","viridis",
  "reshape2","VennDiagram","grid","UpSetR","WGCNA",
  "msigdbr","clusterProfiler","fgsea","org.Hs.eg.db",
  "patchwork","here")

for (p in required_pkgs){
  if (!requireNamespace(p, quietly = TRUE)){
    message(sprintf("• installing %s …", p))
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

## helper-sourcing -----------------------------------------------------
source_if_present <- function(...) {
  path <- here::here(...)
  if (file.exists(path)) {
    source(path, echo = FALSE)
  } else {
    warning("helper not found → ", path)
  }
}

## GSEA module helpers (keep list centralised)
gsea_helpers <- c(
  "scripts/custom_minimal_theme.R",
  "scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R",
  "scripts/GSEA/GSEA_plotting/gsea_dotplot.R",
  "scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R",
  "scripts/GSEA/GSEA_plotting/gsea_barplot.R",
  "scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R",
  "scripts/GSEA/GSEA_plotting/gsea_heatmap.R",
  "scripts/GSEA/GSEA_processing/run_gsea.R",
  "scripts/GSEA/GSEA_processing/run_gsea_analysis.R",
  "scripts/DE/volcano_helpers.R"
)

for (h in gsea_helpers){
  source_if_present(config$helper_root, h)
}

## other single helpers
source_if_present("01_Scripts/R_scripts/read_count_matrix.R")
source_if_present("01_Scripts/GSEA_module/scripts/DE/plot_standard_volcano.R")
source_if_present("/workspaces/GVDRP1/01_Scripts/GSEA_module/scripts/DE/create_fc_b_plot.R")
source_if_present("/workspaces/GVDRP1/01_Scripts/GSEA_module/scripts/DE/create_MD_plot.R")
source_if_present("01_Scripts/GSEA_module/scripts/DE/plotPCA.R")
source_if_present("/workspaces/GVDRP1/01_Scripts/R_scripts/generate_vertical_volcanos.R")
source_if_present("/workspaces/GVDRP1/01_Scripts/R_scripts/run_syngo_gsea.R")
source_if_present("/workspaces/GVDRP1/01_Scripts/GSEA_module/scripts/GSEA/GSEA_plotting/plot_all_gsea_results.R")


# -------------------------------------------------------------------- #
# 2.  Read & pre-process data                                          #
# -------------------------------------------------------------------- #
DGE <- process_rnaseq_data(config$counts_file, config$metadata_file, annotate = FALSE)

## design / contrasts --------------------------------------------------
DGE$samples$group <- with(DGE$samples,
                          interaction(days, genotype, sep = "_", drop = TRUE))
design <- model.matrix(~0 + group, data = DGE$samples)
colnames(design) <- levels(DGE$samples$group)

contrasts <- makeContrasts(
    # 1. Pair-wise comparisons
  # Question: What is the baseline difference between the mutant vs control?
  ## mutation vs ctrl
  G32A_vs_Ctrl_D35   = D35_G32A - D35_Control,
  R403C_vs_Ctrl_D35  = D35_R403C - D35_Control,
  G32A_vs_Ctrl_D65   = D65_G32A - D65_Control,
  R403C_vs_Ctrl_D65  = D65_R403C - D65_Control,

    # 2. Maturation effects within each genotype
  # Question: How does maturation affect each genotype?
  ## maturation
  Time_Ctrl          = D65_Control - D35_Control,
  Time_G32A          = D65_G32A   - D35_G32A,
  Time_R403C         = D65_R403C  - D35_R403C,

    # 3. Interaction effects (difference-in-difference)
  # Question: Do the mutations alter the normal maturation trajectory?
  ## interaction
  Maturation_G32A_specific = (D65_G32A  - D35_G32A)  - (D65_Control - D35_Control),
  Maturation_R403C_specific = (D65_R403C - D35_R403C) - (D65_Control - D35_Control),
  levels = design)

# -------------------------------------------------------------------- #
# 3.  DEG modelling (voomLmFit preferred)                              #
# -------------------------------------------------------------------- #
## keep.EList=TRUE ➜ the returned object contains $EList with log-CPM + weights

fit <- edgeR::voomLmFit(
         DGE,
         design,
         sample.weights = FALSE,   # or TRUE if you want quality weights
         keep.EList     = TRUE)

## unpack the voom log-CPM for QC plots
v <- fit$EList                    # same structure as `voom()` output

## add contrasts & empirical Bayes
fit <- contrasts.fit(fit, contrasts) |>
       eBayes(robust = TRUE)

de_results <- decideTests(fit, p.value = config$p_cutoff)

# ------------------------------------------------------------
# save DEGs lists
# ------------------------------------------------------------
deg_dir <- "/workspaces/GVDRP1/03_Results/02_Analysis/DE_results"
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2.  Loop over every contrast in your 'contrasts' matrix
# ------------------------------------------------------------
for (co in colnames(contrasts)) {

  ## limma::topTable – all rows, sorted by t
  tt <- limma::topTable(fit,
                        coef      = co,
                        number    = Inf,     # keep everything
                        sort.by   = "t")     # descending t-statistic

  ## write with gene IDs as the first column
  fname <- file.path(deg_dir, paste0(co, "_DE_results.csv"))
  write.csv(tt, file = fname, row.names = TRUE)

  message("saved: ", fname)
}

# -------------------------------------------------------------------- #
# 4.  QC: sample correlation heat-map & MDS                            #
# -------------------------------------------------------------------- #
qc_dir <- here::here(config$out_root, "Plots/General")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

## ordered sample vector ------------------------------------------------
ord  <- with(DGE$samples,
             order(days, factor(genotype, c("Control","G32A","R403C"))))
ordered_samples <- rownames(DGE$samples)[ord]

logCPM <- cpm(DGE, log = TRUE)
ordered_cor <- cor(logCPM)[ordered_samples, ordered_samples]

annot <- DGE$samples[ordered_samples, c("genotype","days")]
ann_colors <- list(
  genotype = c(Control = "#1B9E77", G32A = "#D95F02", R403C = "#7570B3"),
  days     = c(D35     = "#E7298A", D65 = "#66A61E"))

pdf(file.path(qc_dir, "sample_correlation_heatmap_ordered.pdf"), 12, 10)
pheatmap(ordered_cor,
         annotation_row  = annot, annotation_col = annot,
         annotation_colors = ann_colors,
         color = colorRampPalette(brewer.pal(9,"YlOrBr"))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = NA)
dev.off()

## MDS plot -------------------------------------------------------------
pdf(file.path(qc_dir,"MDS_plot.pdf"), 10, 8)
plotMDS(v, col = ann_colors$genotype[annot$genotype],
        pch = ifelse(annot$days=="D35", 16, 17),
        cex = 1.4)
leg_lbl <- unique(paste(annot$genotype, annot$days))
leg_col <- ann_colors$genotype[ sub(" .*", "", leg_lbl) ]   # map back to genotype
leg_pch <- ifelse(grepl("D35$", leg_lbl), 16, 17)

legend("topright", legend = leg_lbl,
       col = leg_col, pch = leg_pch,
       bty = "n", cex = 1.0, pt.cex = 1.4)
dev.off()

# -------------------------------------------------------------------- #
# 5.  Volcanoes & DEG numbers                                          #
# -------------------------------------------------------------------- #
# P-value volcanoes
deg_dir      <- here::here(config$out_root, "DE_results")
volcano_dir  <- here::here(config$out_root, "Plots/Volcano/p")
dir.create(deg_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

## helper to make + save volcanoes -------------------------------------
make_all_volcanoes <- function(res_table, name){
  pdf(file.path(volcano_dir, paste0(name,"_standard.pdf")), 8, 7)
  print(create_standard_volcano(
          res_table, 
          p_cutoff = config$p_cutoff, 
          fc_cutoff = config$fc_cutoff,
          decision_by = "p",
          label_method = "top", 
          highlight_gene = config$calcium_genes,
          x_breaks = 2,
          title = name))
  dev.off()
}

for (co in colnames(contrasts)){
  tbl <- topTable(fit, coef = co, number = Inf)
  write.csv(tbl, file = file.path(deg_dir, paste0(co,"_results.csv")))
  make_all_volcanoes(tbl, co)
}

# FDR-value volcanos
deg_dir      <- here::here(config$out_root, "DE_results")
volcano_dir  <- here::here(config$out_root, "Plots/Volcano/fdr")
dir.create(deg_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

## helper to make + save volcanoes -------------------------------------
make_all_volcanoes <- function(res_table, name){
  pdf(file.path(volcano_dir, paste0(name,"_standard.pdf")), 8, 7)
  print(create_standard_volcano(
          res_table, 
          p_cutoff = 0.1, 
          fc_cutoff = config$fc_cutoff,
          decision_by = "fdr",
          label_method = "sig", 
          highlight_gene = config$calcium_genes,
          x_breaks = 2,
          title = name))
  dev.off()
}

for (co in colnames(contrasts)){
  tbl <- topTable(fit, coef = co, number = Inf)
  write.csv(tbl, file = file.path(deg_dir, paste0(co,"_results.csv")))
  make_all_volcanoes(tbl, co)
}

## Vertical volcanoes -------------------------------------
volcano_dir_vert <- here::here(config$out_root, "Plots/Volcano")
dir.create(volcano_dir_vert, recursive = TRUE, showWarnings = FALSE)

# Create a list to store all contrast tables
contrast_tables <- list()
for (co in colnames(contrasts)) {
  contrast_tables[[co]] <- topTable(fit, coef = co, number = Inf)
}

# Generate all the vertical volcano plot combinations
# This will create both versions - with and without calcium gene highlighting
generate_vertical_volcano_sets(contrast_tables, config, highlight_calcium = TRUE)

## ---- MD plot -----------------------------------------------------------
deg_dir      <- here::here(config$out_root, "DE_results")
volcano_dir  <- here::here(config$out_root, "Plots/Volcano/MD")
dir.create(deg_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

# Enhanced MD plot function
make_md_plot <- function(fit, coef, name, status = NULL, highlight_genes = NULL) {
  pdf(file.path(volcano_dir, paste0(name, "_MDplot.pdf")), 7, 6)
  
  # Get the DE results once
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  
  # Create the plot with our new function
  gg <- create_MD_plot(
    fit = fit,
    coef = coef,
    de_results = tt,
    fc_cutoff = 1,
    fdr_cutoff = 0.05,
    top_n = 5,
    highlight_gene = highlight_genes,
    label_method = "top",
    title = paste("MD plot:", name),
    show_grid = FALSE
  )
  
  print(gg)
  dev.off()
}

# Use the enhanced function in your analysis
for (co in colnames(contrasts)) {
  make_md_plot(fit, coef = co, name = co, highlight_genes = config$calcium_genes)
}

## ---- FC vs B -----------------------------------------------------------
deg_dir      <- here::here(config$out_root, "DE_results")
volcano_dir  <- here::here(config$out_root, "Plots/Volcano/FC-B")
dir.create(deg_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

make_fc_vs_B <- function(top, name, fc_cutoff = config$fc_cutoff, highlight_genes = NULL) {
  pdf(file.path(volcano_dir, paste0(name, "_FC_vs_B.pdf")), 7, 6)
  
  # Create the B vs FC plot using our new function
  gg <- create_B_FC_plot(
    top,
    fc_cutoff = fc_cutoff,
    B_cutoff = 0,  # Traditional threshold for B-statistic
    top_n = 5,
    highlight_gene = highlight_genes,
    label_method = "top",
    title = paste("log2FC vs B:", name),
    show_grid = FALSE  # No grid lines as per your request
  )
  
  print(gg)
  dev.off()
}

for (co in colnames(contrasts)) {
  top <- limma::topTable(fit, coef = co, number = Inf)
  make_fc_vs_B(top, name = co, highlight_genes = config$calcium_genes)
}


## DEG counts ----------------------------------------------------------
## -------- prepare data ---------------------------------------------
contrast_order <- c("G32A_vs_Ctrl_D35","R403C_vs_Ctrl_D35",
                    "G32A_vs_Ctrl_D65","R403C_vs_Ctrl_D65",
                    "Time_Ctrl","Time_G32A","Time_R403C",
                    "Maturation_G32A_specific","Maturation_R403C_specific")

deg_counts <- data.frame(
  Contrast = colnames(contrasts),
  Up   = colSums(de_results  > 0),
  Down = colSums(de_results  < 0)) |>
  transform(Contrast = factor(Contrast, levels = contrast_order))

deg_long <- tidyr::pivot_longer(deg_counts,
                                c(Up,Down),
                                names_to = "Direction",
                                values_to = "Count")

## -------- plot ------------------------------------------------------
pdf(file.path(qc_dir,"DE_gene_counts.pdf"), 10, 6, onefile = TRUE)
ggplot(deg_long,
       aes(Contrast, Count, fill = Direction)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_text(aes(label = Count),
            vjust = -.2,
            position = position_dodge(width = .8),
            size = 3) +
  scale_fill_manual(values = c(Up = "#D55E00", Down = "#0072B2")) + # Okabe-Ito
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()) +
  labs(y = "gene count", x = NULL, fill = "")
dev.off()



# Create a list of significant genes for each contrast
sig_genes_list <- list()
for (contrast in colnames(contrasts)) {
  sig_genes_list[[contrast]] <- rownames(fit$coefficients)[which(de_results[, contrast] != 0)]
}

# Create and save UpSet plot
upset_plot <- upset(fromList(sig_genes_list), order.by = "freq", nsets = length(sig_genes_list))
pdf(file.path(qc_dir, "UpSet_plot_all_contrasts.pdf"), width = 12, height = 8)
print(upset_plot)
dev.off()



# -------------------------------------------------------------------- #
# 6.  Generic MSigDB-based GSEA                                        #
# -------------------------------------------------------------------- #
gsea_root <- here::here(config$out_root,"Plots/GSEA")
dir.create(gsea_root, recursive = TRUE, showWarnings = FALSE)

## wrapper  -------------------------
run_gsea_hsmm <- function(tbl, contrast, species){
  results <- run_gsea_analysis(
      de_table     = tbl,
      analysis_name= contrast,
      rank_metric  = "t",
      species      = species,
      n_pathways   = 30,
      padj_cutoff  = 0.05,
      output_dir   = file.path(gsea_root, contrast),
      sample_annotation = DGE$samples[,c("genotype","days"),drop=FALSE],
      sample_order = ordered_samples,
      helper_root  = config$helper_root,
      save_plots   = TRUE)
  
  return(results)  # Explicitly return the results
}

# Create a list to store all GSEA results
all_gsea_results <- list()

# Run GSEA for each contrast and save results
for (co in colnames(contrasts)) {
  message("==== Working on contrast: ", co, " ====")
  tbl <- topTable(fit, coef = co, number = Inf)

  safe_res <- tryCatch({
    run_gsea_hsmm(tbl, co, species = "Homo sapiens")
  }, error = function(e) {
    message("❗ ERROR in GSEA for contrast ", co, ": ", e$message)
    message("👀 Calling traceback():")
    traceback()
    NULL
  })

  all_gsea_results[[co]] <- safe_res
}

# Then call plot_all_gsea_results with these parameters for each contrast
for (co in colnames(contrasts)) {
  message("==== GSEA: ", co, " ====")
  this_res_list <- all_gsea_results[[co]]
  if (is.null(this_res_list) || !length(this_res_list)) {
    message("Skipping ", co, " (no results).")
    next
  }
  
  plot_all_gsea_results(
  gsea_list = this_res_list,
  analysis_name = co,
  out_root = gsea_root,
  n_pathways = 30,
  padj_cutoff = 0.05
 )

}


# -------------------------------------------------------------------- #
# 7.  SynGO GSEA                                                       #
# -------------------------------------------------------------------- #
source("/workspaces/GVDRP1/01_Scripts/R_scripts/run_syngo_gsea.R")

# function to prepare a gmt list 
syngo_gmt <- function(syngo_dir, namespace = "CC") {
  requireNamespace("readxl", quietly = TRUE)
  # read xlsx files
  ann <- readxl::read_xlsx(file.path(syngo_dir, "syngo_annotations.xlsx"))
  ont <- readxl::read_xlsx(file.path(syngo_dir, "syngo_ontologies.xlsx"))

  ann <- subset(ann, go_domain == namespace)

  ## ── TERM2GENE (required) ───────────────────────────────────────────
  term2gene <- ann[, c("go_id", "hgnc_symbol")]
  colnames(term2gene) <- c("gs_name", "gene_symbol")

  ## ── TERM2NAME (optional but recommended) ───────────────────────────
  term2name <- unique(ont[, c("id", "name")])
  colnames(term2name) <- c("gs_name", "description")

  ## strip trailing “ (GO:########)” → cleaner axis labels
  term2name$description <- sub(" \\(GO:[0-9]+\\)$", "", term2name$description)

  list(T2G = term2gene, T2N = term2name)
}

# execute gmt list preparation 
syngo_lists <- syngo_gmt(config$syngo_dir, config$syngo_ns)

# Create a list to store all GSEA results
syngo_gsea_results <- list()

# Run GSEA for each contrast
for (co in colnames(contrasts)) {
  tbl <- topTable(fit, coef = co, number = Inf)
  
  # Close any lingering graphic devices before starting new contrast
  while (dev.cur() > 1) dev.off()
  
  # Capture the result
  syngo_gsea_results[[co]] <- run_syngo_gsea(
    tbl, 
    co,
    T2G = syngo_lists$T2G,
    T2N = syngo_lists$T2N,
    sample_annotation = annot  # Add if you have it
  )
  
  # Ensure all devices are closed after each iteration
  while (dev.cur() > 1) dev.off()
}


###############################################################################
##  SynGO enrichment of intersecting DEG signatures
##  – baseline-9 genes (shared mutant baseline)
##  – mat-38 genes    (shared mutant maturation response)
###############################################################################
# -------------------------------------------------------------------------
#  synGO_enrich_symbol()  –  overlap-aware ORA with synonym resolution
# -------------------------------------------------------------------------
synGO_enrich_symbol <- function(query, name, universe,
                                syngo_T2G, syngo_T2N, syngo_genes_xlsx,
                                namespace_tag = "CC", qcut = .2,
                                out_dir = out_root) {

  requireNamespace("readxl",           quietly = TRUE)
  requireNamespace("clusterProfiler",  quietly = TRUE)

  # 1. synonym table --------------------------------------------------
  canon_tbl <- readxl::read_xlsx(syngo_genes_xlsx,
                                 col_types = "text")[, c("hgnc_symbol",
                                                         "hgnc_synonyms")] |>
               tidyr::separate_rows(hgnc_synonyms, sep = ",\\s*") |>
               dplyr::transmute(canonical = hgnc_symbol,
                                alias      = dplyr::coalesce(hgnc_synonyms,
                                                             hgnc_symbol)) |>
               dplyr::distinct()

  to_canonical <- function(x)
    dplyr::left_join(data.frame(alias = x), canon_tbl,
                     by = "alias")$canonical

  # 2. TERM2GENE fixed-length join -----------------------------------
  syngo_T2G_can <- syngo_T2G |>
    dplyr::left_join(canon_tbl, by = c(gene_symbol = "alias")) |>
    dplyr::transmute(gs_name,
                     gene_symbol = canonical) |>
    dplyr::filter(!is.na(gene_symbol)) |>
    dplyr::distinct()

  print('Intersections with SynGo')
  print(table(query %in% syngo_T2G_can$gene_symbol))

  # 3. canonicalise query / universe ---------------------------------
  canon_query    <- unique(na.omit(to_canonical(query)))
  canon_universe <- unique(na.omit(to_canonical(universe)))

  canon_query    <- intersect(canon_query,    syngo_T2G_can$gene_symbol)
  canon_universe <- intersect(canon_universe, syngo_T2G_can$gene_symbol)

  if (length(canon_query) < 3) {
    warning(name, ": only ", length(canon_query),
            " genes overlap SynGO after canonicalisation – skipping.")
    return(invisible(NULL))
  }

  # 4. ORA ------------------------------------------------------------
  enr <- clusterProfiler::enricher(
           gene          = canon_query,
           universe      = canon_universe,
           TERM2GENE     = syngo_T2G_can,
           TERM2NAME     = syngo_T2N,
           pAdjustMethod = "BH",
           qvalueCutoff  = qcut,
           minGSSize     = 2)

  if (is.null(enr) || nrow(enr@result) == 0) {
    message(name, ": no SynGO term < ", qcut)
    return(invisible(NULL))
  }

  # 5. output ---------------------------------------------------------
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(enr@result,
            file = file.path(out_dir, paste0(name, "_SynGO_enrich.csv")),
            row.names = FALSE)

  pdf(file.path(out_dir, paste0(name, "_SynGO_dotplot.pdf")), 7, 5)
  print(
    clusterProfiler::dotplot(enr, showCategory = 15) +
      ggplot2::labs(title = sprintf("%s – SynGO (%s)", name, namespace_tag)) +
      ggplot2::theme_minimal(base_size = 11)
  )
  dev.off()

  invisible(enr)
}

## ------------------------------------------------------------------ ##
## 4  run enrichment                                                  ##
## ------------------------------------------------------------------ ##
universe_all <- rownames(fit)          # every gene tested in limma

## path to the original SynGO gene table
syngo_gene_file <- file.path(config$syngo_dir, "syngo_genes.xlsx")

enr_base9 <- synGO_enrich_symbol(baseline9, "baseline9",
                                 universe_all,
                                 syngo_lists$T2G, syngo_lists$T2N,
                                 syngo_gene_file)

enr_mat38 <- synGO_enrich_symbol(mat38,     "mat38",
                                 universe_all,
                                 syngo_lists$T2G, syngo_lists$T2N,
                                 syngo_gene_file)


# -------------------------------------------------------------------- #
# 8.  Calcium genes focus                                              #
# -------------------------------------------------------------------- #
calc_dir <- here::here(config$out_root,"Calcium_genes")
dir.create(calc_dir, recursive = TRUE, showWarnings = FALSE)

present <- intersect(config$calcium_genes, rownames(logCPM))
if (length(present)){
  mat <- logCPM[present, ordered_samples]

  ## clean NA / Inf for clustering
  mat[!is.finite(mat)] <- NA
  mat <- t(scale(t(mat)))    # row-scale
  mat[is.na(mat)] <- 0

  pdf(file.path(calc_dir,"calcium_expr_heatmap.pdf"), 12, .4*length(present)+3)
  pheatmap(mat, annotation_col = annot, annotation_colors = ann_colors,
           color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
           cluster_cols = FALSE, border_color = NA)
  dev.off()
}



# -------------------------------------------------------------------- #
# 9.  Summary tables & markdown report                                 #
# -------------------------------------------------------------------- #
sum_dir <- here::here(config$out_root,"Summary")
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(deg_counts, file = file.path(sum_dir,"DE_summary.csv"), row.names = FALSE)

## simple html report (rmarkdown) -------------------------------------
if (requireNamespace("rmarkdown", quietly = TRUE)){
  rmd <- file.path(sum_dir,"analysis_report.Rmd")
  cat("
---
title: \"RNA-seq analysis (bulk, human)\"
output: html_document
toc: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(readr); library(knitr)
