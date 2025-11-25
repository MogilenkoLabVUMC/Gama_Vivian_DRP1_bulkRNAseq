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
  helper_root   = "01_Scripts/RNAseq-toolkit",   # <<-- git submodule
  ## analysis parameters
  p_cutoff      = 0.05,
  fc_cutoff     = 2,
  calcium_genes = c(
      "NNAT","CACNG3","CACNA1S","ATP2A1",  # NNAT for Neuronatin gene
      "RYR1","MYLK3","VDR","STIM1","STIM2",
      "ORAI1_1","CALB1","CALR","PNPO"),  # ORAI1_1 (alternative symbol), PNPO added per collaborator request
  # Note: CACNA1C, CASR, ORAI1 not present in dataset. ORAI1_1 is the alternative symbol for ORAI1.
  ## SynGO data (relative to repo root)
  syngo_dir     = "00_Data/SynGO_bulk_20231201",
  syngo_ns      = "CC",                       # GO cellular-component ontology
  ## MitoCarta data (relative to repo root)
  mitocarta_file = "00_Data/MitoCarta_3.0/MitoPathways3.0.gmx",
  ## checkpoint caching
  force_recompute = FALSE   # Set TRUE to ignore cached checkpoints and recompute all
)

## create checkpoint directory for saving expensive computations -------
checkpoint_dir <- here::here(config$out_root, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

## checkpoint caching helper -------------------------------------------
#' Load cached result or compute if needed
#' @param checkpoint_file Path to checkpoint RDS file
#' @param compute_fn Function that computes the result
#' @param force_recompute Ignore cache and recompute
#' @param description Human-readable description for messages
load_or_compute <- function(checkpoint_file, compute_fn,
                            force_recompute = FALSE,
                            description = "computation") {
  if (!force_recompute && file.exists(checkpoint_file)) {
    message("📦 Loading cached ", description, " from ", basename(checkpoint_file))
    return(readRDS(checkpoint_file))
  }
  message("🔬 Computing ", description, "...")
  result <- compute_fn()
  message("💾 Saving ", description, " to ", basename(checkpoint_file))
  saveRDS(result, checkpoint_file)
  return(result)
}

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
  "scripts/GSEA/GSEA_plotting/format_pathway_names.R",
  "scripts/GSEA/GSEA_plotting/gsea_dotplot.R",
  "scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R",
  "scripts/GSEA/GSE
  A_plotting/gsea_barplot.R",
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
source_if_present(config$helper_root, "scripts/DE/plot_standard_volcano.R")
source_if_present(config$helper_root, "scripts/DE/create_fc_b_plot.R")
source_if_present(config$helper_root, "scripts/DE/create_MD_plot.R")
source_if_present(config$helper_root, "scripts/DE/plotPCA.R")
source_if_present("01_Scripts/R_scripts/generate_vertical_volcanos.R")

## GSEA-specific helpers (SynGO & MitoCarta)
source_if_present("01_Scripts/R_scripts/syngo_running_sum_plot.R")  # Specialized running sum plot
source_if_present("01_Scripts/R_scripts/run_syngo_gsea.R")
source_if_present("01_Scripts/R_scripts/parse_mitocarta_gmx.R")
source_if_present("01_Scripts/R_scripts/run_mitocarta_gsea.R")
source_if_present(config$helper_root, "scripts/GSEA/GSEA_plotting/plot_all_gsea_results.R")

## Load shared utilities (DRY helpers) - now in toolkit for reusability
source_if_present(config$helper_root, "scripts/utils_plotting.R")

# -------------------------------------------------------------------- #
# 2.  Read & pre-process data                                          #
# -------------------------------------------------------------------- #
DGE <- process_rnaseq_data(
  here::here(config$counts_file),
  here::here(config$metadata_file),
  annotate = FALSE
)

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
  # Question: How does maturation affect each genotype?I
  ## maturation
  Time_Ctrl          = D65_Control - D35_Control,
  Time_G32A          = D65_G32A - D35_G32A,
  Time_R403C         = D65_R403C - D35_R403C,

    # 3. Interaction effects (difference-in-difference)
  # Question: Do the mutations alter the normal maturation trajectory?
  ## interaction
  Maturation_G32A_specific = (D65_G32A - D35_G32A)  - (D65_Control - D35_Control),
  Maturation_R403C_specific = (D65_R403C - D35_R403C) - (D65_Control - D35_Control),
  levels = design)

# -------------------------------------------------------------------- #
# 3.  DEG modelling (voomLmFit preferred)                              #
# -------------------------------------------------------------------- #
## Use checkpoint caching to avoid re-running expensive voomLmFit + limma

model_objs <- load_or_compute(
  checkpoint_file = file.path(checkpoint_dir, "model_objects.rds"),
  force_recompute = config$force_recompute,
  description = "voomLmFit + limma differential expression model",
  compute_fn = function() {
    ## keep.EList=TRUE ➜ the returned object contains $EList with log-CPM + weights
    fit <- edgeR::voomLmFit(
             DGE,
             design,
             sample.weights = FALSE,   # turn off
             keep.EList     = TRUE)

    ## unpack the voom log-CPM for QC plots
    v <- fit$EList                    # same structure as `voom()` output

    ## add contrasts & empirical Bayes
    fit <- contrasts.fit(fit, contrasts) |>
           eBayes(robust = TRUE)

    de_results <- decideTests(fit, p.value = config$p_cutoff)

    ## return all objects as a list
    list(fit = fit, v = v, de_results = de_results,
         DGE = DGE, contrasts = contrasts)
  }
)

## unpack loaded/computed objects
fit        <- model_objs$fit
v          <- model_objs$v
de_results <- model_objs$de_results
DGE        <- model_objs$DGE
contrasts  <- model_objs$contrasts

# ------------------------------------------------------------
# Process all contrasts ONCE (DRY principle - eliminates duplicate topTable calls)
# ------------------------------------------------------------
message("\n🔬 Processing all contrasts...")
contrast_tables <- process_all_contrasts(fit, contrasts, number = Inf, sort_by = "t")

# ------------------------------------------------------------
# Save DEG results to CSV
# ------------------------------------------------------------
deg_dir <- here::here(config$out_root, "DE_results")
save_contrast_csvs(contrast_tables, deg_dir, suffix = "_DE_results.csv")

# ------------------------------------------------------------
# CHECKPOINT: Save contrast_tables for downstream modules
# ------------------------------------------------------------
saveRDS(contrast_tables, file.path(checkpoint_dir, "contrast_tables.rds"))
message("💾 Saved contrast_tables checkpoint for downstream analysis")

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

# ------------------------------------------------------------
# CHECKPOINT: Save QC variables for downstream visualization
# ------------------------------------------------------------
saveRDS(list(ordered_samples = ordered_samples,
             annot = annot,
             ann_colors = ann_colors,
             logCPM = logCPM),
        file.path(checkpoint_dir, "qc_variables.rds"))

# -------------------------------------------------------------------- #
# 5.  Volcanoes & DEG numbers                                          #
# -------------------------------------------------------------------- #
# 5.  Volcano plots (horizontal standard - p-value and FDR modes)     #
# -------------------------------------------------------------------- #
# NOTE: Using pre-computed contrast_tables from earlier (DRY principle)

message("\n🌋 Generating volcano plots...")

## Helper function to generate volcano set (eliminates code duplication)
generate_volcano_set <- function(mode = c("p", "fdr")) {
  mode <- match.arg(mode)

  volcano_dir <- here::here(config$out_root, "Plots/Volcano", mode)
  ensure_dir(volcano_dir)

  # Mode-specific parameters
  params <- list(
    p = list(
      p_cutoff = config$p_cutoff,
      decision_by = "p",
      label_method = "top",
      title_suffix = sprintf("(Threshold: p ≤ %.2f)", config$p_cutoff),
      subtitle = "Highlighting by unadjusted p-value"
    ),
    fdr = list(
      p_cutoff = 0.1,
      decision_by = "fdr",
      label_method = "top",
      title_suffix = "(Threshold: FDR ≤ 0.1)",
      subtitle = "Highlighting by FDR (adj.P.Val), displayed on p-value scale for resolution"
    )
  )[[mode]]

  message("  Mode: ", toupper(mode))

  for (co in names(contrast_tables)) {
    plot <- create_standard_volcano(
      contrast_tables[[co]],
      p_cutoff = params$p_cutoff,
      fc_cutoff = config$fc_cutoff,
      decision_by = params$decision_by,
      label_method = params$label_method,
      highlight_gene = config$calcium_genes,
      x_breaks = 2,
      title = paste(co, params$title_suffix),
      subtitle = params$subtitle
    )

    save_plot(plot,
              file.path(volcano_dir, paste0(co, "_standard.pdf")),
              width = 8, height = 7)
  }

  message("  ✓ Saved ", length(contrast_tables), " ", mode, " volcano plots")
}

# Generate both p-value and FDR volcano sets
generate_volcano_set("p")
generate_volcano_set("fdr")

## Vertical volcanoes (uses pre-computed contrast_tables)
message("\n🌋 Generating vertical volcano plots...")
generate_vertical_volcano_sets(contrast_tables, config, highlight_calcium = TRUE)
message("✓ Vertical volcano plots complete")

## ---- MD plots (uses pre-computed contrast_tables) ---------------------
message("\n📊 Generating MD plots...")
md_dir <- here::here(config$out_root, "Plots/Volcano/MD")
ensure_dir(md_dir)

for (co in names(contrast_tables)) {
  plot <- create_MD_plot(
    fit = fit,
    coef = co,
    de_results = contrast_tables[[co]],
    fc_cutoff = 1,
    fdr_cutoff = 0.05,
    top_n = 5,
    highlight_gene = config$calcium_genes,
    label_method = "top",
    title = paste("MD plot:", co),
    show_grid = FALSE
  )

  save_plot(plot, file.path(md_dir, paste0(co, "_MDplot.pdf")),
            width = 7, height = 6)
}
message("✓ MD plots complete")

## ---- FC vs B plots (uses pre-computed contrast_tables) -----------------
message("\n📊 Generating FC-vs-B plots...")
fc_b_dir <- here::here(config$out_root, "Plots/Volcano/FC-B")
ensure_dir(fc_b_dir)

for (co in names(contrast_tables)) {
  plot <- create_B_FC_plot(
    contrast_tables[[co]],
    fc_cutoff = config$fc_cutoff,
    B_cutoff = 0,
    top_n = 5,
    highlight_gene = config$calcium_genes,
    label_method = "top",
    title = paste("log2FC vs B:", co),
    show_grid = FALSE
  )

  save_plot(plot, file.path(fc_b_dir, paste0(co, "_FC_vs_B.pdf")),
            width = 7, height = 6)
}
message("✓ FC-vs-B plots complete")


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
      padj_cutoff  = 1, # to save all pathways for downstream viz
      output_dir   = file.path(gsea_root, contrast),
      sample_annotation = DGE$samples[,c("genotype","days"),drop=FALSE],
      sample_order = ordered_samples,
      helper_root  = config$helper_root,
      save_plots   = TRUE)
  
  return(results)  # Explicitly return the results
}

# Use checkpoint caching for expensive GSEA computation
all_gsea_results <- load_or_compute(
  checkpoint_file = file.path(checkpoint_dir, "all_gsea_results.rds"),
  force_recompute = config$force_recompute,
  description = "Generic MSigDB GSEA for all contrasts",
  compute_fn = function() {
    # Create a list to store all GSEA results
    results <- list()

    # Run GSEA for each contrast (uses pre-computed contrast_tables)
    message("Running GSEA for all contrasts...")
    for (co in names(contrast_tables)) {
      message("==== Working on contrast: ", co, " ====")

      safe_res <- tryCatch({
        run_gsea_hsmm(contrast_tables[[co]], co, species = "Homo sapiens")
      }, error = function(e) {
        message("❗ ERROR in GSEA for contrast ", co, ": ", e$message)
        message("👀 Calling traceback():")
        traceback()
        NULL
      })

      results[[co]] <- safe_res
    }

    return(results)
  }
)

# Generate additional GSEA plots for each contrast
message("\n📊 Generating additional GSEA plots...")
for (co in names(all_gsea_results)) {
  this_res_list <- all_gsea_results[[co]]
  if (is.null(this_res_list) || !length(this_res_list)) {
    message("  Skipping ", co, " (no results)")
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
message("✓ GSEA analysis complete")

# -------------------------------------------------------------------- #
# 7.  SynGO GSEA                                                       #
# -------------------------------------------------------------------- #
source_if_present("01_Scripts/R_scripts/run_syngo_gsea.R")

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
syngo_lists <- syngo_gmt(here::here(config$syngo_dir), config$syngo_ns)

# Use checkpoint caching for SynGO GSEA
syngo_gsea_results <- load_or_compute(
  checkpoint_file = file.path(checkpoint_dir, "syngo_gsea_results.rds"),
  force_recompute = config$force_recompute,
  description = "SynGO GSEA for all contrasts",
  compute_fn = function() {
    # Create a list to store all GSEA results
    results <- list()

    # Run SynGO GSEA for each contrast (uses pre-computed contrast_tables)
    message("Running SynGO GSEA for all contrasts...")
    for (co in names(contrast_tables)) {
      message("==== Working on contrast: ", co, " ====")

      # Close any lingering graphic devices before starting new contrast
      close_all_devices()

      # Capture the result
      results[[co]] <- run_syngo_gsea(
        contrast_tables[[co]],
        co,
        T2G = syngo_lists$T2G,
        T2N = syngo_lists$T2N,
        sample_annotation = annot
      )

      # Ensure all devices are closed after each iteration
      close_all_devices()
    }

    return(results)
  }
)
message("✓ SynGO GSEA complete")

# -------------------------------------------------------------------- #
# 8.  MitoCarta GSEA - MOVED TO SEPARATE SCRIPT                        #
# -------------------------------------------------------------------- #
# NOTE: MitoCarta GSEA analysis has been moved to 2.add_MitoCarta.R
#       Run that script separately after this main pipeline completes.
#       This separation improves modularity and debugging.
message("\n📝 MitoCarta GSEA analysis moved to 2.add_MitoCarta.R")
message("   Run that script separately after this pipeline completes.")

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
## 4  define/load gene intersections for enrichment                   ##
## ------------------------------------------------------------------ ##
## baseline9: genes DE in BOTH G32A_vs_Ctrl AND R403C_vs_Ctrl at D35
## mat38:     genes DE in BOTH Maturation_G32A_specific AND Maturation_R403C_specific

# Try loading existing gene lists from 1stRun
baseline9_file <- here::here(config$out_root, "1stRun/SynGO_intersections/baseline9_genes.txt")
mat38_file     <- here::here(config$out_root, "1stRun/SynGO_intersections/mat38_genes.txt")

# Compute intersections from current DE results for verification
baseline9_computed <- intersect(
  rownames(fit)[de_results[, "G32A_vs_Ctrl_D35"] != 0],
  rownames(fit)[de_results[, "R403C_vs_Ctrl_D35"] != 0]
)

mat38_computed <- intersect(
  rownames(fit)[de_results[, "Maturation_G32A_specific"] != 0],
  rownames(fit)[de_results[, "Maturation_R403C_specific"] != 0]
)

# Load existing files if they exist, otherwise use computed
if (file.exists(baseline9_file)) {
  baseline9_loaded <- readLines(baseline9_file)
  message("📂 Loaded baseline9 from 1stRun: ", length(baseline9_loaded), " genes")
  message("🔬 Computed from current analysis: ", length(baseline9_computed), " genes")
  baseline9 <- baseline9_loaded  # Use loaded version
} else {
  message("⚠️  baseline9_genes.txt not found, using computed version")
  baseline9 <- baseline9_computed
}

if (file.exists(mat38_file)) {
  mat38_loaded <- readLines(mat38_file)
  message("📂 Loaded mat38 from 1stRun: ", length(mat38_loaded), " genes")
  message("🔬 Computed from current analysis: ", length(mat38_computed), " genes")
  mat38 <- mat38_loaded  # Use loaded version
} else {
  message("⚠️  mat38_genes.txt not found, using computed version")
  mat38 <- mat38_computed
}

# Save gene lists to checkpoints
saveRDS(list(baseline9 = baseline9,
             mat38 = mat38,
             baseline9_computed = baseline9_computed,
             mat38_computed = mat38_computed),
        file.path(checkpoint_dir, "gene_intersections.rds"))

## ------------------------------------------------------------------ ##
## 5  run enrichment                                                  ##
## ------------------------------------------------------------------ ##
universe_all <- rownames(fit)          # every gene tested in limma

## path to the original SynGO gene table
syngo_gene_file <- here::here(config$syngo_dir, "syngo_genes.xlsx")

enr_base9 <- synGO_enrich_symbol(baseline9, "baseline9",
                                 universe_all,
                                 syngo_lists$T2G, syngo_lists$T2N,
                                 syngo_gene_file)

enr_mat38 <- synGO_enrich_symbol(mat38,     "mat38",
                                 universe_all,
                                 syngo_lists$T2G, syngo_lists$T2N,
                                 syngo_gene_file)


# -------------------------------------------------------------------- #
# 9.  Calcium genes focus                                              #
# -------------------------------------------------------------------- #
message("🔬 Running calcium gene analysis (separate script)...")
source(here::here("02_Analysis/viz_calcium_genes.R"))


# -------------------------------------------------------------------- #
# 10.  Summary tables                                                  #
# -------------------------------------------------------------------- #
sum_dir <- here::here(config$out_root,"Summary")
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(deg_counts, file = file.path(sum_dir,"DE_summary.csv"), row.names = FALSE)
