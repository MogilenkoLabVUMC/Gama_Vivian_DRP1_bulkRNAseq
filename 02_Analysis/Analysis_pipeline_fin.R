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
  "scripts/DE/volcano_helpers.R"
)

for (h in gsea_helpers){
  source_if_present(config$helper_root, h)
}

## other single helpers
source_if_present("01_Scripts/R_scripts/read_count_matrix.R")
source_if_present("01_Scripts/GSEA_module/scripts/DE/plot_standard_volcano.R")
source_if_present("01_Scripts/GSEA_module/scripts/DE/plotPCA.R")

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
  Int_G32A           = (D65_G32A  - D35_G32A)  - (D65_Control - D35_Control),
  Int_R403C          = (D65_R403C - D35_R403C) - (D65_Control - D35_Control),
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
deg_dir      <- here::here(config$out_root, "DE_results")
volcano_dir  <- here::here(config$out_root, "Plots/Volcano")
dir.create(deg_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

## helper to make + save volcanoes -------------------------------------
make_all_volcanoes <- function(res_table, name){
  pdf(file.path(volcano_dir, paste0(name,"_standard.pdf")), 8, 7)
  print(create_standard_volcano(
          res_table, p_cutoff = config$p_cutoff, fc_cutoff = config$fc_cutoff,
          label_method = "sig", highlight_gene = config$calcium_genes,
          title = name))
  dev.off()
}

for (co in colnames(contrasts)){
  tbl <- topTable(fit, coef = co, number = Inf)
  write.csv(tbl, file = file.path(deg_dir, paste0(co,"_results.csv")))
  make_all_volcanoes(tbl, co)
}
# the are a couple of issues with volcano plots 
# the legend color dots are intersected with symbol "a"
# the labeled genes significant genes have same yellow color as its corresponding yellow significant color, which is fine, but sometimes makes them hard to read when yellow intersects yellow 
# the highlighted calcium genes have the color of their group, depending if they are significant or not, but most of them are not, intersecting with the background genes dots of same color, making them impossible to read 
# the colorscheme must be colorblind friendly

# the volcano plot should have two visual styles: 
# style 1: all three colors for various groups + highlighted genes color
# style 2: only highlighted genes color

# the volcano plots should have two options
# y axis = -log10(p)
# y axis = -log10(FDR)

# all combinations of styles with all y axis styles should be saved

# the combined vertical volcano plots in the old style should be preserved



## DEG counts ----------------------------------------------------------
deg_counts <- data.frame(
  Contrast = colnames(contrasts),
  Up   = colSums(de_results > 0),
  Down = colSums(de_results < 0))
deg_counts$Total <- deg_counts$Up + deg_counts$Down
write.csv(deg_counts, file = file.path(deg_dir,"DE_gene_counts.csv"))

## bar plot (first-page bug fixed: call grid.newpage()) ----------------
pdf(file.path(qc_dir,"DE_gene_counts.pdf"), 10, 6)
grid::grid.newpage()
ggplot(reshape2::melt(deg_counts[,c("Contrast","Up","Down")],
                      id = "Contrast"),
       aes(Contrast, value, fill = variable))+
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values = c(Up="#D95F02",Down="#1B9E77"))+
  theme_minimal(base_size = 12)+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  labs(y="gene count", fill = "")
dev.off()
# for some reason generate a two page pdf, with first page empty, second with the plot 
# need to have numbers of DE genes on top of the bars 
# the order should be 
# G32A_vs_Ctrl_D35
# R403C_vs_Ctrl_D35
# G32A_vs_Ctrl_D65
# R403C_vs_Ctrl_D65
# Time_Ctrl
# Time_G32A
# Time_R403C
# Int_G32A
# Int_R403C

# -------------------------------------------------------------------- #
# 6.  Generic MSigDB-based GSEA                                        #
# -------------------------------------------------------------------- #
gsea_root <- here::here(config$out_root,"Plots/GSEA")
dir.create(gsea_root, recursive = TRUE, showWarnings = FALSE)

## wrapper that chooses HS vs MM automatically -------------------------
run_gsea_hsmm <- function(tbl, contrast, species){
  run_gsea_analysis(
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
}

for (co in colnames(contrasts)){
  tbl <- topTable(fit, coef = co, number = Inf)
  run_gsea_hsmm(tbl, co, species = "Homo sapiens")  # db_species picked up as HS
}
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/custom_minimal_theme.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/custom_minimal_theme.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_barplot.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_barplot.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_heatmap.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_heatmap.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_processing/run_gsea.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_processing/run_gsea.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/DE/volcano_helpers.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/DE/volcano_helpers.R
# ▶  Hallmark
# Fetching MSigDB sets for species='Homo sapiens', db_species='MM', collection='H', subcollection=''
# Running clusterProfiler::GSEA...
# GSEA completed successfully.
# ▶  GO BP
# Fetching MSigDB sets for species='Homo sapiens', db_species='MM', collection='C5', subcollection='GO:BP'
# Running clusterProfiler::GSEA...
#...
# ▶  GO CC
# Fetching MSigDB sets for species='Homo sapiens', db_species='MM', collection='C5', subcollection='GO:CC'
# Running clusterProfiler::GSEA...
# GSEA completed successfully.
# ▶  KEGG
# Fetching MSigDB sets for species='Homo sapiens', db_species='MM', collection='C2', subcollection='CP:KEGG_MEDICUS'
# Running clusterProfiler::GSEA...
# GSEA completed successfully.
# ▶  Reactome
# Fetching MSigDB sets for species='Homo sapiens', db_species='MM', collection='C2', subcollection='CP:REACTOME'
# Running clusterProfiler::GSEA...
# GSEA completed successfully.
# ▶  WikiPath
# Fetching MSigDB sets for species='Homo sapiens', db_species='MM', collection='C2', subcollection='CP:WIKIPATHWAYS'
# Running clusterProfiler::GSEA...
# GSEA completed successfully.
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/custom_minimal_theme.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/custom_minimal_theme.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_barplot.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_barplot.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_heatmap.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_heatmap.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_processing/run_gsea.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/GSEA/GSEA_processing/run_gsea.R
# [DEBUG] sourcing 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/DE/volcano_helpers.R
# [run_gsea_analysis] helper not found → 01_Scripts/GSEA_module/R_GSEA_visualisations/scripts/DE/volcano_helpers.R
# ▶  Hallmark
# Fetching MSigDB sets for species='Homo sapiens', db_species='MM', collection='H', subcollection=''
# Running clusterProfiler::GSEA...
# GSEA completed successfully.
# ▶  GO BP
# Fetching MSigDB sets for species='Homo sapiens', db_species='MM', collection='C5', subcollection='GO:BP'
# Running clusterProfiler::GSEA...
# GSEA completed successfully.
# Error in if (abs(max.ES) > abs(min.ES)) { : 
#   missing value where TRUE/FALSE needed
# In addition: There were 31 warnings (use warnings() to see them)

# -------------------------------------------------------------------- #
# 7.  SynGO GSEA (new)                                                 #
# -------------------------------------------------------------------- #
source("/workspaces/GVDRP1/01_Scripts/R_scripts/run_syngo_gsea.R")

# syngo_gmt <- function(syngo_dir, namespace = "CC"){
#   require(readxl)
#   ann  <- readxl::read_xlsx(file.path(syngo_dir,"syngo_annotations.xlsx"))
#   ont  <- readxl::read_xlsx(file.path(syngo_dir,"syngo_ontologies.xlsx"))
#   ann  <- subset(ann, go_domain == namespace)
#   ont  <- unique(ont[,c("id","name")])
#   term2gene <- merge(ann[,c("go_id","hgnc_symbol")],
#                      ont, by.x="go_id", by.y="id")
#   colnames(term2gene)[1:2] <- c("gs_name","gene_symbol")
#   return(term2gene)
# }

# term2gene_syngo <- syngo_gmt(config$syngo_dir, config$syngo_ns)
# head(term2gene_syngo)

syngo_gmt <- function(syngo_dir, namespace = "CC") {
  requireNamespace("readxl", quietly = TRUE)

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

syngo_lists <- syngo_gmt(config$syngo_dir, config$syngo_ns)

for (co in colnames(contrasts)) {
  tbl <- topTable(fit, coef = co, number = Inf)
  run_syngo_gsea(tbl, co,
                 T2G = syngo_lists$T2G,
                 T2N = syngo_lists$T2N)
}

# run_syngo_gsea <- function(tbl, contrast, term2gene){
#   gene_vec <- sort(setNames(tbl$t, rownames(tbl)), decreasing = TRUE)
#   set.seed(123)
#   res <- GSEA(geneList = gene_vec,
#               TERM2GENE = term2gene,
#               pvalueCutoff = 1,
#               pAdjustMethod = "fdr",
#               eps = 0,
#               by = "fgsea",
#               nPermSimple = 100000)
#   out_dir <- file.path(gsea_root, contrast, "SynGO")
#   dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#   saveRDS(res, file = file.path(out_dir, "GSEA_SynGO_result.rds"))
#   # quick dotplot
#   pdf(file.path(out_dir, "SynGO_dotplot.pdf"), 7, 5)
#   print( gsea_dotplot(res, title = paste0(contrast," – SynGO"),
#                       pos_color = "#fc8d59", neg_color = "#91bfdb") )
#   dev.off()
#   invisible(res)
# }

# for (co in colnames(contrasts)){
#   tbl <- topTable(fit, coef = co, number = Inf)
#   run_syngo_gsea(tbl, co, term2gene_syngo)
# }
# r$> for (co in colnames(contrasts)){                                                                               
#       tbl <- topTable(fit, coef = co, number = Inf)                                                                                                                
#       run_syngo_gsea(tbl, co, term2gene_syngo)                                                             
#     }                                                                           

# using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).                                                                         
# preparing geneSet collections...                                                                  
# GSEA analysis...                                                                     
# leading edge analysis...                                                                     
# done...                                                                                                            
# using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).                                                                         
# preparing geneSet collections...                                                                  
# GSEA analysis...                                                                     
# leading edge analysis...                                                                     
# done...                                                                         
# using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).                                             
# preparing geneSet collections...                                                                  
# GSEA analysis...                                                                     
# leading edge analysis...                                                                     
# done...                                                                         
# using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).                                                                         
# preparing geneSet collections...                                                                  
# GSEA analysis...                                                                     
# leading edge analysis...                                                                     
# done...                                                                         
# using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).                                     
#....
# Warning messages:                                                                       
# 1: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.                                      
# ℹ Please use the `linewidth` argument instead.                                                                     
# This warning is displayed once every 8 hours.                                                                      
# Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.                               
# 2: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :                                           
#   for 'G32A_vs_Ctrl_D35 – SynGO' in 'mbcsToSbcs': - substituted for – (U+2013)                                     
# 3: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :                                           
#   for 'R403C_vs_Ctrl_D35 – SynGO' in 'mbcsToSbcs': - substituted for – (U+2013)                                    
# 4: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :                                           
#   for 'G32A_vs_Ctrl_D65 – SynGO' in 'mbcsToSbcs': - substituted for – (U+2013)                                     
# 5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :                                           
#   for 'R403C_vs_Ctrl_D65 – SynGO' in 'mbcsToSbcs': - substituted for – (U+2013)                                    
# 6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :                                           
#   for 'Time_Ctrl – SynGO' in 'mbcsToSbcs': - substituted for – (U+2013)                                            
# 7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :                                           
#   for 'Time_G32A – SynGO' in 'mbcsToSbcs': - substituted for – (U+2013)                                            
# 8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :                                           
#   for 'Time_R403C – SynGO' in 'mbcsToSbcs': - substituted for – (U+2013)                                           
# 9: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :                                           
#   for 'Int_G32A – SynGO' in 'mbcsToSbcs': - substituted for – (U+2013)                                             
# 10: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :                                          
#   for 'Int_R403C – SynGO' in 'mbcsToSbcs': - substituted for – (U+2013)                                            
                                                                                

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
# invalid or corrupted pdf file saved


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
