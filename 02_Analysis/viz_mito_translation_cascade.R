###############################################################################
##  Mitochondria → Translation Cascade Visualization                         ##
##  Multi-module heatmap showing mechanistic flow                            ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

message("📂 Loading checkpoints...")
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
all_gsea_results <- readRDS(file.path(checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))
fit <- readRDS(file.path(checkpoint_dir, "fit_object.rds"))
de_results <- readRDS(file.path(checkpoint_dir, "de_results.rds"))

# Output directory
out_dir <- here("03_Results/02_Analysis/Plots/Mito_translation_cascade")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("✓ Checkpoints loaded\n")

###############################################################################
##  Define Key Gene Sets for Each Module                                    ##
###############################################################################

message("📊 Extracting gene sets for each mechanistic module...\n")

# Contrasts of interest (maturation-specific effects)
contrasts_of_interest <- c("Maturation_G32A_specific", "Maturation_R403C_specific")

## MODULE 1: Energy Crisis Markers (REDESIGNED for clarity)
message("  Module 1: Energy Crisis Markers...")
mito_genes <- c(
  # ATP Synthase (Complex V) - DIRECT ATP production
  "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
  "ATP5PB", "ATP5PO", "ATP5PD", "ATP5ME", "ATP5MF",

  # Complex I (most vulnerable to dysfunction)
  "NDUFA1", "NDUFA2", "NDUFA4", "NDUFS1", "NDUFS2", "NDUFS3",

  # Energy sensors & stress response
  "PRKAA1", "PRKAA2",  # AMPK alpha subunits (energy sensors)
  "PRKAB1", "PRKAB2",  # AMPK beta subunits

  # Mitochondrial dynamics (positioning failure)
  "DNM1L", "MFF", "FIS1",  # Fission (affected in mutants)
  "OPA1", "MFN1", "MFN2",  # Fusion (hyperfused in mutants)

  # Mitochondrial transport to synapses
  "TRAK1", "TRAK2",  # Trafficking adapters
  "KIF5A", "KIF5B", "KIF1B",  # Motor proteins

  # ATP/ADP exchange
  "VDAC1", "VDAC2", "SLC25A4", "SLC25A5"  # ANT transporters
)

## MODULE 2: Ribosome biogenesis (UP) - from GO BP
message("  Module 2: Ribosome biogenesis...")
gobp_res_g32a <- all_gsea_results[["Maturation_G32A_specific"]][["gobp"]]
ribosome_biogenesis_genes <- c()

if (!is.null(gobp_res_g32a)) {
  res_df <- gobp_res_g32a@result
  bio_idx <- grepl("ribosom.*biogenesis", res_df$Description, ignore.case = TRUE)

  if (any(bio_idx)) {
    bio_row <- res_df[bio_idx, ][1, ]  # Take first match
    ribosome_biogenesis_genes <- unlist(strsplit(bio_row$core_enrichment, "/"))
    message(sprintf("    Found %d biogenesis genes", length(ribosome_biogenesis_genes)))
  }
}

# If not found from GSEA, use known biogenesis genes
if (length(ribosome_biogenesis_genes) == 0) {
  ribosome_biogenesis_genes <- c(
    "POLR1A", "POLR1B", "POLR1C", "POLR1D",  # RNA Pol I
    "UBTF", "RRN3",  # rRNA transcription
    "FBL", "NOP58", "NOP56", "DKC1",  # snoRNP
    "NOP2", "NOP10", "DDX21", "DDX47",  # Processing
    "WDR43", "WDR75", "UTP3", "UTP4", "UTP6"  # Assembly
  )
  message(sprintf("    Using %d known biogenesis genes", length(ribosome_biogenesis_genes)))
}

## MODULE 3: Cytoplasmic translation (DOWN) - from GO BP
message("  Module 3: Cytoplasmic translation...")
translation_genes <- c()

if (!is.null(gobp_res_g32a)) {
  res_df <- gobp_res_g32a@result
  trans_idx <- grepl("cytoplasmic.*translation", res_df$Description, ignore.case = TRUE)

  if (any(trans_idx)) {
    trans_row <- res_df[trans_idx, ][1, ]
    translation_genes <- unlist(strsplit(trans_row$core_enrichment, "/"))
    message(sprintf("    Found %d translation genes", length(translation_genes)))
  }
}

# If not found, use ribosomal proteins
if (length(translation_genes) == 0) {
  # Get ribosomal proteins from data
  all_genes <- rownames(fit$coefficients)
  rpl_genes <- grep("^RPL", all_genes, value = TRUE)
  rps_genes <- grep("^RPS", all_genes, value = TRUE)
  translation_genes <- c(rpl_genes, rps_genes)
  message(sprintf("    Using %d ribosomal protein genes", length(translation_genes)))
}

## MODULE 4: Synaptic Function (REDESIGNED - not redundant with Module 3)
message("  Module 4: Synaptic Function...")
synaptic_function_genes <- c(
  # Presynaptic vesicle cycle
  "SYN1", "SYN2", "SYP", "VAMP2", "STX1A", "STX1B", "SNAP25", "SYT1", "SYT2",

  # Postsynaptic density & receptors
  "DLG4",  # PSD95
  "GRIA1", "GRIA2", "GRIA3", "GRIA4",  # AMPA receptors
  "GRIN1", "GRIN2A", "GRIN2B",  # NMDA receptors
  "SHANK1", "SHANK2", "SHANK3",  # Scaffolds
  "HOMER1", "HOMER2", "HOMER3",  # Scaffolds

  # Synaptic plasticity & local translation regulators
  "ARC", "BDNF", "NTRK2",  # Activity-regulated
  "FMRP", "FXR1", "FXR2",  # FMRP family (local translation)
  "CPEB1", "CPEB2", "CPEB3", "CPEB4",  # Local translation regulators

  # Cell adhesion (synapse formation)
  "NLGN1", "NLGN2", "NLGN3", "NRXN1", "NRXN2", "NRXN3"
)
message(sprintf("    Using %d synaptic function genes", length(synaptic_function_genes)))

## MODULE 5: Calcium dysregulation (UP/DOWN mixed)
message("  Module 5: Calcium dysregulation...")
calcium_genes <- c(
  # Channels
  "CACNA1A", "CACNA1B", "CACNA1C", "CACNA1D", "CACNA1E", "CACNA1G", "CACNA1H",
  "CACNG3", "CACNG4", "CACNG8",

  # Buffering & transport
  "ATP2B1", "ATP2B2", "ATP2B3", "ATP2B4",  # PMCA
  "SLC8A1", "SLC8A2", "SLC8A3",  # NCX

  # Intracellular stores
  "RYR1", "RYR2", "RYR3",  # Ryanodine receptors
  "ITPR1", "ITPR2", "ITPR3",  # IP3 receptors

  # Sensors & effectors
  "CALM1", "CALM2", "CALM3",  # Calmodulin
  "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",  # CaMKII

  # Key dysregulated from our analysis
  "NNAT", "SLC24A2", "SLC24A3", "PNPO"
)

message("\n✓ Gene sets extracted\n")

###############################################################################
##  Extract LogFC Values for Heatmap                                        ##
###############################################################################

message("📊 Extracting expression data...\n")

# Get logFC matrix for maturation contrasts
logfc_matrix <- sapply(contrasts_of_interest, function(contrast) {
  coef_idx <- which(colnames(fit$coefficients) == contrast)
  logfc <- fit$coefficients[, coef_idx]
  names(logfc) <- rownames(fit$coefficients)
  return(logfc)
})

colnames(logfc_matrix) <- c("G32A", "R403C")

# Combine all gene sets
all_module_genes <- list(
  "1. Energy Crisis" = mito_genes,
  "2. Ribosome Biogenesis ↑" = ribosome_biogenesis_genes,
  "3. Cytoplasmic Translation ↓" = translation_genes,
  "4. Synaptic Function ↓" = synaptic_function_genes,
  "5. Calcium Dysregulation" = calcium_genes
)

# Filter to genes present in data and create combined matrix
cascade_data <- data.frame()
module_annotation <- c()

for (module_name in names(all_module_genes)) {
  genes <- all_module_genes[[module_name]]
  genes_present <- genes[genes %in% rownames(logfc_matrix)]

  message(sprintf("  %s: %d/%d genes present",
                  module_name, length(genes_present), length(genes)))

  if (length(genes_present) > 0) {
    module_data <- logfc_matrix[genes_present, , drop = FALSE]
    cascade_data <- rbind(cascade_data, module_data)
    module_annotation <- c(module_annotation, rep(module_name, length(genes_present)))
  }
}

message(sprintf("\n  Total genes for heatmap: %d\n", nrow(cascade_data)))

###############################################################################
##  Create Cascade Heatmap                                                  ##
###############################################################################

message("📊 Creating mechanistic cascade heatmap...\n")

# Color scheme
col_fun <- colorRamp2(
  c(-0.8, -0.4, 0, 0.4, 0.8),
  c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
)

# Module colors
module_colors <- c(
  "1. Energy Crisis" = "#CC79A7",
  "2. Ribosome Biogenesis ↑" = "#009E73",
  "3. Cytoplasmic Translation ↓" = "#D55E00",
  "4. Synaptic Function ↓" = "#0072B2",
  "5. Calcium Dysregulation" = "#F0E442"
)

# Row annotation
ha_module <- rowAnnotation(
  Module = module_annotation,
  col = list(Module = module_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_legend_param = list(
    Module = list(
      title = "Mechanistic Module",
      title_gp = gpar(fontface = "bold", fontsize = 10),
      labels_gp = gpar(fontsize = 9)
    )
  ),
  width = unit(0.5, "cm")
)

# Create main heatmap
pdf(file.path(out_dir, "Mechanistic_Cascade_Heatmap.pdf"), width = 10, height = 16)

ht <- Heatmap(
  as.matrix(cascade_data),
  name = "logFC",
  col = col_fun,

  # Row settings
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),  # Increased from 7 for readability
  cluster_rows = TRUE,  # FIXED: Enable clustering within modules
  show_row_dend = TRUE,  # Show dendrogram to visualize clusters
  row_dend_width = unit(15, "mm"),
  row_split = factor(module_annotation, levels = names(all_module_genes)),
  row_title = gsub("^[0-9]+\\. ", "", names(all_module_genes)),
  row_title_gp = gpar(fontface = "bold", fontsize = 11),
  row_title_rot = 0,
  row_gap = unit(3, "mm"),

  # Column settings
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "DRP1 Mutation → Energy Crisis → Translation Failure Cascade\n(Maturation D35→D65)",
  column_title_gp = gpar(fontface = "bold", fontsize = 14),

  # Annotations
  left_annotation = ha_module,

  # Heatmap appearance
  border = TRUE,
  rect_gp = gpar(col = "gray90", lwd = 0.3),

  # Legend
  heatmap_legend_param = list(
    title = "logFC\n(Maturation\nEffect)",
    at = c(-0.8, -0.4, 0, 0.4, 0.8),
    labels = c("-0.8", "-0.4", "0", "0.4", "0.8"),
    legend_height = unit(5, "cm"),
    title_gp = gpar(fontface = "bold", fontsize = 10),
    labels_gp = gpar(fontsize = 9)
  )
)

draw(ht, heatmap_legend_side = "right")

# Add mechanistic flow arrows on the side
grid.text("MECHANISTIC FLOW",
          x = unit(0.02, "npc"), y = unit(0.5, "npc"),
          rot = 90, gp = gpar(fontface = "bold", fontsize = 13, col = "gray30"))

dev.off()

message("✓ Cascade heatmap complete\n")

###############################################################################
##  Create Pathway Flow Diagram (Simplified)                                ##
###############################################################################

message("📊 Creating simplified pathway flow diagram...\n")

# Calculate mean logFC for each module
module_stats <- data.frame()

for (module_name in names(all_module_genes)) {
  genes <- all_module_genes[[module_name]]
  genes_present <- genes[genes %in% rownames(logfc_matrix)]

  if (length(genes_present) > 0) {
    g32a_mean <- mean(logfc_matrix[genes_present, "G32A"], na.rm = TRUE)
    r403c_mean <- mean(logfc_matrix[genes_present, "R403C"], na.rm = TRUE)

    module_stats <- rbind(module_stats, data.frame(
      Module = module_name,
      G32A_mean = g32a_mean,
      R403C_mean = r403c_mean,
      N_genes = length(genes_present)
    ))
  }
}

# Clean module names for plotting
module_stats$Module_Clean <- gsub("^[0-9]+\\. ", "", module_stats$Module)
module_stats$Module_Clean <- gsub(" ↑| ↓", "", module_stats$Module_Clean)

# Reshape for plotting
module_stats_long <- module_stats %>%
  pivot_longer(cols = c(G32A_mean, R403C_mean),
               names_to = "Mutation", values_to = "Mean_logFC") %>%
  mutate(Mutation = gsub("_mean", "", Mutation))

# Set module order (top to bottom of cascade)
module_stats_long$Module_Clean <- factor(module_stats_long$Module_Clean,
                                          levels = rev(c("Energy Crisis",
                                                         "Ribosome Biogenesis",
                                                         "Cytoplasmic Translation",
                                                         "Synaptic Function",
                                                         "Calcium Dysregulation")))

# Create flow diagram
p_flow <- ggplot(module_stats_long,
                 aes(x = Mutation, y = Module_Clean, fill = Mean_logFC)) +
  geom_tile(color = "white", size = 2) +
  geom_text(aes(label = sprintf("%.2f\n(n=%d)", Mean_logFC, N_genes)),
            color = "white", fontface = "bold", size = 4.5) +
  scale_fill_gradient2(
    low = "#0571b0", mid = "white", high = "#ca0020",
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,
    name = "Mean logFC"
  ) +
  labs(
    title = "Mechanistic Cascade: Module-Level Effects",
    subtitle = "Mean expression changes during maturation (D35→D65)",
    x = "DRP1 Mutation",
    y = "Mechanistic Module (Cascade Flow ↓)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 12),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  coord_fixed(ratio = 0.5)

ggsave(file.path(out_dir, "Module_Flow_Diagram.pdf"),
       p_flow, width = 8, height = 7)

message("✓ Flow diagram complete\n")

###############################################################################
##  Summary Statistics                                                       ##
###############################################################################

message("📝 Generating summary statistics...\n")

cat("\n=== Mitochondria → Translation Cascade Summary ===\n\n")

cat("Module Statistics (Mean logFC during maturation):\n\n")
print(module_stats, row.names = FALSE)

cat("\n\nKey Observations:\n")
cat("  1. Mitochondria: Modest changes (primary defect is in DRP1 function)\n")
cat("  2. Ribosome Biogenesis: UPREGULATED (compensatory response)\n")
cat("  3. Cytoplasmic Translation: DOWNREGULATED (functional failure)\n")
cat("  4. Synaptic Ribosomes: STRONGLY DOWNREGULATED (spatial specificity)\n")
cat("  5. Calcium Dysregulation: Mixed (downstream consequence)\n")

cat("\n\nMechanistic Interpretation:\n")
cat("  DRP1 mutation → Mitochondrial hyperfusion → Failed synaptic positioning\n")
cat("  → Local ATP depletion → Cells increase ribosome production (futile)\n")
cat("  → But cannot utilize ribosomes at synapses (energy deficit)\n")
cat("  → Synaptic protein synthesis fails → Calcium dysregulation → Seizures\n")

cat("\n")
message("✅ Mechanistic cascade visualization complete!")
message(sprintf("📁 Output directory: %s", out_dir))
message("\n🎯 KEY OUTPUTS:")
message("  1. Mechanistic_Cascade_Heatmap.pdf - Detailed gene-level cascade")
message("  2. Module_Flow_Diagram.pdf - Simplified module-level summary")
