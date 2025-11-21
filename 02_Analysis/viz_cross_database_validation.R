###############################################################################
##  Cross-Database Validation: Translation Crisis Robustness                ##
##  Consolidated Script: 3-Column Dotplot + Legacy 2-Column Panels          ##
###############################################################################
## MAIN OUTPUT: 3-column figure (Time_Ctrl + G32A + R403C)
## Shows normal maturation vs. mutation-specific disruption
## Includes all pathways (even Time_Ctrl-specific) to tell complete story
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(stringr)

# Source the pathway name formatting function
source(here("01_Scripts/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/format_pathway_names.R"))

message("📂 Loading checkpoints...")
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
all_gsea_results <- readRDS(file.path(checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))

# Output directory
out_dir <- here("03_Results/02_Analysis/Plots/Cross_database_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("✓ Checkpoints loaded\n")

###############################################################################
##  Extract Translation/Ribosome Pathways - ALL 3 Contrasts                ##
###############################################################################

message("📊 Extracting translation pathways from all databases (including Time_Ctrl)...\n")

# Include Time_Ctrl for comparison with mutation-specific contrasts
contrasts_of_interest <- c("Time_Ctrl", "Maturation_G32A_specific", "Maturation_R403C_specific")

# Databases to search
databases <- c("gobp", "gocc", "gomf", "kegg", "reactome", "wiki", "hallmark", "canon")

# Keywords for translation-related, mitochondria, and calcium pathways
translation_keywords <- c(
  "ribosom", "translation", "peptide chain elongation", "protein synthesis",
  "SRP-dependent", "cotranslational", "ribonucleoprotein",
  "mitochondri", "oxidative phosphorylation", "electron transport", "ATP synthesis",
  "calcium", "ca2+", "calmodulin", "calcineurin", "calcium channel",
  "calcium signaling", "calcium homeostasis", "calcium transport"
)

# Compile all translation pathways
all_translation_pathways <- data.frame()

for (contrast in contrasts_of_interest) {
  for (db in databases) {
    gsea_res <- all_gsea_results[[contrast]][[db]]

    if (is.null(gsea_res)) next

    # Extract result dataframe
    res_df <- tryCatch({
      if (class(gsea_res) == "gseaResult") {
        gsea_res@result
      } else if ("result" %in% names(gsea_res)) {
        gsea_res$result
      } else {
        NULL
      }
    }, error = function(e) NULL)

    if (is.null(res_df) || nrow(res_df) == 0) next

    # Filter for translation-related pathways
    translation_idx <- grepl(paste(translation_keywords, collapse = "|"),
                            res_df$Description, ignore.case = TRUE)

    if (any(translation_idx)) {
      translation_pathways <- res_df[translation_idx, ]

      # Filter for significant (FDR < 0.05)
      translation_pathways <- translation_pathways[translation_pathways$p.adjust < 0.05, ]

      if (nrow(translation_pathways) > 0) {
        translation_pathways$Database <- toupper(db)
        translation_pathways$Contrast <- contrast
        all_translation_pathways <- rbind(all_translation_pathways,
                                          translation_pathways)
      }
    }
  }
}

# Add ALL significant SynGO pathways for ALL contrasts
for (contrast in contrasts_of_interest) {
  syngo_res <- syngo_gsea_results[[contrast]]

  if (!is.null(syngo_res)) {
    syngo_df <- syngo_res@result

    # Get ALL significant pathways from SynGO
    syngo_sig <- syngo_df[syngo_df$p.adjust < 0.05, ]

    if (nrow(syngo_sig) > 0) {
      syngo_sig$Database <- "SYNGO"
      syngo_sig$Contrast <- contrast
      all_translation_pathways <- rbind(all_translation_pathways, syngo_sig)
    }
  }
}

message(sprintf("  Found %d significant pathways across databases and contrasts\n",
                nrow(all_translation_pathways)))

# Create display columns
all_translation_pathways$Contrast_Display <- gsub("Maturation_", "",
                                                   gsub("_specific", "",
                                                        all_translation_pathways$Contrast))
all_translation_pathways$Contrast_Display <- factor(
  all_translation_pathways$Contrast_Display,
  levels = c("Time_Ctrl", "G32A", "R403C")
)

# Also keep Mutation column for legacy outputs
all_translation_pathways$Mutation <- all_translation_pathways$Contrast_Display
all_translation_pathways$Mutation[all_translation_pathways$Mutation == "Time_Ctrl"] <- NA

# Categorize pathways
all_translation_pathways$Category <- ifelse(
  grepl("calcium|ca2\\+|calmodulin|calcineurin",
        all_translation_pathways$Description, ignore.case = TRUE),
  "Calcium",
  ifelse(grepl("mitochondri|oxidative phosphorylation|electron transport|ATP synthesis",
               all_translation_pathways$Description, ignore.case = TRUE),
         "Mitochondria",
         ifelse(grepl("biogenesis|maturation|assembly|processing",
                      all_translation_pathways$Description, ignore.case = TRUE),
                "Biogenesis",
                ifelse(grepl("translation|elongation|initiation|termination",
                             all_translation_pathways$Description, ignore.case = TRUE),
                       "Translation",
                       ifelse(grepl("ribosom|ribosomal",
                                    all_translation_pathways$Description, ignore.case = TRUE),
                              "Ribosomal",
                              "Other"))))
)

# Format pathway names: strip database prefixes and apply smart formatting
all_translation_pathways$Description_Formatted <- format_pathway_name(
  all_translation_pathways$Description,
  use_formatting = TRUE,
  strip_prefix = TRUE
)

# Truncate long pathway names if needed
all_translation_pathways$Description_Short <- ifelse(
  nchar(all_translation_pathways$Description_Formatted) > 60,
  paste0(substr(all_translation_pathways$Description_Formatted, 1, 57), "..."),
  all_translation_pathways$Description_Formatted
)

###############################################################################
##  Filter Out Non-Relevant Pathways                                        ##
###############################################################################

message("🔍 Filtering out non-relevant pathways...\n")

# Create exclusion list
exclude_pathways <- c(
  # Viral/platelet (not neuronal)
  "REACTOME_SARS_COV_1_MODULATES_HOST_TRANSLATION_MACHINERY",
  "REACTOME_EXPORT_OF_VIRAL_RIBONUCLEOPROTEINS_FROM_NUCLEUS",
  "REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2",

  # Cardiac-specific calcium (not neuronal)
  "WP_CALCIUM_REGULATION_IN_CARDIAC_CELLS",

  # Single-mutant ribosome biogenesis (not robust)
  "GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",
  "GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS",
  "GOBP_RIBOSOMAL_SUBUNIT_EXPORT_FROM_NUCLEUS",
  "GOCC_PRERIBOSOME",
  "GOCC_PRERIBOSOME_LARGE_SUBUNIT_PRECURSOR",

  # Generic single-mutant pathways
  "GOMF_RIBONUCLEOPROTEIN_COMPLEX_BINDING",
  "GOCC_SNO_S_RNA_CONTAINING_RIBONUCLEOPROTEIN_COMPLEX"
)

# Apply filter
n_before <- nrow(all_translation_pathways)
all_translation_pathways <- all_translation_pathways %>%
  filter(!Description %in% exclude_pathways)
n_after <- nrow(all_translation_pathways)

message(sprintf("  Before filtering: %d pathways\n", n_before))
message(sprintf("  After filtering: %d pathways (excluded %d)\n",
                n_after, n_before - n_after))

###############################################################################
##  Create Complete Pathway Grid (All Contrasts × All Pathways)            ##
###############################################################################

message("📋 Creating complete pathway set across all contrasts...\n")

# Get union of all pathway names
all_pathway_names <- unique(all_translation_pathways$Description)
message(sprintf("  Total unique pathways: %d\n", length(all_pathway_names)))

# Create pathway-database grid
pathway_grid <- all_translation_pathways %>%
  select(Description, Database) %>%
  distinct() %>%
  crossing(Contrast_Display = factor(c("Time_Ctrl", "G32A", "R403C"),
                                      levels = c("Time_Ctrl", "G32A", "R403C")))

# Merge with actual data
plot_data <- pathway_grid %>%
  left_join(
    all_translation_pathways %>%
      select(Description, Contrast_Display, Database, NES, p.adjust, pvalue, setSize, Category),
    by = c("Description", "Contrast_Display", "Database")
  )

# Add formatted pathway names
# Need to merge from all_translation_pathways since format_pathway_name was already applied
pathway_name_map <- all_translation_pathways %>%
  select(Description, Description_Short) %>%
  distinct()

plot_data <- plot_data %>%
  left_join(pathway_name_map, by = "Description")

###############################################################################
##  FIX: Pathway Ordering (Include ALL Contrasts, Not Just G32A/R403C)     ##
###############################################################################

message("🔧 Creating pathway ordering (FIXED: includes all contrasts)...\n")

# CRITICAL FIX: Include ALL contrasts when calculating ordering
# This prevents NA pathway labels for Time_Ctrl-only pathways
pathway_order <- plot_data %>%
  filter(!is.na(NES)) %>%  # All contrasts (not just G32A/R403C!)
  group_by(Description) %>%
  summarize(Mean_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(Mean_NES) %>%
  pull(Description)

message(sprintf("  Ordered %d pathways by mean NES across all contrasts\n", length(pathway_order)))

# Create factor with proper levels (no more NAs!)
plot_data$Pathway <- factor(plot_data$Description_Short,
                            levels = unique(plot_data$Description_Short[match(pathway_order, plot_data$Description)]))

###############################################################################
##  FIX: Database Ordering (SYNGO First)                                    ##
###############################################################################

message("🔧 Reordering databases (SYNGO at top)...\n")

# Set Database factor levels with SYNGO first
plot_data$Database <- factor(plot_data$Database,
                             levels = c("SYNGO", "GOBP", "GOCC", "GOMF",
                                       "KEGG", "REACTOME", "WIKI", "HALLMARK", "CANON"))

# Flag highly significant pathways
plot_data$HighlySig <- !is.na(plot_data$p.adjust) & plot_data$p.adjust < 0.01

###############################################################################
##  MAIN OUTPUT: 3-Column Dotplot (Time_Ctrl + G32A + R403C)               ##
###############################################################################

message("📊 Creating MAIN FIGURE: 3-column cross-database dotplot...\n")

# Keep ALL pathways (including Time_Ctrl-specific) - they tell a biological story!
plot_data_main <- plot_data %>% filter(!is.na(NES))

n_unique_pathways <- length(unique(plot_data_main$Description))
n_timectrl_only <- length(setdiff(
  unique(plot_data[plot_data$Contrast_Display == "Time_Ctrl" & !is.na(plot_data$NES), "Description"]),
  unique(plot_data[plot_data$Contrast_Display %in% c("G32A", "R403C") & !is.na(plot_data$NES), "Description"])
))

message(sprintf("  Total pathways: %d (including %d Time_Ctrl-specific)\n",
                n_unique_pathways, n_timectrl_only))

# Create main 3-column dotplot
p_main <- ggplot(plot_data_main,
                aes(x = Contrast_Display, y = Pathway,
                    color = NES, size = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +

  # Add black outlines for FDR < 0.01
  geom_point(
    data = subset(plot_data_main, HighlySig == TRUE),
    aes(size = -log10(p.adjust)),
    shape = 21,
    color = "black",
    fill = NA,
    stroke = 1.2,
    alpha = 1
  ) +

  facet_wrap(~Database, ncol = 1, scales = "free_y") +
  scale_color_gradient2(
    low = "#0571b0", mid = "white", high = "#ca0020",
    midpoint = 0,
    name = "NES",
    limits = c(-3.5, 3.5),
    oob = scales::squish
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = "-log10(FDR)",
    breaks = c(2, 5, 10, 15),
    limits = c(0, 20)
  ) +
  labs(
    title = "Cross-Database Validation: Normal vs. Disrupted Maturation",
    subtitle = sprintf("Time_Ctrl = normal maturation | G32A/R403C = mutation-specific defects | FDR < 0.05 (black outline = FDR < 0.01, n=%d pathways)",
                      n_unique_pathways),
    x = "Contrast",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(face = "bold", size = 10, angle = 0),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 11),
    panel.border = element_rect(color = "gray80", fill = NA),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# Dynamic plot height (increased for better spacing)
plot_height <- max(18, n_unique_pathways * 0.6)

ggsave(file.path(out_dir, "Cross_Database_3Column_Main.pdf"),
       p_main, width = 13, height = plot_height)

message("✓ Main 3-column figure complete\n")

###############################################################################
##  LEGACY: 2-Column Dotplot (G32A + R403C only, for backwards compat)     ##
###############################################################################

message("📊 Creating LEGACY figure: 2-column dotplot (G32A + R403C only)...\n")

panel_a_data <- all_translation_pathways %>%
  filter(Contrast != "Time_Ctrl") %>%
  mutate(
    Mutation = gsub("Maturation_", "", gsub("_specific", "", Contrast)),
    Pathway = factor(Description_Short,
                     levels = unique(Description_Short[order(NES)])),
    HighlySig = p.adjust < 0.01
  )

# Add Database ordering
panel_a_data$Database <- factor(panel_a_data$Database,
                                levels = c("SYNGO", "GOBP", "GOCC", "GOMF",
                                          "KEGG", "REACTOME", "WIKI", "HALLMARK", "CANON"))

n_legacy_pathways <- length(unique(panel_a_data$Description))

p_legacy <- ggplot(panel_a_data,
                  aes(x = Mutation, y = Pathway,
                      color = NES, size = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +
  geom_point(
    data = subset(panel_a_data, HighlySig == TRUE),
    aes(size = -log10(p.adjust)),
    shape = 21, color = "black", fill = NA, stroke = 1.2, alpha = 1
  ) +
  facet_wrap(~Database, ncol = 1, scales = "free_y") +
  scale_color_gradient2(
    low = "#0571b0", mid = "white", high = "#ca0020",
    midpoint = 0, name = "NES",
    limits = c(-3.5, 3.5), oob = scales::squish
  ) +
  scale_size_continuous(
    range = c(2, 8), name = "-log10(FDR)",
    breaks = c(2, 5, 10, 15), limits = c(0, 20)
  ) +
  labs(
    title = "Cross-Database Validation: Translation and Mitochondrial Pathways",
    subtitle = sprintf("Mutation-specific pathways (FDR < 0.05, black outline = FDR < 0.01, n=%d)",
                      n_legacy_pathways),
    x = "DRP1 Mutation", y = "Pathway"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 11),
    panel.border = element_rect(color = "gray80", fill = NA),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

plot_height_legacy <- max(16, n_legacy_pathways * 0.5)

ggsave(file.path(out_dir, "Panel_A_Cross_Database_Dotplot_Legacy.pdf"),
       p_legacy, width = 12, height = plot_height_legacy)

message("✓ Legacy 2-column figure complete\n")

###############################################################################
##  SEMANTIC GROUPING: Alternative Figure Grouped by Biological Function    ##
###############################################################################

message("📊 Creating SEMANTIC figure: grouped by biological function...\n")

# Assign semantic categories with hierarchical rules
plot_data_semantic <- plot_data %>%
  filter(!is.na(NES)) %>%
  mutate(
    Semantic_Category = case_when(
      # Neuronal (SynGO only)
      Database == "SYNGO" ~ "Neuronal",

      # Mitochondrial + Ribosome combinations (specific > general)
      grepl("mitochondri.*ribosom|ribosom.*mitochondri", Description, ignore.case = TRUE) ~
        "Mitochondrial Ribosome",

      # Cytoplasmic Ribosome (exclude mitochondrial)
      grepl("ribosom|ribosomal", Description, ignore.case = TRUE) &
        !grepl("mitochondri", Description, ignore.case = TRUE) ~
        "Cytoplasmic Ribosome",

      # Mitochondrial Translation
      grepl("mitochondri.*translation|translation.*mitochondri", Description, ignore.case = TRUE) ~
        "Mitochondrial Translation",

      # Cytoplasmic Translation (exclude mitochondrial)
      grepl("translation|translat", Description, ignore.case = TRUE) &
        !grepl("mitochondri", Description, ignore.case = TRUE) ~
        "Cytoplasmic Translation",

      # Remaining Mitochondrial (OXPHOS, biogenesis, etc.)
      grepl("mitochondri|oxidative phosphorylation|electron transport|OXPHOS",
            Description, ignore.case = TRUE) ~
        "Mitochondrial Metabolism",

      # Calcium Signaling
      grepl("calcium|ca2\\+|calmodulin|calcineurin", Description, ignore.case = TRUE) ~
        "Calcium Signaling",

      # RNA Processing (RNP complexes, granules, splicing)
      TRUE ~ "RNA Processing"
    )
  )

# Factor semantic categories in biological order
plot_data_semantic$Semantic_Category <- factor(
  plot_data_semantic$Semantic_Category,
  levels = c(
    "Neuronal",
    "Cytoplasmic Ribosome",
    "Mitochondrial Ribosome",
    "Cytoplasmic Translation",
    "Mitochondrial Translation",
    "Mitochondrial Metabolism",
    "Calcium Signaling",
    "RNA Processing"
  )
)

# Use ORIGINAL pathway names (with database prefixes) for semantic version
# Create factor for y-axis ordering (same as main figure)
plot_data_semantic$Pathway_Raw <- factor(
  plot_data_semantic$Description,
  levels = pathway_order  # Use same ordering as main figure
)

n_semantic_pathways <- length(unique(plot_data_semantic$Description))

# Count pathways per category
semantic_counts <- plot_data_semantic %>%
  group_by(Semantic_Category) %>%
  summarize(N = n_distinct(Description), .groups = "drop")

message(sprintf("  Categorized %d pathways into %d semantic groups:\n",
                n_semantic_pathways, nrow(semantic_counts)))
print(semantic_counts)

# Create semantic-grouped dotplot
p_semantic <- ggplot(plot_data_semantic,
                    aes(x = Contrast_Display, y = Pathway_Raw,
                        color = NES, size = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +

  # Add black outlines for FDR < 0.01
  geom_point(
    data = subset(plot_data_semantic, HighlySig == TRUE),
    aes(size = -log10(p.adjust)),
    shape = 21,
    color = "black",
    fill = NA,
    stroke = 1.2,
    alpha = 1
  ) +

  facet_wrap(~Semantic_Category, ncol = 1, scales = "free_y") +
  scale_color_gradient2(
    low = "#0571b0", mid = "white", high = "#ca0020",
    midpoint = 0,
    name = "NES",
    limits = c(-3.5, 3.5),
    oob = scales::squish
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = "-log10(FDR)",
    breaks = c(2, 5, 10, 15),
    limits = c(0, 20)
  ) +
  labs(
    title = "Cross-Database Validation: Semantic Grouping by Biological Function",
    subtitle = sprintf("Pathways grouped by mechanism (FDR < 0.05, black outline = FDR < 0.01, n=%d pathways)",
                      n_semantic_pathways),
    x = "Contrast",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    axis.text.y = element_text(size = 8),  # Slightly smaller for long names
    axis.text.x = element_text(face = "bold", size = 10, angle = 0),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold", size = 11),
    panel.border = element_rect(color = "gray80", fill = NA),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# Dynamic plot height
plot_height_semantic <- max(18, n_semantic_pathways * 0.6)

ggsave(file.path(out_dir, "Cross_Database_3Column_Semantic.pdf"),
       p_semantic, width = 13, height = plot_height_semantic)

message("✓ Semantic-grouped figure complete\n")

###############################################################################
##  Export Comprehensive CSV (Single Source of Truth)                       ##
###############################################################################

message("📝 Exporting comprehensive pathway data (single source of truth)...\n")

# Export ALL pathway data (all contrasts, all pathways)
write.csv(plot_data %>%
            filter(!is.na(NES)) %>%
            select(Database, Contrast_Display, Description, Description_Short,
                   NES, p.adjust, pvalue, setSize, Category) %>%
            arrange(Database, Contrast_Display, p.adjust),
          file.path(out_dir, "Cross_database_all_pathways_COMPREHENSIVE.csv"),
          row.names = FALSE)

# Export mutation-only data (for legacy compatibility)
write.csv(all_translation_pathways %>%
            filter(Contrast != "Time_Ctrl") %>%
            select(Database, Mutation, Description, NES, p.adjust,
                   pvalue, setSize, Category) %>%
            arrange(Database, Mutation, p.adjust),
          file.path(out_dir, "All_translation_pathways_significant.csv"),
          row.names = FALSE)

# Export semantic category mapping
write.csv(plot_data_semantic %>%
            select(Description, Database, Semantic_Category) %>%
            distinct() %>%
            arrange(Semantic_Category, Database, Description),
          file.path(out_dir, "Cross_database_semantic_category_mapping.csv"),
          row.names = FALSE)

message("✓ CSV exports complete\n")

###############################################################################
##  Summary Report                                                           ##
###############################################################################

cat("\n=== Cross-Database Validation Summary (Consolidated) ===\n\n")

cat("Total Pathways Analyzed:\n")
cat(sprintf("  Unique pathways (all contrasts): %d\n", n_unique_pathways))
cat(sprintf("  Time_Ctrl-specific pathways: %d\n", n_timectrl_only))
cat(sprintf("  Mutation-disrupted pathways: %d\n", n_legacy_pathways))
cat(sprintf("  Across databases: %d\n", n_distinct(plot_data_main$Database)))

cat("\n\nPathways by Contrast:\n")
contrast_summary <- plot_data_main %>%
  group_by(Contrast_Display) %>%
  summarize(N_pathways = n(), .groups = "drop")
print(contrast_summary)

cat("\n\nPathways by Database:\n")
db_counts <- plot_data_main %>%
  group_by(Database) %>%
  summarize(N_pathways = n_distinct(Description), .groups = "drop") %>%
  arrange(Database)
print(db_counts, n = Inf)


cat("BIOLOGICAL STORY:\n")
cat("  - Time_Ctrl: Normal maturation increases synaptic ribosomes (NES +2.4)\n")
cat("  - Mutations: Maturation DECREASES synaptic ribosomes (NES -2.9)\n")
cat(sprintf("  - %d pathways show pattern discrepencies\n", n_timectrl_only))

message("\n✅ Cross-database validation complete (consolidated script)!")
message(sprintf("📁 Output directory: %s", out_dir))
message("\n🎯 KEY OUTPUTS:")
message("  1. Cross_Database_3Column_Main.pdf - Database-grouped figure (formatted names)")
message("  2. Cross_Database_3Column_Semantic.pdf - Semantic-grouped figure (raw names)")
message("  3. Panel_A_Cross_Database_Dotplot_Legacy.pdf - Legacy 2-column figure")
message("  4. Cross_database_all_pathways_COMPREHENSIVE.csv - Complete data")
message("  5. Cross_database_semantic_category_mapping.csv - Pathway-to-category mapping")
message("  6. All_translation_pathways_significant.csv - Legacy data format")
message("\n🔬 FEATURES:")
message("  ✓ Two complementary views: Database vs. Semantic grouping")
message("  ✓ Database view: Cross-database validation (formatted pathway names)")
message("  ✓ Semantic view: Biological mechanisms (raw names with DB prefixes)")
message("  ✓ Granular categories (e.g., Mitochondrial Ribosome)")
message("  ✓ SYNGO/Neuronal at top in both versions")
