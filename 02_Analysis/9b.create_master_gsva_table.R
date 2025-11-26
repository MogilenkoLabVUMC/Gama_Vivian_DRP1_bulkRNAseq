###############################################################################
##  Master GSVA Table Generator                                              ##
###############################################################################
##  PURPOSE: Create a comprehensive master table of GSVA trajectory analysis ##
##           similar to the master_pathway_table.csv for GSEA results        ##
##                                                                            ##
##  INPUTS:  - GSVA checkpoint (gsva_module_scores.rds)                      ##
##           - QC variables checkpoint (for sample metadata)                 ##
##                                                                            ##
##  OUTPUTS: - master_gsva_table.csv (comprehensive table for exploration)   ##
##           - gsva_statistics_summary.txt (summary statistics)              ##
###############################################################################

library(here)
library(dplyr)
library(tidyr)
library(tibble)

message("=" |> rep(78) |> paste(collapse = ""))
message("Master GSVA Table Generator")
message("=" |> rep(78) |> paste(collapse = ""))

# ============================================================================ #
# 0. Configuration                                                             #
# ============================================================================ #

config <- list(
  checkpoint_dir = here("03_Results/02_Analysis/checkpoints"),
  out_root = here("03_Results/02_Analysis"),

  # Significance thresholds
  p_cutoff = 0.05,

  # Module metadata
  module_info = list(
    list(key = "Ribosome_Biogenesis", display = "Ribosome Biogenesis",
         source = "GO:BP", panel = "A"),
    list(key = "Cytoplasmic_Translation", display = "Cytoplasmic Translation",
         source = "GO:BP", panel = "B"),
    list(key = "Synaptic_Ribosomes", display = "Synaptic Ribosomes",
         source = "SynGO", panel = "C"),
    list(key = "Mitochondrial_Ribosome", display = "Mitochondrial Ribosome",
         source = "MitoCarta", panel = "D"),
    list(key = "Mito_Ribosome_Assembly", display = "Mito Ribosome Assembly",
         source = "MitoCarta", panel = "E"),
    list(key = "mtDNA_Maintenance", display = "mtDNA Maintenance",
         source = "MitoCarta", panel = "F"),
    list(key = "OXPHOS", display = "OXPHOS",
         source = "MitoCarta", panel = "G")
  )
)

# ============================================================================ #
# 1. Load Checkpoint Data                                                      #
# ============================================================================ #

message("\n📂 Loading checkpoints...")

# Load GSVA scores
gsva_checkpoint <- readRDS(file.path(config$checkpoint_dir, "gsva_module_scores.rds"))
gsva_scores <- gsva_checkpoint$scores
gene_modules_filtered <- gsva_checkpoint$gene_modules_filtered

# Load QC variables for sample metadata
qc_vars <- readRDS(file.path(config$checkpoint_dir, "qc_variables.rds"))
annot <- qc_vars$annot

message("✓ Checkpoints loaded")

# ============================================================================ #
# 2. Transform GSVA Scores to Long Format with Metadata                       #
# ============================================================================ #

message("\n📊 Transforming GSVA scores to long format...")

# Convert to long format
gsva_long <- gsva_scores %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Module") %>%
  tidyr::pivot_longer(
    cols = -Module,
    names_to = "Sample",
    values_to = "GSVA_Score"
  )

# Add sample metadata
annot_df <- annot %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Genotype = gsub("Control", "Ctrl", genotype),
    Genotype = factor(Genotype, levels = c("Ctrl", "G32A", "R403C")),
    Day = ifelse(days == "D35", 35, 65),
    Timepoint = days
  )

gsva_long <- gsva_long %>%
  left_join(annot_df, by = "Sample")

message(sprintf("  ✓ Long format: %d observations (7 modules × %d samples)",
                nrow(gsva_long), length(unique(gsva_long$Sample))))

# ============================================================================ #
# 3. Calculate Group-Level Statistics                                         #
# ============================================================================ #

message("\n📊 Calculating group-level statistics...")

# Calculate mean, SD, SE for each module-genotype-timepoint combination
gsva_group_stats <- gsva_long %>%
  group_by(Module, Genotype, Day, Timepoint) %>%
  summarize(
    Mean_GSVA = mean(GSVA_Score, na.rm = TRUE),
    SD_GSVA = sd(GSVA_Score, na.rm = TRUE),
    SE_GSVA = SD_GSVA / sqrt(n()),
    N = n(),
    .groups = "drop"
  )

# Calculate Ctrl D35 baseline for trajectory analysis
ctrl_d35_baseline <- gsva_group_stats %>%
  filter(Genotype == "Ctrl", Day == 35) %>%
  select(Module, Baseline_CtrlD35 = Mean_GSVA)

# Calculate control mean at each timepoint for divergence analysis
ctrl_by_day <- gsva_group_stats %>%
  filter(Genotype == "Ctrl") %>%
  select(Module, Day, Ctrl_Mean_Same_Day = Mean_GSVA)

# Add trajectory and divergence metrics
gsva_group_stats <- gsva_group_stats %>%
  left_join(ctrl_d35_baseline, by = "Module") %>%
  left_join(ctrl_by_day, by = c("Module", "Day")) %>%
  mutate(
    # Trajectory: relative to Ctrl D35
    Expression_vs_CtrlD35 = Mean_GSVA - Baseline_CtrlD35,

    # Divergence: relative to control at same timepoint
    Divergence_vs_Ctrl = ifelse(Genotype == "Ctrl", 0, Mean_GSVA - Ctrl_Mean_Same_Day)
  )

message(sprintf("  ✓ Group statistics: %d groups", nrow(gsva_group_stats)))

# ============================================================================ #
# 4. Statistical Testing (T-tests)                                            #
# ============================================================================ #

message("\n📊 Performing statistical tests...")

# Function to perform t-test for a given module-genotype-timepoint
perform_ttest <- function(module, genotype, day, gsva_data) {

  # Skip control (nothing to compare against)
  if (genotype == "Ctrl") {
    return(data.frame(
      Module = module,
      Genotype = genotype,
      Day = day,
      t_statistic = NA,
      p_value = NA,
      p_adjusted = NA,
      significant = FALSE
    ))
  }

  # Get mutant samples
  mutant_samples <- gsva_data %>%
    filter(Module == !!module, Genotype == !!genotype, Day == !!day) %>%
    pull(GSVA_Score)

  # Get control samples at same timepoint
  ctrl_samples <- gsva_data %>%
    filter(Module == !!module, Genotype == "Ctrl", Day == !!day) %>%
    pull(GSVA_Score)

  # Perform t-test if we have enough samples
  if (length(mutant_samples) >= 2 && length(ctrl_samples) >= 2) {
    test_result <- t.test(mutant_samples, ctrl_samples)

    return(data.frame(
      Module = module,
      Genotype = genotype,
      Day = day,
      t_statistic = test_result$statistic,
      p_value = test_result$p.value,
      p_adjusted = NA,  # Will adjust later
      significant = FALSE  # Will update later
    ))
  } else {
    return(data.frame(
      Module = module,
      Genotype = genotype,
      Day = day,
      t_statistic = NA,
      p_value = NA,
      p_adjusted = NA,
      significant = FALSE
    ))
  }
}

# Perform tests for all module-genotype-timepoint combinations
test_combinations <- gsva_group_stats %>%
  filter(Genotype != "Ctrl") %>%
  select(Module, Genotype, Day) %>%
  distinct()

ttest_results <- purrr::pmap_dfr(
  test_combinations,
  function(Module, Genotype, Day) {
    perform_ttest(Module, Genotype, Day, gsva_long)
  }
)

# Adjust p-values for multiple testing (Benjamini-Hochberg)
ttest_results <- ttest_results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = !is.na(p_adjusted) & p_adjusted < config$p_cutoff
  )

message(sprintf("  ✓ T-tests performed: %d comparisons", nrow(ttest_results)))
message(sprintf("  ✓ Significant comparisons (p.adj < 0.05): %d",
                sum(ttest_results$significant, na.rm = TRUE)))

# ============================================================================ #
# 5. Add Module Metadata                                                      #
# ============================================================================ #

message("\n📊 Adding module metadata...")

# Create module metadata dataframe
module_metadata <- data.frame(
  Module = sapply(config$module_info, function(x) x$key),
  Display_Name = sapply(config$module_info, function(x) x$display),
  Source_Database = sapply(config$module_info, function(x) x$source),
  Panel_ID = sapply(config$module_info, function(x) x$panel),
  N_genes = sapply(names(gene_modules_filtered), function(m) length(gene_modules_filtered[[m]]))
)

message(sprintf("  ✓ Module metadata created for %d modules", nrow(module_metadata)))

# ============================================================================ #
# 6. Create Master Table                                                      #
# ============================================================================ #

message("\n📊 Creating master GSVA table...")

# Merge all components
master_table <- gsva_group_stats %>%
  left_join(module_metadata, by = "Module") %>%
  left_join(ttest_results, by = c("Module", "Genotype", "Day")) %>%
  arrange(Module, Day, Genotype)

# Add trajectory categories (similar to GSEA contrasts)
master_table <- master_table %>%
  mutate(
    Trajectory_Category = case_when(
      Day == 35 & Genotype != "Ctrl" ~ "Early",
      Day == 65 & Genotype != "Ctrl" ~ "Late",
      Genotype == "Ctrl" ~ "Reference",
      TRUE ~ "Other"
    ),

    Contrast_Label = case_when(
      Day == 35 & Genotype == "G32A" ~ "Early_G32A",
      Day == 35 & Genotype == "R403C" ~ "Early_R403C",
      Day == 65 & Genotype == "G32A" ~ "Late_G32A",
      Day == 65 & Genotype == "R403C" ~ "Late_R403C",
      Day == 35 & Genotype == "Ctrl" ~ "Ctrl_D35",
      Day == 65 & Genotype == "Ctrl" ~ "Ctrl_D65",
      TRUE ~ NA_character_
    )
  )

# Reorder columns for clarity
master_table <- master_table %>%
  select(
    # Module identification
    Module,
    Display_Name,
    Panel_ID,
    Source_Database,
    N_genes,

    # Sample grouping
    Genotype,
    Day,
    Timepoint,
    Trajectory_Category,
    Contrast_Label,

    # GSVA scores
    Mean_GSVA,
    SD_GSVA,
    SE_GSVA,
    N,

    # Trajectory metrics
    Baseline_CtrlD35,
    Expression_vs_CtrlD35,

    # Divergence metrics
    Ctrl_Mean_Same_Day,
    Divergence_vs_Ctrl,

    # Statistical tests
    t_statistic,
    p_value,
    p_adjusted,
    significant
  )

message(sprintf("  ✓ Master table created: %d rows × %d columns",
                nrow(master_table), ncol(master_table)))

# ============================================================================ #
# 7. Create Wide-Format Summary for Pattern Analysis                          #
# ============================================================================ #

message("\n📊 Creating pattern summary table...")

# Pivot to wide format for easier pattern analysis
pattern_summary <- master_table %>%
  filter(Genotype != "Ctrl") %>%
  select(Module, Display_Name, Source_Database, Genotype, Day,
         Expression_vs_CtrlD35, Divergence_vs_Ctrl, significant) %>%
  pivot_wider(
    id_cols = c(Module, Display_Name, Source_Database),
    names_from = c(Genotype, Day),
    values_from = c(Expression_vs_CtrlD35, Divergence_vs_Ctrl, significant),
    names_sep = "_"
  )

# Add pattern classifications based on trajectory
pattern_summary <- pattern_summary %>%
  mutate(
    # G32A pattern
    Pattern_G32A = case_when(
      is.na(Expression_vs_CtrlD35_G32A_35) | is.na(Expression_vs_CtrlD35_G32A_65) ~ "Insufficient_data",

      # Compensation: starts down/neutral, ends up
      Expression_vs_CtrlD35_G32A_35 <= 0 & Expression_vs_CtrlD35_G32A_65 > 0.2 ~ "Compensation",

      # Progressive: starts up, gets more up
      Expression_vs_CtrlD35_G32A_35 > 0 & Expression_vs_CtrlD35_G32A_65 > Expression_vs_CtrlD35_G32A_35 ~ "Progressive",

      # Persistent: consistently dysregulated (same direction)
      sign(Expression_vs_CtrlD35_G32A_35) == sign(Expression_vs_CtrlD35_G32A_65) &
        abs(Expression_vs_CtrlD35_G32A_65 - Expression_vs_CtrlD35_G32A_35) < 0.1 ~ "Persistent",

      # Transient: starts dysregulated, returns to baseline
      abs(Expression_vs_CtrlD35_G32A_35) > 0.2 & abs(Expression_vs_CtrlD35_G32A_65) < 0.1 ~ "Transient",

      # Natural worsening: starts up, ends down (or vice versa with worsening)
      Expression_vs_CtrlD35_G32A_35 > 0 & Expression_vs_CtrlD35_G32A_65 < 0 ~ "Natural_worsening",

      TRUE ~ "Complex"
    ),

    # R403C pattern (same logic)
    Pattern_R403C = case_when(
      is.na(Expression_vs_CtrlD35_R403C_35) | is.na(Expression_vs_CtrlD35_R403C_65) ~ "Insufficient_data",

      Expression_vs_CtrlD35_R403C_35 <= 0 & Expression_vs_CtrlD35_R403C_65 > 0.2 ~ "Compensation",

      Expression_vs_CtrlD35_R403C_35 > 0 & Expression_vs_CtrlD35_R403C_65 > Expression_vs_CtrlD35_R403C_35 ~ "Progressive",

      sign(Expression_vs_CtrlD35_R403C_35) == sign(Expression_vs_CtrlD35_R403C_65) &
        abs(Expression_vs_CtrlD35_R403C_65 - Expression_vs_CtrlD35_R403C_35) < 0.1 ~ "Persistent",

      abs(Expression_vs_CtrlD35_R403C_35) > 0.2 & abs(Expression_vs_CtrlD35_R403C_65) < 0.1 ~ "Transient",

      Expression_vs_CtrlD35_R403C_35 > 0 & Expression_vs_CtrlD35_R403C_65 < 0 ~ "Natural_worsening",

      TRUE ~ "Complex"
    ),

    # Change consistency
    Change_Consistency = case_when(
      Pattern_G32A == Pattern_R403C ~ paste0("Consistent_", Pattern_G32A),
      Pattern_G32A != Pattern_R403C ~ paste0("Inconsistent_", Pattern_G32A, "_vs_", Pattern_R403C),
      TRUE ~ "Unknown"
    )
  )

message(sprintf("  ✓ Pattern summary created: %d modules", nrow(pattern_summary)))

# ============================================================================ #
# 8. Export Master Tables                                                     #
# ============================================================================ #

message("\n💾 Exporting master tables...")

# Export main master table (long format - one row per module-genotype-timepoint)
master_file <- file.path(config$out_root, "master_gsva_table.csv")
write.csv(master_table, master_file, row.names = FALSE)
message(sprintf("  ✓ master_gsva_table.csv (%d rows × %d columns)",
                nrow(master_table), ncol(master_table)))

# Export pattern summary (wide format - one row per module)
pattern_file <- file.path(config$out_root, "gsva_pattern_summary.csv")
write.csv(pattern_summary, pattern_file, row.names = FALSE)
message(sprintf("  ✓ gsva_pattern_summary.csv (%d rows × %d columns)",
                nrow(pattern_summary), ncol(pattern_summary)))

# ============================================================================ #
# 9. Generate Summary Statistics                                              #
# ============================================================================ #

message("\n📊 Generating summary statistics...")

summary_text <- c(
  "=" |> rep(78) |> paste(collapse = ""),
  "MASTER GSVA TABLE SUMMARY",
  "=" |> rep(78) |> paste(collapse = ""),
  "",
  sprintf("Generated: %s", Sys.time()),
  "",
  "--- MASTER TABLE (master_gsva_table.csv) ---",
  sprintf("Total rows: %d", nrow(master_table)),
  sprintf("Columns: %d", ncol(master_table)),
  sprintf("Modules: %d", length(unique(master_table$Module))),
  sprintf("Genotypes: %s", paste(unique(master_table$Genotype), collapse = ", ")),
  sprintf("Timepoints: %s", paste(unique(master_table$Timepoint), collapse = ", ")),
  "",
  "--- STATISTICAL TESTS ---",
  sprintf("Total comparisons (mutants vs controls): %d", nrow(ttest_results)),
  sprintf("Significant comparisons (p.adj < 0.05): %d", sum(ttest_results$significant, na.rm = TRUE)),
  "",
  "Significant comparisons by module:",
  tapply(ttest_results$significant, ttest_results$Module, sum, na.rm = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("Module") %>%
    setNames(c("Module", "N_significant")) %>%
    arrange(desc(N_significant)) %>%
    apply(1, function(x) sprintf("  %s: %s", x[1], x[2])) %>%
    paste(collapse = "\n"),
  "",
  "--- PATTERN SUMMARY (gsva_pattern_summary.csv) ---",
  sprintf("Modules analyzed: %d", nrow(pattern_summary)),
  "",
  "G32A patterns:",
  table(pattern_summary$Pattern_G32A) %>%
    as.data.frame() %>%
    setNames(c("Pattern", "Count")) %>%
    apply(1, function(x) sprintf("  %s: %s", x[1], x[2])) %>%
    paste(collapse = "\n"),
  "",
  "R403C patterns:",
  table(pattern_summary$Pattern_R403C) %>%
    as.data.frame() %>%
    setNames(c("Pattern", "Count")) %>%
    apply(1, function(x) sprintf("  %s: %s", x[1], x[2])) %>%
    paste(collapse = "\n"),
  "",
  "Change consistency:",
  table(pattern_summary$Change_Consistency) %>%
    as.data.frame() %>%
    setNames(c("Consistency", "Count")) %>%
    apply(1, function(x) sprintf("  %s: %s", x[1], x[2])) %>%
    paste(collapse = "\n"),
  "",
  "--- MODULE GENE COUNTS ---",
  module_metadata %>%
    arrange(desc(N_genes)) %>%
    apply(1, function(x) sprintf("  %s (%s): %s genes", x[2], x[3], x[5])) %>%
    paste(collapse = "\n"),
  "",
  "=" |> rep(78) |> paste(collapse = ""),
  "KEY COLUMNS IN MASTER TABLE:",
  "=" |> rep(78) |> paste(collapse = ""),
  "",
  "Identification:",
  "  - Module: Internal module key",
  "  - Display_Name: Human-readable module name",
  "  - Source_Database: Gene set source (GO:BP, SynGO, MitoCarta)",
  "",
  "GSVA Scores:",
  "  - Mean_GSVA: Mean GSVA enrichment score for this group",
  "  - SD_GSVA: Standard deviation",
  "  - SE_GSVA: Standard error",
  "",
  "Trajectory Metrics (relative to Ctrl D35 baseline):",
  "  - Baseline_CtrlD35: Control D35 baseline score",
  "  - Expression_vs_CtrlD35: Deviation from baseline",
  "",
  "Divergence Metrics (relative to matched control):",
  "  - Ctrl_Mean_Same_Day: Control mean at same timepoint",
  "  - Divergence_vs_Ctrl: Difference from control",
  "",
  "Statistical Tests:",
  "  - t_statistic: T-test statistic (mutant vs control)",
  "  - p_value: Raw p-value",
  "  - p_adjusted: Benjamini-Hochberg adjusted p-value",
  "  - significant: TRUE if p.adj < 0.05",
  "",
  "=" |> rep(78) |> paste(collapse = ""),
  "USAGE EXAMPLES:",
  "=" |> rep(78) |> paste(collapse = ""),
  "",
  "# Load in R:",
  "df <- read.csv('03_Results/02_Analysis/master_gsva_table.csv')",
  "",
  "# Find significant G32A changes at D65:",
  "df %>% filter(Genotype == 'G32A', Day == 65, significant == TRUE)",
  "",
  "# Compare trajectory patterns:",
  "df %>% filter(Module == 'OXPHOS') %>% ",
  "  select(Genotype, Day, Expression_vs_CtrlD35, Divergence_vs_Ctrl)",
  "",
  "=" |> rep(78) |> paste(collapse = "")
)

summary_file <- file.path(config$out_root, "gsva_statistics_summary.txt")
writeLines(summary_text, summary_file)
message(sprintf("  ✓ gsva_statistics_summary.txt"))

# ============================================================================ #
# 10. Display Summary                                                         #
# ============================================================================ #

message("\n" |> paste0(rep("=", 78) |> paste(collapse = "")))
message("✅ MASTER GSVA TABLE GENERATION COMPLETE")
message(rep("=", 78) |> paste(collapse = ""))

message("\n📁 Output files:")
message(sprintf("  1. %s", basename(master_file)))
message(sprintf("     - %d rows (module × genotype × timepoint combinations)", nrow(master_table)))
message(sprintf("     - %d columns", ncol(master_table)))
message(sprintf("\n  2. %s", basename(pattern_file)))
message(sprintf("     - %d modules with pattern classifications", nrow(pattern_summary)))
message(sprintf("\n  3. %s", basename(summary_file)))
message("     - Summary statistics and usage guide")

message("\n📊 Key findings:")
message(sprintf("  - Significant comparisons: %d/%d (%.1f%%)",
                sum(ttest_results$significant, na.rm = TRUE),
                nrow(ttest_results),
                100 * sum(ttest_results$significant, na.rm = TRUE) / nrow(ttest_results)))

message("\n📈 Pattern distribution:")
message("  G32A:")
for (pattern in names(table(pattern_summary$Pattern_G32A))) {
  count <- sum(pattern_summary$Pattern_G32A == pattern)
  message(sprintf("    - %s: %d modules", pattern, count))
}
message("  R403C:")
for (pattern in names(table(pattern_summary$Pattern_R403C))) {
  count <- sum(pattern_summary$Pattern_R403C == pattern)
  message(sprintf("    - %s: %d modules", pattern, count))
}

message("\n✨ Done!")
message("Use these tables to explore GSVA trajectory patterns and identify")
message("modules with interesting developmental dynamics!\n")
