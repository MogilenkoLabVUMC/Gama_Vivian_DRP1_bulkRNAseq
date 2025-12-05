###############################################################################
##  Export GSEA results to CSV for Python visualization                      ##
###############################################################################
##  Purpose: Convert all GSEA results from RDS checkpoints to CSV format
##           for comprehensive cross-database trajectory analysis in Python
###############################################################################

library(here)
library(dplyr)
library(tidyr)

# -------------------------------------------------------------------- #
# 1. Configuration                                                     #
# -------------------------------------------------------------------- #
config <- list(
  out_root = "03_Results/02_Analysis",
  checkpoint_dir = "03_Results/02_Analysis/checkpoints"
)

# -------------------------------------------------------------------- #
# 2. Load GSEA results from checkpoints                                #
# -------------------------------------------------------------------- #
message("📂 Loading GSEA results from checkpoints...")

# Load all GSEA results
all_gsea <- readRDS(here::here(config$checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea <- readRDS(here::here(config$checkpoint_dir, "syngo_gsea_results.rds"))
mitocarta_gsea <- readRDS(here::here(config$checkpoint_dir, "mitocarta_gsea_results.rds"))

# -------------------------------------------------------------------- #
# 3. Define contrast mapping to new framework                          #
# -------------------------------------------------------------------- #
contrast_mapping <- data.frame(
  old_name = c(
    "G32A_vs_Ctrl_D35",
    "R403C_vs_Ctrl_D35",
    "G32A_vs_Ctrl_D65",
    "R403C_vs_Ctrl_D65",
    "Maturation_G32A_specific",
    "Maturation_R403C_specific",
    "Time_Ctrl",
    "Time_G32A",
    "Time_R403C"
  ),
  new_name = c(
    "Early_G32A",
    "Early_R403C",
    "Late_G32A",
    "Late_R403C",
    "TrajDev_G32A",
    "TrajDev_R403C",
    "Time_Ctrl",
    "Time_G32A",
    "Time_R403C"
  ),
  category = c(
    "Early",
    "Early",
    "Late",
    "Late",
    "TrajDev",
    "TrajDev",
    "Reference",
    "Reference",
    "Reference"
  ),
  mutation = c(
    "G32A",
    "R403C",
    "G32A",
    "R403C",
    "G32A",
    "R403C",
    "Control",
    "G32A",
    "R403C"
  ),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------- #
# 4. Function to extract GSEA results to dataframe                     #
# -------------------------------------------------------------------- #
extract_gsea_to_df <- function(gsea_obj, contrast_name, database_name) {
  # Handle NULL results
  if (is.null(gsea_obj)) {
    return(NULL)
  }

  # Extract result dataframe
  if (inherits(gsea_obj, "list") && "result" %in% names(gsea_obj)) {
    result_df <- gsea_obj$result
  } else if (inherits(gsea_obj, "enrichResult") || inherits(gsea_obj, "gseaResult")) {
    result_df <- as.data.frame(gsea_obj@result)
  } else {
    warning("Unknown GSEA object structure for ", contrast_name, " - ", database_name)
    return(NULL)
  }

  # Check if result is empty
  if (is.null(result_df) || nrow(result_df) == 0) {
    return(NULL)
  }

  # Add metadata columns
  result_df$contrast <- contrast_name
  result_df$database <- database_name

  # Select key columns (adjust based on actual structure)
  key_cols <- c("ID", "Description", "setSize", "enrichmentScore", "NES",
                "pvalue", "p.adjust", "qvalue", "contrast", "database")

  # Only keep columns that exist
  existing_cols <- intersect(key_cols, colnames(result_df))
  result_df <- result_df[, existing_cols, drop = FALSE]

  return(result_df)
}

# -------------------------------------------------------------------- #
# 5. Extract all MSigDB GSEA results                                   #
# -------------------------------------------------------------------- #
message("\n🔬 Extracting MSigDB GSEA results...")

msigdb_results <- list()

for (contrast_name in names(all_gsea)) {
  contrast_data <- all_gsea[[contrast_name]]

  if (is.null(contrast_data) || !is.list(contrast_data)) {
    next
  }

  for (db_name in names(contrast_data)) {
    db_data <- contrast_data[[db_name]]

    df <- extract_gsea_to_df(db_data, contrast_name, db_name)

    if (!is.null(df)) {
      msigdb_results[[paste(contrast_name, db_name, sep = "_")]] <- df
    }
  }
}

# Combine all MSigDB results
if (length(msigdb_results) > 0) {
  msigdb_combined <- bind_rows(msigdb_results)
  message("✓ Extracted ", nrow(msigdb_combined), " MSigDB pathway results")
} else {
  msigdb_combined <- data.frame()
  warning("No MSigDB results extracted")
}

# -------------------------------------------------------------------- #
# 6. Extract SynGO GSEA results                                        #
# -------------------------------------------------------------------- #
message("\n🔬 Extracting SynGO GSEA results...")

syngo_results <- list()

for (contrast_name in names(syngo_gsea)) {
  df <- extract_gsea_to_df(syngo_gsea[[contrast_name]], contrast_name, "SynGO")

  if (!is.null(df)) {
    syngo_results[[contrast_name]] <- df
  }
}

if (length(syngo_results) > 0) {
  syngo_combined <- bind_rows(syngo_results)
  message("✓ Extracted ", nrow(syngo_combined), " SynGO pathway results")
} else {
  syngo_combined <- data.frame()
  warning("No SynGO results extracted")
}

# -------------------------------------------------------------------- #
# 7. Extract MitoCarta GSEA results                                    #
# -------------------------------------------------------------------- #
message("\n🔬 Extracting MitoCarta GSEA results...")

mitocarta_results <- list()

for (contrast_name in names(mitocarta_gsea)) {
  df <- extract_gsea_to_df(mitocarta_gsea[[contrast_name]], contrast_name, "MitoCarta")

  if (!is.null(df)) {
    mitocarta_results[[contrast_name]] <- df
  }
}

if (length(mitocarta_results) > 0) {
  mitocarta_combined <- bind_rows(mitocarta_results)
  message("✓ Extracted ", nrow(mitocarta_combined), " MitoCarta pathway results")
} else {
  mitocarta_combined <- data.frame()
  warning("No MitoCarta results extracted")
}

# -------------------------------------------------------------------- #
# 8. Combine all databases                                             #
# -------------------------------------------------------------------- #
message("\n📊 Combining all databases...")

all_combined <- bind_rows(
  msigdb_combined,
  syngo_combined,
  mitocarta_combined
)

# Add new framework mapping
all_combined <- all_combined %>%
  left_join(contrast_mapping, by = c("contrast" = "old_name"))

message("✓ Total pathways across all databases: ", nrow(all_combined))
message("  Contrasts: ", paste(unique(all_combined$contrast), collapse = ", "))
message("  Databases: ", paste(unique(all_combined$database), collapse = ", "))

# -------------------------------------------------------------------- #
# 9. Create wide format for visualization                              #
# -------------------------------------------------------------------- #
message("\n📊 Creating wide format pivot table...")

# Define trajectory contrasts (needed for significance flags)
trajectory_contrasts <- c(
  "G32A_vs_Ctrl_D35", "R403C_vs_Ctrl_D35",  # Early
  "Maturation_G32A_specific", "Maturation_R403C_specific",  # TrajDev
  "G32A_vs_Ctrl_D65", "R403C_vs_Ctrl_D65"   # Late
)

# Create a pathway identifier that's unique across databases
all_combined <- all_combined %>%
  mutate(pathway_id = paste(database, ID, sep = "::"))

# Add significance flags for each pathway
all_combined <- all_combined %>%
  group_by(pathway_id) %>%
  mutate(
    ever_significant = any(p.adjust < 0.05, na.rm = TRUE),
    ever_significant_trajectory = any(
      p.adjust < 0.05 & contrast %in% trajectory_contrasts,
      na.rm = TRUE
    )
  ) %>%
  ungroup()

message("  Added ever_significant flags to all pathways")

# Pivot to wide format: rows = pathways, columns = contrasts
wide_format <- all_combined %>%
  select(pathway_id, database, Description, contrast, NES, p.adjust) %>%
  pivot_wider(
    id_cols = c(pathway_id, database, Description),
    names_from = contrast,
    values_from = c(NES, p.adjust),
    names_sep = "_"
  )

message("✓ Wide format created: ", nrow(wide_format), " unique pathways")

# -------------------------------------------------------------------- #
# 10. Save to CSV files                                                #
# -------------------------------------------------------------------- #
output_dir <- here::here(config$out_root, "Python_exports")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("\n💾 Saving CSV files...")

# Save long format (complete data)
write.csv(
  all_combined,
  file.path(output_dir, "gsea_results_long.csv"),
  row.names = FALSE
)
message("  ✓ Saved: gsea_results_long.csv")

# Save wide format (for heatmap and p.adjust column merging)
write.csv(
  wide_format,
  file.path(output_dir, "gsea_results_wide.csv"),
  row.names = FALSE
)
message("  ✓ Saved: gsea_results_wide.csv")

# -------------------------------------------------------------------- #
# 11. Summary statistics                                               #
# -------------------------------------------------------------------- #
message("\n", paste(rep("=", 78), collapse = ""))
message("📊 GSEA Export Summary")
message(paste(rep("=", 78), collapse = ""))

summary_stats <- all_combined %>%
  group_by(database, contrast) %>%
  summarise(
    total = n(),
    significant = sum(p.adjust < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

message("\n✨ Export complete!")
message("📁 Output directory: ", output_dir)
message("   - gsea_results_long.csv (complete data with ever_significant flags)")
message("   - gsea_results_wide.csv (pivot table with p.adjust columns)")
message("\n📝 Next step: Run 1.5.create_master_pathway_table.py to create master_gsea_table.csv")
