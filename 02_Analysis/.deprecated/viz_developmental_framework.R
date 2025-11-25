###############################################################################
##  Developmental Framework Visualization                                    ##
##  Early (baseline) → TrajDev (trajectory) → Late (outcome)                 ##
###############################################################################
##  Purpose: Reframe cross-database validation around developmental story
##  - Early: Baseline defect at D35 (G32A_vs_Ctrl_D35, R403C_vs_Ctrl_D35)
##  - TrajDev: Trajectory deviation (Maturation_*_specific)
##  - Late: Final outcome at D65 (G32A_vs_Ctrl_D65, R403C_vs_Ctrl_D65)
##  - Reference: Normal maturation (Time_Ctrl) shown as marker in TrajDev
##  - SuccessScore: -log2(|Late|/|Early|) shown as annotation in Late column
###############################################################################

# -------------------------------------------------------------------- #
# 0. Configuration & Setup                                             #
# -------------------------------------------------------------------- #

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)
library(plotly)

# Match main pipeline config
config <- list(
  out_root      = "03_Results/02_Analysis",
  helper_root   = "01_Scripts/RNAseq-toolkit",
  p_cutoff      = 0.05,
  fc_cutoff     = 2
)

# Output directory
dev_framework_dir <- here::here(config$out_root, "Plots/Developmental_Framework")
dir.create(dev_framework_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dev_framework_dir, "By_Category"), recursive = TRUE, showWarnings = FALSE)

# Checkpoint directory
checkpoint_dir <- here::here(config$out_root, "checkpoints")

message("🚀 Starting Developmental Framework Visualization")
message("📁 Output directory: ", dev_framework_dir)

# -------------------------------------------------------------------- #
# 1. Load GSEA Results from Checkpoints                               #
# -------------------------------------------------------------------- #

message("\n📦 Loading GSEA results from checkpoints...")

# Load all GSEA results
all_gsea_results <- readRDS(file.path(checkpoint_dir, "all_gsea_results.rds"))
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))
mitocarta_gsea_results <- readRDS(file.path(checkpoint_dir, "mitocarta_gsea_results.rds"))

message("✓ Loaded MSigDB results: ", length(all_gsea_results), " contrasts")
message("✓ Loaded SynGO results: ", length(syngo_gsea_results), " contrasts")
message("✓ Loaded MitoCarta results: ", length(mitocarta_gsea_results), " contrasts")

# -------------------------------------------------------------------- #
# 2. Define Contrast Mapping for Developmental Framework              #
# -------------------------------------------------------------------- #

# Developmental framework mapping
dev_mapping <- list(
  # Early baseline defect at D35
  Early = list(
    G32A  = "G32A_vs_Ctrl_D35",
    R403C = "R403C_vs_Ctrl_D35"
  ),
  # Trajectory deviation (interaction effects)
  TrajDev = list(
    G32A  = "Maturation_G32A_specific",
    R403C = "Maturation_R403C_specific",
    Normal = "Time_Ctrl"  # Reference for normal maturation
  ),
  # Late outcome at D65
  Late = list(
    G32A  = "G32A_vs_Ctrl_D65",
    R403C = "R403C_vs_Ctrl_D65"
  )
)

message("\n🗺️  Developmental framework mapping:")
message("  Early:   G32A_vs_Ctrl_D35, R403C_vs_Ctrl_D35")
message("  TrajDev: Maturation_G32A_specific, Maturation_R403C_specific")
message("  Late:    G32A_vs_Ctrl_D65, R403C_vs_Ctrl_D65")
message("  Reference: Time_Ctrl (normal maturation)")

# -------------------------------------------------------------------- #
# 3. Helper Function: Extract NES Values from GSEA Results            #
# -------------------------------------------------------------------- #

#' Extract NES values for a specific pathway from GSEA results
#' @param pathway_name Name of the pathway to extract
#' @param contrast_name Name of the contrast (e.g., "G32A_vs_Ctrl_D35")
#' @param database_name Database name (e.g., "gobp", "SynGO", "MitoCarta")
#' @param gsea_results GSEA results list
#' @return Named vector with NES, padj, or NULL if not found
extract_pathway_nes <- function(pathway_name, contrast_name, database_name, gsea_results) {
  # Handle different database structures
  if (database_name == "SynGO") {
    # SynGO results are stored directly
    result <- gsea_results[[contrast_name]]
  } else if (database_name == "MitoCarta") {
    # MitoCarta results are stored directly
    result <- gsea_results[[contrast_name]]
  } else {
    # MSigDB results are nested by database
    result <- gsea_results[[contrast_name]][[database_name]]
  }

  # Check if result exists and has data
  if (is.null(result) || nrow(result@result) == 0) {
    return(c(NES = NA, padj = NA, pvalue = NA))
  }

  # Find pathway in results (handle partial matching)
  pathway_idx <- which(result@result$Description == pathway_name)
  if (length(pathway_idx) == 0) {
    # Try partial match
    pathway_idx <- grep(pathway_name, result@result$Description, fixed = TRUE)
  }

  if (length(pathway_idx) == 0) {
    return(c(NES = NA, padj = NA, pvalue = NA))
  }

  # Return first match if multiple
  pathway_idx <- pathway_idx[1]

  return(c(
    NES = result@result$NES[pathway_idx],
    padj = result@result$p.adjust[pathway_idx],
    pvalue = result@result$pvalue[pathway_idx]
  ))
}

#' Extract developmental trajectory NES values for a pathway
#' @param pathway_info List with pathway name, database, and display name
#' @return Data frame with developmental trajectory
extract_developmental_nes <- function(pathway_info) {
  pathway_name <- pathway_info$pathway
  database <- pathway_info$database
  display_name <- pathway_info$display_name

  # Select appropriate GSEA results based on database
  if (database == "SynGO") {
    gsea_res <- syngo_gsea_results
  } else if (database == "MitoCarta") {
    gsea_res <- mitocarta_gsea_results
  } else {
    gsea_res <- all_gsea_results
  }

  # Extract NES for each developmental stage and mutation
  results <- list()

  # G32A trajectory
  g32a_early <- extract_pathway_nes(pathway_name, dev_mapping$Early$G32A, database, gsea_res)
  g32a_trajdev <- extract_pathway_nes(pathway_name, dev_mapping$TrajDev$G32A, database, gsea_res)
  g32a_late <- extract_pathway_nes(pathway_name, dev_mapping$Late$G32A, database, gsea_res)

  # R403C trajectory
  r403c_early <- extract_pathway_nes(pathway_name, dev_mapping$Early$R403C, database, gsea_res)
  r403c_trajdev <- extract_pathway_nes(pathway_name, dev_mapping$TrajDev$R403C, database, gsea_res)
  r403c_late <- extract_pathway_nes(pathway_name, dev_mapping$Late$R403C, database, gsea_res)

  # Normal maturation reference (Time_Ctrl)
  time_ctrl <- extract_pathway_nes(pathway_name, dev_mapping$TrajDev$Normal, database, gsea_res)

  # Create tidy data frame
  data.frame(
    pathway = display_name,
    database = database,
    category = pathway_info$category,

    # G32A
    mutation = "G32A",
    Early_NES = g32a_early["NES"],
    Early_padj = g32a_early["padj"],
    TrajDev_NES = g32a_trajdev["NES"],
    TrajDev_padj = g32a_trajdev["padj"],
    Late_NES = g32a_late["NES"],
    Late_padj = g32a_late["padj"],
    TimeCtrl_NES = time_ctrl["NES"],
    TimeCtrl_padj = time_ctrl["padj"],
    stringsAsFactors = FALSE
  ) |>
  bind_rows(
    data.frame(
      pathway = display_name,
      database = database,
      category = pathway_info$category,

      # R403C
      mutation = "R403C",
      Early_NES = r403c_early["NES"],
      Early_padj = r403c_early["padj"],
      TrajDev_NES = r403c_trajdev["NES"],
      TrajDev_padj = r403c_trajdev["padj"],
      Late_NES = r403c_late["NES"],
      Late_padj = r403c_late["padj"],
      TimeCtrl_NES = time_ctrl["NES"],
      TimeCtrl_padj = time_ctrl["padj"],
      stringsAsFactors = FALSE
    )
  )
}

# -------------------------------------------------------------------- #
# 4. Define Pathway Selection (Comprehensive + Mechanistic Focus)     #
# -------------------------------------------------------------------- #

message("\n🎯 Defining pathway selection...")

# Comprehensive pathway list combining mechanistic story + cross-validation
# NOTE: Pathway names must EXACTLY match GSEA result Description field
pathway_list <- list(
  # ===== NEURONAL (SynGO) =====
  list(
    pathway = "postsynaptic ribosome",
    database = "SynGO",
    display_name = "Postsynaptic Ribosome",
    category = "Neuronal"
  ),
  list(
    pathway = "presynaptic ribosome",
    database = "SynGO",
    display_name = "Presynaptic Ribosome",
    category = "Neuronal"
  ),
  list(
    pathway = "synapse",
    database = "SynGO",
    display_name = "Synapse",
    category = "Neuronal"
  ),

  # ===== CYTOPLASMIC TRANSLATION =====
  list(
    pathway = "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
    database = "reactome",
    display_name = "Translation Initiation (Reactome)",
    category = "Cytoplasmic Translation"
  ),
  list(
    pathway = "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
    database = "reactome",
    display_name = "Translation Elongation (Reactome)",
    category = "Cytoplasmic Translation"
  ),
  list(
    pathway = "GOCC_CYTOSOLIC_RIBOSOME",
    database = "gocc",
    display_name = "Cytosolic Ribosome (GOCC)",
    category = "Cytoplasmic Translation"
  ),
  list(
    pathway = "GOBP_TRANSLATIONAL_INITIATION",
    database = "gobp",
    display_name = "Translational Initiation (GOBP)",
    category = "Cytoplasmic Translation"
  ),
  list(
    pathway = "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION",
    database = "kegg",
    display_name = "Translation Initiation (KEGG)",
    category = "Cytoplasmic Translation"
  ),

  # ===== MITOCHONDRIAL TRANSLATION =====
  list(
    pathway = "Mitochondrial_central_dogma",
    database = "MitoCarta",
    display_name = "Mitochondrial Central Dogma",
    category = "Mitochondrial Translation"
  ),
  list(
    pathway = "GOBP_MITOCHONDRIAL_TRANSLATION",
    database = "gobp",
    display_name = "Mitochondrial Translation (GOBP)",
    category = "Mitochondrial Translation"
  ),
  list(
    pathway = "GOCC_ORGANELLAR_RIBOSOME",
    database = "gocc",
    display_name = "Organellar Ribosome (GOCC)",
    category = "Mitochondrial Translation"
  ),
  list(
    pathway = "Mitochondrial_ribosome",
    database = "MitoCarta",
    display_name = "Mitochondrial Ribosome (MitoCarta)",
    category = "Mitochondrial Translation"
  ),

  # ===== MITOCHONDRIAL METABOLISM =====
  list(
    pathway = "Complex_V",
    database = "MitoCarta",
    display_name = "Complex V - ATP Synthase",
    category = "Mitochondrial Metabolism"
  ),
  list(
    pathway = "OXPHOS",
    database = "MitoCarta",
    display_name = "OXPHOS",
    category = "Mitochondrial Metabolism"
  ),
  list(
    pathway = "Complex_I",
    database = "MitoCarta",
    display_name = "Complex I",
    category = "Mitochondrial Metabolism"
  ),

  # ===== CALCIUM SIGNALING =====
  list(
    pathway = "GOBP_REGULATION_OF_CALCIUM_ION_TRANSPORT",
    database = "gobp",
    display_name = "Calcium Ion Transport Regulation",
    category = "Calcium Signaling"
  ),
  list(
    pathway = "GOBP_CALCIUM_ION_HOMEOSTASIS",
    database = "gobp",
    display_name = "Calcium Ion Homeostasis",
    category = "Calcium Signaling"
  )
)

message("✓ Selected ", length(pathway_list), " pathways across 5 categories:")
message("  - Neuronal: ", sum(sapply(pathway_list, function(x) x$category == "Neuronal")))
message("  - Cytoplasmic Translation: ", sum(sapply(pathway_list, function(x) x$category == "Cytoplasmic Translation")))
message("  - Mitochondrial Translation: ", sum(sapply(pathway_list, function(x) x$category == "Mitochondrial Translation")))
message("  - Mitochondrial Metabolism: ", sum(sapply(pathway_list, function(x) x$category == "Mitochondrial Metabolism")))
message("  - Calcium Signaling: ", sum(sapply(pathway_list, function(x) x$category == "Calcium Signaling")))

message("\n⚠️  NOTE: Pathway names must match GSEA results exactly (case-sensitive)")
message("   If pathways show as NA, check GSEA result Description field")

# -------------------------------------------------------------------- #
# 5. Extract Developmental NES Values for All Pathways                #
# -------------------------------------------------------------------- #

message("\n🔬 Extracting developmental NES values...")

# Extract data for all pathways
dev_data <- bind_rows(lapply(pathway_list, extract_developmental_nes))

# Clean up numeric columns
numeric_cols <- c("Early_NES", "Early_padj", "TrajDev_NES", "TrajDev_padj",
                  "Late_NES", "Late_padj", "TimeCtrl_NES", "TimeCtrl_padj")
for (col in numeric_cols) {
  dev_data[[col]] <- as.numeric(dev_data[[col]])
}

message("✓ Extracted data for ", nrow(dev_data) / 2, " pathways × 2 mutations = ", nrow(dev_data), " rows")

# -------------------------------------------------------------------- #
# 6. Calculate SuccessScore                                            #
# -------------------------------------------------------------------- #

#' Calculate SuccessScore: -log2(|Late|/|Early|)
#' >0 = improvement, <0 = worsening, ≈0 = unchanged
calculate_success_score <- function(early_nes, late_nes) {
  # Handle edge cases
  if (is.na(early_nes) || is.na(late_nes)) return(NA)
  if (abs(early_nes) < 0.1) return(NA)  # Early effect too small to interpret

  score <- -log2(abs(late_nes) / abs(early_nes))
  return(score)
}

#' Convert SuccessScore to annotation symbol
#' @param score SuccessScore value
#' @return Character: ↑ (improvement), ↓ (worsening), − (unchanged), or "" (NA)
score_to_symbol <- function(score) {
  if (is.na(score)) return("")
  if (score > 0.5) return("↑")      # Improvement (>1.4× better)
  if (score < -0.5) return("↓")     # Worsening (>1.4× worse)
  return("−")                        # Unchanged
}

message("\n📊 calculating SuccessScores...")

dev_data <- dev_data |>
  mutate(
    SuccessScore = mapply(calculate_success_score, Early_NES, Late_NES),
    SuccessSymbol = sapply(SuccessScore, score_to_symbol)
  )

message("✓ SuccessScore calculated")
message("  ↑ Improvement: ", sum(dev_data$SuccessSymbol == "↑", na.rm = TRUE))
message("  ↓ Worsening:   ", sum(dev_data$SuccessSymbol == "↓", na.rm = TRUE))
message("  − Unchanged:   ", sum(dev_data$SuccessSymbol == "−", na.rm = TRUE))

# -------------------------------------------------------------------- #
# 7. Reshape Data for Plotting                                         #
# -------------------------------------------------------------------- #

message("\n🔄 Reshaping data for visualization...")

# Create long format for ggplot
dev_data_long <- dev_data |>
  select(pathway, database, category, mutation,
         Early_NES, Early_padj,
         TrajDev_NES, TrajDev_padj,
         Late_NES, Late_padj,
         TimeCtrl_NES, SuccessSymbol) |>
  pivot_longer(
    cols = c(Early_NES, TrajDev_NES, Late_NES),
    names_to = "stage_var",
    values_to = "NES"
  ) |>
  mutate(
    stage = sub("_NES", "", stage_var),
    stage = factor(stage, levels = c("Early", "TrajDev", "Late"))
  )

# Add corresponding padj values
dev_data_long <- dev_data_long |>
  mutate(
    padj = case_when(
      stage == "Early" ~ Early_padj,
      stage == "TrajDev" ~ TrajDev_padj,
      stage == "Late" ~ Late_padj,
      TRUE ~ NA_real_
    )
  )

# Set category order (for faceting)
category_order <- c("Neuronal", "Cytoplasmic Translation",
                   "Mitochondrial Translation", "Mitochondrial Metabolism",
                   "Calcium Signaling")
dev_data_long$category <- factor(dev_data_long$category, levels = category_order)

# Set pathway order within categories (alphabetical within category)
pathway_order <- dev_data |>
  distinct(pathway, category) |>
  arrange(category, pathway) |>
  pull(pathway)
dev_data_long$pathway <- factor(dev_data_long$pathway, levels = rev(pathway_order))

message("✓ Data reshaped: ", nrow(dev_data_long), " rows")

# -------------------------------------------------------------------- #
# 8. Save CSV Data Export                                              #
# -------------------------------------------------------------------- #

message("\n💾 Exporting CSV data table...")

csv_output <- dev_data |>
  select(pathway, database, category, mutation,
         Early_NES, Early_padj,
         TrajDev_NES, TrajDev_padj,
         Late_NES, Late_padj,
         TimeCtrl_NES, TimeCtrl_padj,
         SuccessScore, SuccessSymbol) |>
  arrange(category, pathway, mutation)

write_csv(csv_output, file.path(dev_framework_dir, "Developmental_Framework_Data.csv"))
message("✓ Saved CSV: Developmental_Framework_Data.csv")

# -------------------------------------------------------------------- #
# 9. Create Main Developmental Framework Dotplot                      #
# -------------------------------------------------------------------- #

message("\n🎨 Creating developmental framework dotplot...")

# Custom theme
theme_dev_framework <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "gray95", color = "gray70"),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10, face = "bold"),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30")
    )
}

# Create the plot
create_developmental_dotplot <- function(data, show_time_ctrl = TRUE) {
  # Calculate dot size from padj
  data <- data |>
    mutate(
      sig_size = ifelse(is.na(padj), 0, -log10(padj)),
      sig_size = pmin(sig_size, 5),  # Cap at 5 for visualization
      is_sig = !is.na(padj) & padj < 0.05
    )

  # Base plot
  p <- ggplot(data, aes(x = stage, y = pathway)) +
    # Add NES colored dots
    geom_point(aes(color = NES, size = sig_size), alpha = 0.9) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = c(-3.5, 3.5),
      na.value = "gray90",
      name = "NES"
    ) +
    scale_size_continuous(
      range = c(1, 8),
      limits = c(0, 5),
      breaks = c(0, 1.3, 2, 3),
      labels = c("NS", "0.05", "0.01", "0.001"),
      name = "FDR"
    ) +
    # Add black outline for significant hits
    geom_point(
      data = filter(data, is_sig),
      aes(x = stage, y = pathway, size = sig_size),
      shape = 21, color = "black", fill = NA, stroke = 0.8
    ) +
    facet_wrap(~ mutation, ncol = 1, scales = "free_y") +
    theme_dev_framework() +
    labs(
      x = NULL,
      y = NULL,
      title = "Developmental Framework: Early → TrajDev → Late",
      subtitle = "Mutation effects across developmental stages"
    )

  # Add Time_Ctrl markers in TrajDev column (small gray triangles)
  if (show_time_ctrl) {
    time_ctrl_data <- data |>
      filter(stage == "TrajDev", !is.na(TimeCtrl_NES)) |>
      distinct(pathway, mutation, TimeCtrl_NES)

    if (nrow(time_ctrl_data) > 0) {
      p <- p +
        geom_point(
          data = time_ctrl_data,
          aes(x = "TrajDev", y = pathway),
          color = "gray40",
          shape = 17,  # Triangle
          size = 2,
          alpha = 0.6
        )
    }
  }

  # Add SuccessScore annotations in Late column
  success_data <- data |>
    filter(stage == "Late", SuccessSymbol != "") |>
    distinct(pathway, mutation, SuccessSymbol)

  if (nrow(success_data) > 0) {
    p <- p +
      geom_text(
        data = success_data,
        aes(x = "Late", y = pathway, label = SuccessSymbol),
        nudge_x = 0.35,
        size = 4,
        fontface = "bold",
        color = "gray20"
      )
  }

  return(p)
}

# Generate main plot
main_plot <- create_developmental_dotplot(dev_data_long, show_time_ctrl = TRUE)

# Save main PDF
pdf(file.path(dev_framework_dir, "Cross_database_Developmental_Framework_MAIN.pdf"),
    width = 12, height = 10)
print(main_plot)
dev.off()

message("✓ Saved main figure: Cross_database_Developmental_Framework_MAIN.pdf")

# -------------------------------------------------------------------- #
# 10. Create Separate PDFs by Semantic Category                       #
# -------------------------------------------------------------------- #

message("\n📚 Creating category-specific figures...")

for (cat in category_order) {
  cat_data <- dev_data_long |> filter(category == cat)

  if (nrow(cat_data) == 0) {
    message("  ⚠️  No data for category: ", cat)
    next
  }

  cat_plot <- create_developmental_dotplot(cat_data, show_time_ctrl = TRUE) +
    labs(
      title = paste("Developmental Framework:", cat),
      subtitle = "Early → TrajDev → Late"
    )

  # Save category PDF
  cat_filename <- gsub(" ", "_", cat)
  pdf(file.path(dev_framework_dir, "By_Category", paste0(cat_filename, "_Developmental_Framework.pdf")),
      width = 10, height = 6)
  print(cat_plot)
  dev.off()

  message("  ✓ Saved: ", cat_filename, "_Developmental_Framework.pdf")
}

# -------------------------------------------------------------------- #
# 11. Create Interactive HTML Version with Plotly                     #
# -------------------------------------------------------------------- #

message("\n🌐 Creating interactive HTML version...")

# Prepare data for plotly with hover text
dev_data_plotly <- dev_data_long |>
  mutate(
    hover_text = paste0(
      "<b>", as.character(pathway), "</b><br>",
      "Database: ", database, "<br>",
      "Category: ", category, "<br>",
      "Mutation: ", mutation, "<br>",
      "Stage: ", as.character(stage), "<br>",
      "NES: ", round(NES, 2), "<br>",
      "FDR: ", ifelse(is.na(padj), "NA", format(padj, digits = 3, scientific = TRUE)),
      ifelse(as.character(stage) == "TrajDev" & !is.na(TimeCtrl_NES),
             paste0("<br>Time_Ctrl NES: ", round(TimeCtrl_NES, 2)),
             "")
    ),
    sig_size = ifelse(is.na(padj), 1, -log10(padj)),
    sig_size = pmin(sig_size, 5),
    sig_size = pmax(sig_size, 1),
    # Create stage numeric for x-axis
    stage_num = as.numeric(stage),
    # Create pathway numeric for y-axis
    pathway_num = as.numeric(pathway),
    # Marker border color
    border_color = ifelse(!is.na(padj) & padj < 0.05, "black", "transparent")
  )

# Create separate traces for each mutation to get proper faceting
fig_g32a <- dev_data_plotly |>
  filter(mutation == "G32A") |>
  plot_ly() |>
  add_trace(
    type = "scatter",
    mode = "markers",
    x = ~stage_num,
    y = ~pathway_num,
    color = ~NES,
    colors = c("blue", "white", "red"),
    cmin = -3.5,
    cmax = 3.5,
    size = ~sig_size,
    sizes = c(10, 100),
    marker = list(line = list(width = 1)),
    text = ~hover_text,
    hoverinfo = "text",
    name = "G32A",
    showlegend = FALSE
  )

fig_r403c <- dev_data_plotly |>
  filter(mutation == "R403C") |>
  plot_ly() |>
  add_trace(
    type = "scatter",
    mode = "markers",
    x = ~stage_num,
    y = ~pathway_num,
    color = ~NES,
    colors = c("blue", "white", "red"),
    cmin = -3.5,
    cmax = 3.5,
    size = ~sig_size,
    sizes = c(10, 100),
    marker = list(line = list(width = 1)),
    text = ~hover_text,
    hoverinfo = "text",
    name = "R403C",
    showlegend = FALSE
  )

# Combine into subplot
plotly_fig <- subplot(
  fig_g32a, fig_r403c,
  nrows = 2,
  shareX = TRUE,
  shareY = TRUE,
  titleY = TRUE,
  titleX = FALSE
) |>
  layout(
    title = list(text = "Developmental Framework: Early -> TrajDev -> Late (Interactive)"),
    xaxis = list(
      title = "",
      tickmode = "array",
      tickvals = c(1, 2, 3),
      ticktext = c("Early", "TrajDev", "Late")
    ),
    yaxis = list(
      title = "",
      tickmode = "array",
      tickvals = seq_along(pathway_order),
      ticktext = rev(pathway_order)
    ),
    hovermode = "closest"
  )

# Save HTML with error handling
tryCatch({
  htmlwidgets::saveWidget(
    plotly_fig,
    file.path(dev_framework_dir, "Developmental_Framework_Interactive.html"),
    selfcontained = TRUE
  )
  message("✓ Saved interactive HTML: Developmental_Framework_Interactive.html")
}, error = function(e) {
  message("⚠️  Could not create interactive HTML: ", e$message)
  message("   PDFs are available as main output")
})

# -------------------------------------------------------------------- #
# 12. Create Legend Guide PDF                                         #
# -------------------------------------------------------------------- #

message("\n📖 Creating legend guide...")

# Create a simple guide explaining annotations
pdf(file.path(dev_framework_dir, "Legend_Guide.pdf"), width = 8, height = 6)

# Create blank plot for text
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

# Title
text(0.5, 0.95, "Developmental Framework - Legend Guide", cex = 1.5, font = 2, adj = 0.5)

# Column explanations
text(0.1, 0.85, "Column Interpretation:", cex = 1.2, font = 2, adj = 0)
text(0.15, 0.80, "• Early: Baseline defect at D35 (mutant vs control)", cex = 1, adj = 0)
text(0.15, 0.75, "• TrajDev: Trajectory deviation (mutant maturation vs normal)", cex = 1, adj = 0)
text(0.15, 0.70, "• Late: Final outcome at D65 (mutant vs control)", cex = 1, adj = 0)

# Time_Ctrl marker
text(0.1, 0.60, "Time_Ctrl Marker (gray triangle ▼):", cex = 1.2, font = 2, adj = 0)
text(0.15, 0.55, "• Shows normal maturation direction in TrajDev column", cex = 1, adj = 0)
text(0.15, 0.50, "• Reference: How control cells change from D35 to D65", cex = 1, adj = 0)

# SuccessScore symbols
text(0.1, 0.40, "SuccessScore Annotations (right of Late column):", cex = 1.2, font = 2, adj = 0)
text(0.15, 0.35, "↑ = Improvement: |Late| < |Early| (defect reduced)", cex = 1, adj = 0, col = "darkgreen")
text(0.15, 0.30, "↓ = Worsening: |Late| > |Early| (defect increased)", cex = 1, adj = 0, col = "darkred")
text(0.15, 0.25, "− = Unchanged: |Late| ≈ |Early| (persistent defect)", cex = 1, adj = 0, col = "gray30")

# Significance
text(0.1, 0.15, "Significance Indicators:", cex = 1.2, font = 2, adj = 0)
text(0.15, 0.10, "• Dot size = significance (-log10 FDR)", cex = 1, adj = 0)
text(0.15, 0.05, "• Black outline = FDR < 0.05", cex = 1, adj = 0)

dev.off()

message("✓ Saved legend guide: Legend_Guide.pdf")

# -------------------------------------------------------------------- #
# 13. Summary Report                                                   #
# -------------------------------------------------------------------- #

message("\n", paste(rep("=", 78), collapse = ""))
message("✓ Developmental Framework Visualization Complete!")
message(paste(rep("=", 78), collapse = ""))

message("\n📊 Summary:")
message("  Pathways analyzed: ", length(pathway_list))
message("  Categories: ", length(category_order))
message("  Mutations: 2 (G32A, R403C)")
message("  Developmental stages: 3 (Early, TrajDev, Late)")

message("\n📁 Output files:")
message("  1. Cross_database_Developmental_Framework_MAIN.pdf (publication figure)")
message("  2. Developmental_Framework_Data.csv (data table)")
message("  3. Developmental_Framework_Interactive.html (interactive version)")
message("  4. Legend_Guide.pdf (annotation explanations)")
message("  5. By_Category/*.pdf (", length(category_order), " category-specific figures)")

message("\n💡 Key insights:")
# Count patterns
improvement_count <- sum(dev_data$SuccessSymbol == "↑", na.rm = TRUE)
worsening_count <- sum(dev_data$SuccessSymbol == "↓", na.rm = TRUE)
unchanged_count <- sum(dev_data$SuccessSymbol == "−", na.rm = TRUE)

message("  Pathway trajectories:")
message("    ↑ Improvement: ", improvement_count, " instances")
message("    ↓ Worsening:   ", worsening_count, " instances")
message("    − Unchanged:   ", unchanged_count, " instances")

message("\n✨ Done! Check output directory:")
message("   ", dev_framework_dir)
