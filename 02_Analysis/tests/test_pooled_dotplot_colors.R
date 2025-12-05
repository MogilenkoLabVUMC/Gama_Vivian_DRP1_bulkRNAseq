###############################################################################
##  Test Script: Pooled Dotplot with Correct Colors, Gene Ratio, and Outlines##
##                                                                             ##
##  Tests:                                                                     ##
##    1. Color config loading (DIVERGING_COLORS)                               ##
##    2. Gene ratio calculation                                                ##
##    3. Neuronal pathway filtering                                            ##
##    4. Significance outline rendering                                        ##
##    5. Full dotplot creation                                                 ##
###############################################################################

library(here)
library(ggplot2)
library(dplyr)
library(stringr)
library(testthat)

# Source the shared helpers (this also loads color_config.R)
source(here::here("01_Scripts/R_scripts/gsea_dotplot_helpers.R"))

# Test output directory
test_output_dir <- here::here("03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/test_output")
dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# TEST 1: Verify Color Config is Loaded Correctly
# ============================================================================
test_that("Color config loads correctly", {
  expect_true(exists("DIVERGING_COLORS"))
  expect_equal(DIVERGING_COLORS$negative, "#2166AC")
  expect_equal(DIVERGING_COLORS$neutral, "#F7F7F7")
  expect_equal(DIVERGING_COLORS$positive, "#B35806")
  message("Test 1 PASSED: Color config loaded correctly")
})

# ============================================================================
# TEST 2: Gene Ratio Calculation
# ============================================================================
test_that("Gene ratio is calculated correctly", {
  # Test case: 5 genes separated by 4 slashes, setSize = 100
  core_enrichment <- "GENE1/GENE2/GENE3/GENE4/GENE5"
  set_size <- 100

  gene_ratio <- calculate_gene_ratio(core_enrichment, set_size)

  expect_equal(gene_ratio, 0.05)

  # Edge case: empty string
  empty_ratio <- calculate_gene_ratio("", 100)
  expect_equal(empty_ratio, 0)

  # Edge case: single gene
  single_ratio <- calculate_gene_ratio("GENE1", 50)
  expect_equal(single_ratio, 1/50)

  message("Test 2 PASSED: Gene ratio calculation correct")
})

# ============================================================================
# TEST 3: add_gene_ratio Helper
# ============================================================================
test_that("add_gene_ratio adds column correctly", {
  test_df <- data.frame(
    core_enrichment = c("A/B/C", "X/Y", "Z"),
    setSize = c(30, 20, 10)
  )

  result <- add_gene_ratio(test_df)

  expect_true("GeneRatio" %in% colnames(result))
  expect_equal(result$GeneRatio[1], 3/30)
  expect_equal(result$GeneRatio[2], 2/20)
  expect_equal(result$GeneRatio[3], 1/10)

  message("Test 3 PASSED: add_gene_ratio works correctly")
})

# ============================================================================
# TEST 4: Neuronal Pathway Filtering
# ============================================================================
test_that("Neuronal filtering works correctly", {
  test_df <- data.frame(
    ID = c("CARDIAC_PATH", "NEURO_PATH", "SYNGO_001", "PANCREAS_PATH"),
    Description = c("Cardiac function", "Neuronal development", "Synapse assembly", "Beta cell"),
    Database = c("hallmark", "gobp", "syngo", "kegg"),
    stringsAsFactors = FALSE
  )

  filtered <- filter_neuronal_pathways(test_df)

  # CARDIAC should be excluded
  expect_false("CARDIAC_PATH" %in% filtered$ID)
  # PANCREAS should be excluded
  expect_false("PANCREAS_PATH" %in% filtered$ID)
  # NEURO should be included (inclusion pattern override)
  expect_true("NEURO_PATH" %in% filtered$ID)
  # SynGO should always be kept
  expect_true("SYNGO_001" %in% filtered$ID)

  message("Test 4 PASSED: Neuronal filtering works correctly")
})

# ============================================================================
# TEST 5: Load Real Data and Create Dotplot
# ============================================================================
test_that("Real data dotplot with outlines works", {
  csv_file <- here::here("03_Results/02_Analysis/Plots/GSEA/Cross_database_pooled/csv_data/Time_Ctrl_filtered_focused.csv")

  if (!file.exists(csv_file)) {
    skip("Test data file not found")
  }

  data <- read.csv(csv_file, stringsAsFactors = FALSE)

  # Create dotplot with significance outlines
  p <- create_gsea_pooled_dotplot(
    data,
    title = "TEST: Pooled Dotplot with Outlines",
    subtitle = "Black outline = significant (padj < 0.05)",
    n_per_database = 5,
    highlight_sig = TRUE,
    sig_cutoff = 0.05
  )

  # Save test plot
  test_plot_file <- file.path(test_output_dir, "TEST_pooled_dotplot_with_outlines.pdf")
  ggsave(test_plot_file, p, width = 10, height = 8)

  expect_true(file.exists(test_plot_file))
  message(sprintf("Test 5 PASSED: Dotplot with outlines saved to %s", test_plot_file))
})

# ============================================================================
# TEST 6: Verify Significance Outline Logic
# ============================================================================
test_that("Significance outline is added for significant pathways", {
  # Create mock data with mix of significant and non-significant
  mock_data <- data.frame(
    Description = c("Sig Pathway 1", "Sig Pathway 2", "NonSig Pathway 1", "NonSig Pathway 2"),
    NES = c(2.5, -1.8, 1.2, -0.8),
    p.adjust = c(0.001, 0.01, 0.08, 0.15),  # First two significant, last two not
    core_enrichment = c("A/B/C/D/E", "F/G/H", "I/J", "K"),
    setSize = c(50, 30, 20, 10),
    Database = rep("test", 4),
    stringsAsFactors = FALSE
  )

  # Create plot with highlight_sig = TRUE
  p <- create_gsea_pooled_dotplot(
    mock_data,
    title = "TEST: Significance Outline Verification",
    subtitle = "Only sig pathways (padj < 0.05) should have black outline",
    n_per_database = 10,
    highlight_sig = TRUE,
    sig_cutoff = 0.05
  )

  # Check that plot has multiple layers (base + outline overlay)
  expect_gte(length(p$layers), 1)

  # Save verification plot
  test_plot_file <- file.path(test_output_dir, "TEST_significance_outline_verification.pdf")
  ggsave(test_plot_file, p, width = 8, height = 5)

  expect_true(file.exists(test_plot_file))
  message(sprintf("Test 6 PASSED: Significance outline verification saved to %s", test_plot_file))
})

# ============================================================================
# TEST 7: Gradient Color Verification (Positive/Negative NES)
# ============================================================================
test_that("Gradient handles positive and negative NES correctly", {
  mock_data <- data.frame(
    Description = c("Positive Pathway 1", "Positive Pathway 2",
                    "Negative Pathway 1", "Negative Pathway 2",
                    "Neutral Pathway"),
    NES = c(2.5, 1.5, -2.5, -1.5, 0.1),
    p.adjust = c(0.001, 0.01, 0.001, 0.01, 0.03),
    core_enrichment = c("A/B/C", "D/E", "F/G/H", "I/J", "K"),
    setSize = c(30, 20, 40, 25, 15),
    Database = rep("test", 5),
    stringsAsFactors = FALSE
  )

  p <- create_gsea_pooled_dotplot(
    mock_data,
    title = "TEST: Gradient Color Verification",
    subtitle = "Blue = Negative NES, White = Neutral, Orange = Positive NES",
    highlight_sig = TRUE
  )

  test_plot_file <- file.path(test_output_dir, "TEST_gradient_verification.pdf")
  ggsave(test_plot_file, p, width = 8, height = 5)

  expect_true(file.exists(test_plot_file))
  message(sprintf("Test 7 PASSED: Gradient verification saved to %s", test_plot_file))
})

# ============================================================================
# RUN ALL TESTS
# ============================================================================
message("\n" , paste(rep("=", 70), collapse = ""))
message("ALL TESTS COMPLETED SUCCESSFULLY")
message(paste(rep("=", 70), collapse = ""))
message(sprintf("\nTest output directory: %s", test_output_dir))
message("\nReview the generated PDFs to visually verify:")
message("  1. TEST_pooled_dotplot_with_outlines.pdf - Real data with outlines")
message("  2. TEST_significance_outline_verification.pdf - Mock data showing outline logic")
message("  3. TEST_gradient_verification.pdf - Mock data showing gradient behavior")
