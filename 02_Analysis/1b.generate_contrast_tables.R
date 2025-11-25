###############################################################################
##  Generate contrast_tables checkpoint from existing model objects          ##
###############################################################################
##  Purpose: Create contrast_tables.rds checkpoint for downstream modules
##           without re-running the entire main pipeline
###############################################################################

library(here)

# Configuration
checkpoint_dir <- here::here("03_Results/02_Analysis/checkpoints")
helper_root <- "01_Scripts/RNAseq-toolkit"

# Source helper for process_all_contrasts
source(here::here(helper_root, "scripts/utils_plotting.R"))

# Load existing model objects
message("📂 Loading model objects from checkpoint...")
model_objs <- readRDS(file.path(checkpoint_dir, "model_objects.rds"))

fit <- model_objs$fit
contrasts <- model_objs$contrasts

# Process all contrasts
message("🔬 Processing all contrasts...")
contrast_tables <- process_all_contrasts(fit, contrasts, number = Inf, sort_by = "t")

# Save checkpoint
message("💾 Saving contrast_tables checkpoint...")
saveRDS(contrast_tables, file.path(checkpoint_dir, "contrast_tables.rds"))

message("✓ contrast_tables.rds created successfully")
message("  Location: ", file.path(checkpoint_dir, "contrast_tables.rds"))
message("  Contrasts: ", paste(names(contrast_tables), collapse = ", "))
