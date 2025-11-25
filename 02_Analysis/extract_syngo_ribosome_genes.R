###############################################################################
##  Extract SynGO Ribosome Gene Lists                                       ##
##  Script to extract presynaptic and postsynaptic ribosome genes           ##
###############################################################################

library(here)
library(dplyr)

# Load checkpoints
checkpoint_dir <- here("03_Results/02_Analysis/checkpoints")
syngo_lists <- readRDS(file.path(checkpoint_dir, "syngo_lists.rds"))
syngo_gsea_results <- readRDS(file.path(checkpoint_dir, "syngo_gsea_results.rds"))

# Output directory
out_dir <- here("03_Results/02_Analysis/Verification_reports")

# Extract gene lists from T2G (Term to Gene mapping)
T2G <- syngo_lists$T2G

message("Extracting SynGO ribosome gene lists...\n")

## Presynaptic ribosome genes
presyn_genes <- T2G %>%
  filter(gs_name == "SYNGO:presyn_ribosome") %>%
  pull(gene_symbol) %>%
  unique() %>%
  sort()

message(sprintf("✓ Presynaptic ribosome genes: %d", length(presyn_genes)))
writeLines(presyn_genes, file.path(out_dir, "syngo_presyn_ribosome_genes.txt"))

## Postsynaptic ribosome genes
postsyn_genes <- T2G %>%
  filter(gs_name == "SYNGO:postsyn_ribosome") %>%
  pull(gene_symbol) %>%
  unique() %>%
  sort()

message(sprintf("✓ Postsynaptic ribosome genes: %d", length(postsyn_genes)))
writeLines(postsyn_genes, file.path(out_dir, "syngo_postsyn_ribosome_genes.txt"))

## Combined (unique)
all_synaptic_ribosome_genes <- unique(c(presyn_genes, postsyn_genes)) %>% sort()
message(sprintf("✓ Total unique synaptic ribosome genes: %d", length(all_synaptic_ribosome_genes)))
writeLines(all_synaptic_ribosome_genes, file.path(out_dir, "syngo_all_synaptic_ribosome_genes.txt"))

## Overlap analysis
shared_genes <- intersect(presyn_genes, postsyn_genes)
presyn_only <- setdiff(presyn_genes, postsyn_genes)
postsyn_only <- setdiff(postsyn_genes, presyn_genes)

message(sprintf("\nOverlap analysis:"))
message(sprintf("  Shared between pre and post: %d genes", length(shared_genes)))
message(sprintf("  Presynaptic only: %d genes", length(presyn_only)))
message(sprintf("  Postsynaptic only: %d genes", length(postsyn_only)))

## Create summary data frame
gene_membership <- data.frame(
  Gene = all_synaptic_ribosome_genes,
  Presynaptic = all_synaptic_ribosome_genes %in% presyn_genes,
  Postsynaptic = all_synaptic_ribosome_genes %in% postsyn_genes,
  Both = all_synaptic_ribosome_genes %in% shared_genes,
  stringsAsFactors = FALSE
)

write.csv(gene_membership,
          file.path(out_dir, "syngo_ribosome_gene_membership.csv"),
          row.names = FALSE)

## Extract core enrichment genes from maturation contrasts
message("\n\nExtracting core enrichment (leading edge) genes from maturation contrasts...")

for (contrast in c("Maturation_G32A_specific", "Maturation_R403C_specific")) {
  message(sprintf("\n%s:", contrast))

  syngo_res <- syngo_gsea_results[[contrast]]
  if (!is.null(syngo_res)) {
    res_df <- syngo_res@result

    # Presynaptic
    presyn_row <- res_df[res_df$ID == "SYNGO:presyn_ribosome", ]
    if (nrow(presyn_row) > 0) {
      presyn_core <- unlist(strsplit(presyn_row$core_enrichment, "/"))
      message(sprintf("  Presynaptic core enrichment: %d genes", length(presyn_core)))
      writeLines(presyn_core, file.path(out_dir, sprintf("%s_presyn_core_genes.txt", contrast)))
    }

    # Postsynaptic
    postsyn_row <- res_df[res_df$ID == "SYNGO:postsyn_ribosome", ]
    if (nrow(postsyn_row) > 0) {
      postsyn_core <- unlist(strsplit(postsyn_row$core_enrichment, "/"))
      message(sprintf("  Postsynaptic core enrichment: %d genes", length(postsyn_core)))
      writeLines(postsyn_core, file.path(out_dir, sprintf("%s_postsyn_core_genes.txt", contrast)))
    }
  }
}

message("\n\n✓ Gene list extraction complete!")
message(sprintf("Output directory: %s", out_dir))
