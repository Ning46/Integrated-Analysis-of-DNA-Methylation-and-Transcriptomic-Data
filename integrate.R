# Gene-level result overlap, not sample-level correlation or causal inference.
source("config.R")
source("R/helpers.R")

integration_sets <- function(methylation, rna, config) {
  require_columns(methylation, c("CpG", "Gene_Symbol", "adj.P.Val", "deltaBeta"),
                  "Methylation results")
  require_columns(rna, c("Gene_Symbol", "FDR", "logFC"), "RNA results")
  background <- intersect(clean_symbols(methylation$Gene_Symbol), clean_symbols(rna$Gene_Symbol))
  methylated <- intersect(background, clean_symbols(methylation$Gene_Symbol[
    significant_methylation(methylation, config)]))
  expressed <- intersect(background, clean_symbols(rna$Gene_Symbol[significant_rna(rna, config)]))
  list(background = background, methylation = methylated, RNA = expressed,
       overlap = intersect(methylated, expressed))
}

run_integration <- function(config) {
  rna <- read_tsv(file.path(config$results_dir, "rnaseq", "res_edgeR.tsv"))
  out <- output_dir(config, "integration")
  for (contrast in c("FA", "interaction")) {
    methylation <- read_tsv(file.path(config$results_dir, "methylation",
                                      paste0("DMP_", contrast, "_genes.tsv")))
    sets <- integration_sets(methylation, rna, config)
    if (!length(sets$background)) stop("No common tested gene symbols; check annotation.")
    prefix <- file.path(out, contrast)
    write_tsv(data.frame(Gene_Symbol = sets$background), paste0(prefix, "_shared_tested_genes.tsv"))
    write_tsv(data.frame(Gene_Symbol = sets$overlap), paste0(prefix, "_overlap_genes.tsv"))
    methyl_sig <- methylation[significant_methylation(methylation, config), , drop = FALSE]
    rna_sig <- rna[significant_rna(rna, config), , drop = FALSE]
    m <- methyl_sig[methyl_sig$Gene_Symbol %in% sets$overlap, , drop = FALSE]
    r <- rna_sig[rna_sig$Gene_Symbol %in% sets$overlap, , drop = FALSE]
    # Join explicitly by symbol; preserve many CpGs and multiple Ensembl IDs.
    joined <- merge(m, r, by = "Gene_Symbol", suffixes = c("_methylation", "_RNA"))
    write_tsv(joined, paste0(prefix, "_overlap_details.tsv"))
    write_tsv(summarize_methylation_genes(m), paste0(prefix, "_methylation_gene_summary.tsv"))
    n <- length(sets$background); a <- length(sets$overlap)
    b <- length(setdiff(sets$methylation, sets$RNA))
    c <- length(setdiff(sets$RNA, sets$methylation)); d <- n - a - b - c
    contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                           dimnames = list(Methylation = c("yes", "no"), RNA = c("yes", "no")))
    # Descriptive only: CpGs/genes and annotation are dependent.
    p <- if (min(rowSums(contingency)) > 0 && min(colSums(contingency)) > 0) {
      stats::fisher.test(contingency, alternative = "greater")$p.value
    } else NA_real_
    write_tsv(data.frame(contrast = contrast, background = n,
                         methylation = length(sets$methylation), RNA = length(sets$RNA),
                         overlap = a, descriptive_fisher_p = p), paste0(prefix, "_overlap_summary.tsv"))
    grDevices::pdf(paste0(prefix, "_overlap_counts.pdf"))
    graphics::barplot(c(Methylation_only = b, Both = a, RNA_only = c),
                      main = paste(contrast, "shared-universe gene overlap"), ylab = "Unique genes")
    grDevices::dev.off()
  }
  record_run(config, out)
  message("Integration complete. This compares gene lists from different biological contrasts.")
}
if (sys.nframe() == 0L) run_integration(config)
