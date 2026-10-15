# Offline regression checks using synthetic data only. Run from project root.
source("all_for_DNA.R")
source("all_for_rnaseq.R")
source("integrate.R")
set.seed(2026)
expect_error <- function(expr) {
  ok <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  stopifnot(ok)
}

# Metadata must follow explicit identifiers; preserve the literal "NA" group.
meta <- data.frame(sample_id = paste0("s", 1:12),
                   individual = rep(paste0("p", 1:6), each = 2),
                   group = rep(c("NA", "FA"), each = 6),
                   status = rep(c("U", "S"), 6), sex = rep(c("F", "M"), 6))
stopifnot(identical(align_samples(meta[12:1, ], meta$sample_id)$sample_id, meta$sample_id))
expect_error(align_samples(meta[-1, ], meta$sample_id))
expect_error(paired_design(meta[-1, ]))
design <- paired_design(meta)

# A synthetic FA-specific activation effect has a positive interaction.
m <- matrix(rnorm(100 * 12, sd = 0.1), 100, 12,
             dimnames = list(paste0("cg", 1:100), meta$sample_id))
m[1, ] <- m[1, ] + 2 * (meta$status == "S" & meta$group == "FA")
fit <- limma::eBayes(limma::lmFit(m, design))
stopifnot(fit$coefficients[1, "FA_activation"] > 1.5)
beta <- plogis(m)
effects <- beta_effects(beta, meta)
tab <- data.frame(logFC = c(1, 2), row.names = c("cg2", "cg1"))
attached <- attach_effects(tab, effects, "FA")
stopifnot(identical(attached$CpG, c("cg2", "cg1")),
          isTRUE(all.equal(attached$deltaBeta, effects$deltaBeta_FA[c(2, 1)])))

# Annotation expands genes and paired features, never the scalar island field.
ann <- data.frame(Name = c("cg1", "cg2"), UCSC_RefGene_Name = c("A;B", "C"),
                   UCSC_RefGene_Group = c("TSS200;Body", "Body"),
                   Relation_to_Island = c("Island", "OpenSea"))
expanded <- expand_probe_annotation(attached, ann)
stopifnot(nrow(expanded) == 3L,
          expanded$CpG_Feature[expanded$Gene_Symbol == "B"] == "Body")
test_config <- config
test_config$run_enrichment <- FALSE
test_config$fdr <- 0.05
test_config$min_delta_beta <- 0.05
test_config$min_abs_rna_logfc <- 0.5
ms <- data.frame(CpG = c("a", "b", "c"), Gene_Symbol = c("A", "B", "C"),
                  adj.P.Val = c(0.001, 0.8, 0.001), deltaBeta = c(0.1, 0.1, 0.01))
rs <- data.frame(Gene_Symbol = c("A", "B", "D"), FDR = c(0.001, 0.8, 0.001),
                  logFC = c(1, 1, 1))
sets <- integration_sets(ms, rs, test_config)
stopifnot(identical(sets$background, c("A", "B")), identical(sets$overlap, "A"))

# End-to-end RNA and integration smoke test with synthetic counts and local
# annotation. All synthetic inputs/results go to a temporary directory.
require_packages(c("edgeR", "limma", "org.Hs.eg.db", "AnnotationDbi"))
tmp <- tempfile("omics_workflow_")
dir.create(tmp)
test_config$results_dir <- file.path(tmp, "results")
test_config$rna_input <- "counts"
test_config$rna_samples <- file.path(tmp, "samples.tsv")
test_config$rna_counts_file <- file.path(tmp, "counts.tsv")
rna_meta <- data.frame(sample_id = paste0("RNA", 1:12),
                       group = rep(c("NA", "FA"), each = 6),
                       sex = rep(c("F", "M"), 6), individual = paste0("person", 1:12))
gene_ids <- head(AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL"), 300)
counts <- matrix(rnbinom(300 * 12, mu = 100, size = 10), 300, 12,
                  dimnames = list(gene_ids, rna_meta$sample_id))
counts[1:20, 7:12] <- counts[1:20, 7:12] * 6
write_tsv(rna_meta[12:1, ], test_config$rna_samples) # deliberately scrambled metadata
write_tsv(data.frame(Gene = gene_ids, counts, check.names = FALSE), test_config$rna_counts_file)
rna_result <- run_rnaseq(test_config)
stopifnot(nrow(rna_result) > 0, all(c("Gene_Symbol", "Ensembl_ID", "FDR") %in% names(rna_result)))
saved_meta <- read_tsv(file.path(test_config$results_dir, "rnaseq", "sample_metadata.tsv"))
stopifnot(identical(saved_meta$sample_id, colnames(counts)),
          sum(saved_meta$group == "NA") == 6L)

# Exercise the complete methylation model/export route on synthetic normalized
# intensities with real platform annotation (no IDAT or study data required).
require_packages(c("minfi", "ENmix", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))
anno <- minfi::getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
probe_ids <- head(anno$Name[!is.na(anno$UCSC_RefGene_Name) & nzchar(anno$UCSC_RefGene_Name)], 300)
intensity <- matrix(runif(300 * 12, 1500, 8000), 300, 12,
                     dimnames = list(probe_ids, meta$sample_id))
unmeth <- matrix(runif(300 * 12, 1500, 8000), 300, 12, dimnames = dimnames(intensity))
intensity[1:20, meta$group == "FA" & meta$status == "S"] <- 20000
normalized <- minfi::MethylSet(
  Meth = intensity, Unmeth = unmeth,
  annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
)
test_config$methylation_normalized_rds <- file.path(tmp, "synthetic_normalized.rds")
test_config$methylation_samples <- file.path(tmp, "methylation_samples.tsv")
saveRDS(normalized, test_config$methylation_normalized_rds)
write_tsv(meta[12:1, ], test_config$methylation_samples)
methyl_results <- run_methylation(test_config)
stopifnot(all(c("FA", "NA", "interaction") %in% names(methyl_results)),
          all(c("CpG", "deltaBeta", "adj.P.Val") %in% names(methyl_results$interaction)))
symbol <- clean_symbols(rna_result$Gene_Symbol)[1]
for (contrast in c("FA", "interaction")) {
  fake_methyl <- data.frame(CpG = "synthetic_probe", Gene_Symbol = symbol,
                            adj.P.Val = 0.001, deltaBeta = 0.1)
  write_tsv(fake_methyl, file.path(test_config$results_dir, "methylation",
                                  paste0("DMP_", contrast, "_genes.tsv")))
}
run_integration(test_config)
stopifnot(file.exists(file.path(test_config$results_dir, "integration", "FA_overlap_details.tsv")))
cat("All offline synthetic workflow checks passed.\n")
