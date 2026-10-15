# Paired activation responses and FA-vs-NA response interaction.
source("config.R")
source("R/helpers.R")

run_methylation <- function(config) {
  require_packages(c("minfi", "ENmix", "limma", "GEOquery", "Biobase", "S4Vectors"))
  if (config$run_enrichment) require_packages("missMethyl")
  out <- output_dir(config, "methylation")
  metadata <- methylation_metadata(config)
  require_columns(metadata, c("sample_id", "group", "status", "individual"), "Metadata")
  if (nzchar(config$methylation_normalized_rds)) {
    if (!file.exists(config$methylation_normalized_rds)) stop("Normalized RDS not found.")
    normalized <- readRDS(config$methylation_normalized_rds)
    metadata <- align_samples(metadata, Biobase::sampleNames(normalized))
  } else {
    if (!dir.exists(config$idat_dir)) stop("IDAT directory not found: ", config$idat_dir)
    rg <- minfi::read.metharray.exp(config$idat_dir, extended = TRUE, force = TRUE)
    ids <- substr(Biobase::sampleNames(rg), 1, 10)
    metadata <- align_samples(metadata, ids)
    Biobase::sampleNames(rg) <- ids
    Biobase::pData(rg) <- S4Vectors::DataFrame(metadata)
    raw <- minfi::preprocessRaw(rg)
    grDevices::pdf(file.path(out, "raw_qc.pdf"))
    minfi::plotQC(minfi::getQC(raw))
    minfi::densityPlot(minfi::getBeta(raw), main = "Raw beta values")
    grDevices::dev.off()
    grDevices::pdf(file.path(out, "ENmix_qc.pdf"))
    qc <- ENmix::QCinfo(rg)
    grDevices::dev.off()
    capture.output(str(qc), file = file.path(out, "ENmix_QC_summary.txt"))
    meth <- minfi::getMeth(raw); unmeth <- minfi::getUnmeth(raw)
    invalid <- rownames(meth)[rowSums(!is.finite(meth) | !is.finite(unmeth)) > 0L]
    extra <- if (nzchar(config$excluded_cpg_file)) {
      scan(config$excluded_cpg_file, what = character(), quiet = TRUE)
    } else character()
    rm(raw, meth, unmeth)
    preprocessed <- ENmix::preprocessENmix(
      rg, bgParaEst = "oob", dyeCorr = "RELIC", QCinfo = qc,
      exQCsample = TRUE, exQCcpg = TRUE, exCpG = unique(c(invalid, extra)),
      nCores = config$n_cores
    )
    normalized <- ENmix::norm.quantile(preprocessed, method = "quantile2")
    retained <- Biobase::sampleNames(normalized)
    write_tsv(data.frame(sample_id = setdiff(metadata$sample_id, retained)),
              file.path(out, "qc_removed_samples.tsv"))
    metadata <- metadata[metadata$sample_id %in% retained, , drop = FALSE]
    metadata <- align_samples(metadata, retained)
    rm(rg, preprocessed)
  }
  Biobase::pData(normalized) <- S4Vectors::DataFrame(metadata)
  design <- paired_design(metadata) # fail on incomplete pairs after QC
  annotation <- as.data.frame(minfi::getAnnotation(normalized))
  if (!"Name" %in% names(annotation)) annotation$Name <- rownames(annotation)
  beta <- minfi::getMeth(normalized) /
    (minfi::getMeth(normalized) + minfi::getUnmeth(normalized) + 100)
  m <- log2((minfi::getMeth(normalized) + 1) / (minfi::getUnmeth(normalized) + 1))
  keep <- rowSums(!is.finite(m) | !is.finite(beta)) == 0L
  if (nzchar(config$excluded_cpg_file)) {
    excluded <- scan(config$excluded_cpg_file, what = character(), quiet = TRUE)
    keep <- keep & !rownames(m) %in% excluded
  }
  if (config$remove_sex_chromosomes) {
    require_columns(annotation, "chr", "Annotation")
    chr <- annotation$chr[match(rownames(m), annotation$Name)]
    keep <- keep & !is.na(chr) & !chr %in% c("chrX", "chrY", "X", "Y")
  }
  m <- m[keep, , drop = FALSE]; beta <- beta[keep, , drop = FALSE]
  if (nrow(m) < 2L) stop("Too few probes after filtering.")
  write_tsv(metadata, file.path(out, "sample_metadata.tsv"))
  write_tsv(data.frame(sample_id = rownames(design), design, check.names = FALSE),
            file.path(out, "design_matrix.tsv"))
  grDevices::pdf(file.path(out, "normalized_qc.pdf"))
  minfi::densityPlot(beta, main = "Normalized beta values")
  limma::plotMDS(m, labels = metadata$individual,
                 col = ifelse(metadata$group == "FA", "firebrick", "steelblue"))
  grDevices::dev.off()
  effects <- beta_effects(beta, metadata)
  contrasts <- matrix(0, ncol(design), 3,
                       dimnames = list(colnames(design), c("FA", "NA", "interaction")))
  contrasts["activated", c("FA", "NA")] <- 1
  contrasts["FA_activation", c("FA", "interaction")] <- 1
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(m, design), contrasts))
  tables <- genes <- list()
  for (contrast in colnames(contrasts)) {
    tab <- limma::topTable(fit, coef = contrast, number = Inf, adjust.method = "BH")
    tab <- attach_effects(tab, effects, contrast)
    tab$significant <- significant_methylation(tab, config)
    tables[[contrast]] <- tab
    write_tsv(tab, file.path(out, paste0("DMP_", contrast, ".tsv")))
    genes[[contrast]] <- expand_probe_annotation(tab, annotation)
    write_tsv(genes[[contrast]], file.path(out, paste0("DMP_", contrast, "_genes.tsv")))
    write_tsv(summarize_methylation_genes(genes[[contrast]][genes[[contrast]]$significant, ]),
              file.path(out, paste0("significant_", contrast, "_gene_summary.tsv")))
    grDevices::pdf(file.path(out, paste0("diagnostics_", contrast, ".pdf")))
    limma::volcanoplot(fit, coef = contrast)
    limma::plotMA(fit, coef = contrast)
    graphics::hist(tab$P.Value, main = paste(contrast, "P values"), xlab = "P value")
    grDevices::dev.off()
  }
  fa_genes <- clean_symbols(genes$FA$Gene_Symbol[genes$FA$significant])
  na_genes <- clean_symbols(genes[["NA"]]$Gene_Symbol[genes[["NA"]]$significant])
  write_tsv(data.frame(Gene_Symbol = setdiff(fa_genes, na_genes)),
            file.path(out, "exploratory_FA_only_gene_list.tsv"))
  # Direct interaction results are primary. A difference of significance lists
  # is retained as exploratory only; it is not a test of group differences.
  if (config$run_enrichment) {
    for (contrast in c("FA", "interaction")) {
      tab <- tables[[contrast]]
      for (direction in c("positive", "negative")) {
        sig <- tab$significant & if (direction == "positive") tab$deltaBeta > 0 else tab$deltaBeta < 0
        if (!any(sig)) {
          message("No significant CpGs for enrichment: ", contrast, "/", direction)
          next
        }
        for (collection in c("GO", "KEGG")) {
          enrichment <- missMethyl::gometh(
            sig.cpg = tab$CpG[sig], all.cpg = tab$CpG,
            collection = collection, array.type = config$methylation_array,
            anno = minfi::getAnnotation(normalized), prior.prob = TRUE
          )
          enrichment$Term_ID <- rownames(enrichment)
          write_tsv(enrichment, file.path(out, paste(contrast, direction, collection,
                                                       "enrichment.tsv", sep = "_")))
        }
      }
    }
  }
  record_run(config, out)
  invisible(tables)
}
if (sys.nframe() == 0L) run_methylation(config)
