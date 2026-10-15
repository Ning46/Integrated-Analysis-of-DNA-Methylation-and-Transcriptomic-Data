# FA versus NA gene expression, adjusted for sex by default.
source("config.R")
source("R/helpers.R")

annotate_rna <- function(table) {
  require_packages(c("AnnotationDbi", "org.Hs.eg.db"))
  table$Ensembl_ID <- sub("\\.[0-9]+$", "", rownames(table))
  db <- org.Hs.eg.db::org.Hs.eg.db
  # Stable per-key lookup: never reuse an index from another annotation table.
  map <- function(column) {
    x <- AnnotationDbi::mapIds(db, keys = unique(table$Ensembl_ID),
                              keytype = "ENSEMBL", column = column, multiVals = "first")
    unname(x[table$Ensembl_ID])
  }
  table$Gene_Symbol <- map("SYMBOL")
  table$Entrez_ID <- map("ENTREZID")
  table
}

rna_enrichment <- function(table, config, out) {
  # One tested Ensembl representative per Entrez ID, selected by expression,
  # not by P value. Universe and selection both use that same representative.
  table <- table[!is.na(table$Entrez_ID) & nzchar(table$Entrez_ID), , drop = FALSE]
  table <- table[order(-table$logCPM, table$Ensembl_ID), , drop = FALSE]
  table <- table[!duplicated(table$Entrez_ID), , drop = FALSE]
  universe <- as.character(table$Entrez_ID)
  if (length(universe) < 10L) {
    message("Too few mapped background genes for RNA enrichment.")
    return(invisible(NULL))
  }
  for (direction in c("up", "down")) {
    sig <- significant_rna(table, config) &
      if (direction == "up") table$logFC > 0 else table$logFC < 0
    selected <- universe[sig]
    if (!length(selected)) {
      message("No significant RNA genes for ", direction, " enrichment.")
      next
    }
    # covariate is aligned to the deduplicated tested Entrez universe.
    # Use mean log-expression to adjust detectability; no recycled length vector.
    for (collection in c("GO", "KEGG")) {
      res <- if (collection == "GO") {
        limma::goana(selected, universe = universe, species = "Hs",
                     covariate = table$logCPM)
      } else {
        limma::kegga(selected, universe = universe, species.KEGG = "hsa",
                     covariate = table$logCPM)
      }
      res$FDR <- stats::p.adjust(res$P.DE, method = "BH")
      res$Term_ID <- rownames(res)
      res <- res[order(res$FDR, res$P.DE), , drop = FALSE]
      write_tsv(res, file.path(out, paste(direction, collection, "enrichment.tsv", sep = "_")))
    }
  }
}

run_rnaseq <- function(config) {
  require_packages(c("edgeR", "limma", "AnnotationDbi", "org.Hs.eg.db"))
  metadata <- read_tsv(config$rna_samples)
  require_columns(metadata, c("sample_id", "group", "sex", "individual"), "RNA sample sheet")
  if (anyNA(metadata[, c("sample_id", "group", "sex", "individual")]) ||
      any(!nzchar(metadata$individual)) || !setequal(metadata$group, c("NA", "FA"))) {
    stop("RNA sample sheet requires nonmissing sample_id, group (NA/FA), sex, individual.")
  }
  # Never infer sample labels from directory order or hard-coded group vectors.
  metadata <- align_samples(metadata, metadata$sample_id)
  metadata$group <- factor(metadata$group, levels = c("NA", "FA"))
  metadata$sex <- factor(metadata$sex)
  if (anyDuplicated(metadata$individual)) {
    stop("Repeated RNA individuals detected. Select one biologically appropriate sample ",
         "per individual, or implement a design that handles within-person dependence.")
  }
  formula <- stats::as.formula(config$rna_design)
  variables <- all.vars(formula)
  require_columns(metadata, variables, "RNA design metadata")
  if (anyNA(metadata[, variables, drop = FALSE])) stop("Missing design covariates.")
  design <- stats::model.matrix(formula, metadata)
  rownames(design) <- metadata$sample_id
  check_design(design)
  if (!config$rna_coefficient %in% colnames(design)) stop("RNA coefficient not in design: ", config$rna_coefficient)

  if (config$rna_input == "kallisto") {
    require_packages("tximport")
    require_columns(metadata, "quant_file", "Kallisto sample sheet")
    if (anyNA(metadata$quant_file) || !all(file.exists(metadata$quant_file))) {
      stop("Every quant_file must name an existing Kallisto abundance TSV.")
    }
    tx2gene <- read_tsv(config$tx2gene_file)
    require_columns(tx2gene, c("TXNAME", "GENEID"), "Transcript-to-gene mapping")
    if (anyNA(tx2gene[, c("TXNAME", "GENEID")])) stop("Missing transcript/gene identifiers.")
    tx2gene <- unique(tx2gene[, c("TXNAME", "GENEID")])
    # Transcript versions/bars are removed consistently with tximport options.
    tx2gene$TXNAME <- sub("\\.[0-9]+$", "", sub("\\|.*$", "", tx2gene$TXNAME))
    tx2gene <- unique(tx2gene)
    if (anyDuplicated(tx2gene$TXNAME)) stop("A transcript maps to multiple genes.")
    files <- stats::setNames(metadata$quant_file, metadata$sample_id)
    imported <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene,
                                   ignoreAfterBar = TRUE, ignoreTxVersion = TRUE,
                                   countsFromAbundance = "lengthScaledTPM")
    counts <- imported$counts
  } else if (config$rna_input == "counts") {
    raw <- read_tsv(config$rna_counts_file)
    if (ncol(raw) < 2L) stop("Count table needs gene IDs and sample columns.")
    ids <- as.character(raw[[1]])
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop("Count gene IDs must be unique; resolve duplicated identifiers explicitly.")
    }
    if (!all(vapply(raw[, -1, drop = FALSE], is.numeric, logical(1)))) stop("Nonnumeric counts.")
    counts <- as.matrix(raw[, -1, drop = FALSE]); rownames(counts) <- ids
    metadata <- align_samples(metadata, colnames(counts))
    design <- stats::model.matrix(formula, metadata)
    rownames(design) <- metadata$sample_id
    check_design(design)
  } else stop("rna_input must be 'kallisto' or 'counts'.")
  if (!identical(colnames(counts), metadata$sample_id)) stop("RNA count/metadata order mismatch.")
  if (any(!is.finite(counts)) || any(counts < 0) || any(colSums(counts) <= 0)) {
    stop("RNA counts must be finite, nonnegative, with nonzero library sizes.")
  }
  if (anyDuplicated(sub("\\.[0-9]+$", "", rownames(counts)))) {
    stop("Gene IDs become duplicated after stripping Ensembl versions; resolve explicitly.")
  }
  y <- edgeR::DGEList(counts = counts, samples = metadata)
  keep <- edgeR::filterByExpr(y, design = design, min.count = 10)
  if (sum(keep) < 2L) stop("Too few genes after expression filtering.")
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, coef = config$rna_coefficient)
  table <- edgeR::topTags(test, n = Inf)$table
  table <- annotate_rna(table)
  table$significant <- significant_rna(table, config)
  out <- output_dir(config, "rnaseq")
  write_tsv(table, file.path(out, "res_edgeR.tsv"))
  write_tsv(metadata, file.path(out, "sample_metadata.tsv"))
  write_tsv(data.frame(sample_id = rownames(design), design, check.names = FALSE),
            file.path(out, "design_matrix.tsv"))
  grDevices::pdf(file.path(out, "rnaseq_qc.pdf"))
  graphics::boxplot(y$samples$lib.size ~ metadata$group, ylab = "Library size", xlab = "Group")
  limma::plotMDS(y, labels = metadata$sample_id, col = as.integer(metadata$group))
  edgeR::plotBCV(y)
  limma::plotMD(test)
  with(table, graphics::plot(logFC, -log10(pmax(PValue, .Machine$double.xmin)),
                             pch = 16, cex = 0.3, main = "RNA-seq volcano plot"))
  graphics::hist(table$PValue, main = "RNA-seq P values", xlab = "P value")
  cpms <- edgeR::cpm(y, log = TRUE)
  for (id in head(rownames(table), 20L)) {
    graphics::boxplot(as.numeric(cpms[id, ]) ~ metadata$group,
                      main = id, ylab = "Log2 CPM", xlab = "Group")
  }
  grDevices::dev.off()
  if (config$run_enrichment) rna_enrichment(table, config, out)
  record_run(config, out)
  invisible(table)
}
if (sys.nframe() == 0L) run_rnaseq(config)
