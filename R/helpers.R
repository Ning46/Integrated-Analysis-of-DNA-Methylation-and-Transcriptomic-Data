require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "),
                            ". See README.md.", call. = FALSE)
}
require_columns <- function(x, columns, context) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(context, " missing columns: ", paste(missing, collapse = ", "))
}
read_tsv <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                   na.strings = "", quote = "", comment.char = "")
}
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
output_dir <- function(config, stage) {
  path <- file.path(config$results_dir, stage)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
record_run <- function(config, out) {
  dput(config, file = file.path(out, "run_config.R"))
  capture.output(utils::sessionInfo(), file = file.path(out, "sessionInfo.txt"))
}
align_samples <- function(metadata, sample_ids) {
  require_columns(metadata, "sample_id", "Sample metadata")
  if (anyNA(metadata$sample_id) || any(!nzchar(metadata$sample_id)) ||
      anyDuplicated(metadata$sample_id) || anyDuplicated(sample_ids)) {
    stop("Sample IDs must be unique and nonempty.")
  }
  if (!setequal(metadata$sample_id, sample_ids)) {
    stop("Sample identifiers do not match: missing metadata for [",
         paste(setdiff(sample_ids, metadata$sample_id), collapse = ", "),
         "]; unused metadata rows [",
         paste(setdiff(metadata$sample_id, sample_ids), collapse = ", "), "].")
  }
  metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]
  rownames(metadata) <- metadata$sample_id
  metadata
}
check_design <- function(design) {
  if (anyNA(design) || qr(design)$rank < ncol(design)) {
    stop("Design matrix is rank deficient; check confounding and factor levels.")
  }
  if (nrow(design) <= ncol(design)) stop("No residual degrees of freedom in design.")
  invisible(design)
}
methylation_metadata <- function(config) {
  if (nzchar(config$methylation_samples)) return(read_tsv(config$methylation_samples))
  dir.create(config$geo_cache, recursive = TRUE, showWarnings = FALSE)
  geo <- GEOquery::getGEO("GSE189148", destdir = config$geo_cache,
                         GSEMatrix = TRUE)[[1]]
  x <- Biobase::pData(geo)
  require_columns(x, c("indiidual:ch1", "stimulation_status:ch1", "Sex:ch1"),
                  "GSE189148 metadata")
  individual <- as.character(x[["indiidual:ch1"]])
  if (anyNA(individual) || any(!nzchar(individual))) stop("Missing individual IDs.")
  data.frame(sample_id = rownames(x),
             group = ifelse(grepl("NA", individual), "NA", "FA"),
             status = x[["stimulation_status:ch1"]], individual = individual,
             sex = x[["Sex:ch1"]], stringsAsFactors = FALSE)
}
paired_design <- function(metadata) {
  require_columns(metadata, c("sample_id", "group", "status", "individual"),
                  "Methylation metadata")
  if (anyNA(metadata[, c("group", "status", "individual")]) ||
      any(!nzchar(metadata$individual)) ||
      !setequal(metadata$group, c("NA", "FA")) ||
      !setequal(metadata$status, c("U", "S"))) {
    stop("Expected nonmissing individual IDs, groups NA/FA, and statuses U/S.")
  }
  by_subject <- split(metadata, metadata$individual)
  valid <- vapply(by_subject, function(x) {
    nrow(x) == 2L && setequal(x$status, c("U", "S")) && length(unique(x$group)) == 1L
  }, logical(1))
  if (!all(valid)) stop("Each individual must have exactly one U/S pair in one group. Invalid: ",
                        paste(names(valid)[!valid], collapse = ", "),
                        ". Resolve incomplete pairs explicitly after QC.")
  if (any(table(vapply(by_subject, function(x) x$group[1], character(1))) < 2L)) {
    stop("At least two complete individuals per group are required.")
  }
  individual <- factor(metadata$individual)
  activated <- as.numeric(metadata$status == "S")
  FA_activation <- activated * as.numeric(metadata$group == "FA")
  design <- stats::model.matrix(~ 0 + individual + activated + FA_activation)
  rownames(design) <- metadata$sample_id
  check_design(design)
  design
}
beta_effects <- function(beta, metadata) {
  metadata <- align_samples(metadata, colnames(beta))
  means <- function(g, s) rowMeans(beta[, metadata$group == g & metadata$status == s,
                                        drop = FALSE])
  fa_u <- means("FA", "U"); fa_s <- means("FA", "S")
  na_u <- means("NA", "U"); na_s <- means("NA", "S")
  data.frame(CpG = rownames(beta), FA_U = fa_u, FA_S = fa_s, NA_U = na_u, NA_S = na_s,
             deltaBeta_FA = fa_s - fa_u, deltaBeta_NA = na_s - na_u,
             deltaBeta_interaction = (fa_s - fa_u) - (na_s - na_u),
             row.names = NULL)
}
attach_effects <- function(table, effects, contrast) {
  table$CpG <- rownames(table)
  ix <- match(table$CpG, effects$CpG)
  if (anyNA(ix)) stop("CpG effect lookup failed.")
  table <- cbind(table, effects[ix, setdiff(names(effects), "CpG"), drop = FALSE])
  table$deltaBeta <- table[[paste0("deltaBeta_", contrast)]]
  rownames(table) <- NULL
  table
}
expand_probe_annotation <- function(table, annotation) {
  require_columns(table, "CpG", "Probe results")
  require_columns(annotation, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group"),
                  "Probe annotation")
  if (anyDuplicated(annotation$Name)) stop("Duplicate CpG annotation IDs.")
  ix <- match(table$CpG, annotation$Name)
  if (anyNA(ix)) stop("Missing annotation for tested probes.")
  a <- annotation[ix, , drop = FALSE]
  gene_lists <- strsplit(ifelse(is.na(a$UCSC_RefGene_Name), "", a$UCSC_RefGene_Name),
                        ";", fixed = TRUE)
  counts <- lengths(gene_lists)
  idx <- rep(seq_len(nrow(table)), counts)
  out <- table[idx, , drop = FALSE]
  out$Gene_Symbol <- trimws(unlist(gene_lists, use.names = FALSE))
  # Features are paired with gene entries, not with CpG-island annotations.
  feature_lists <- lapply(seq_along(gene_lists), function(i) {
    f <- strsplit(ifelse(is.na(a$UCSC_RefGene_Group[i]), "", a$UCSC_RefGene_Group[i]),
                  ";", fixed = TRUE)[[1]]
    n <- length(gene_lists[[i]])
    if (length(f) == n) return(f)
    rep(paste(f, collapse = ";"), n) # retain ambiguity rather than misassign a feature
  })
  out$CpG_Feature <- unlist(feature_lists, use.names = FALSE)
  for (field in intersect(c("chr", "pos", "Relation_to_Island"), names(a))) {
    out[[field]] <- a[[field]][idx]
  }
  out <- unique(out[nzchar(out$Gene_Symbol), , drop = FALSE])
  rownames(out) <- NULL
  out
}
summarize_methylation_genes <- function(x) {
  x <- unique(x[, c("CpG", "Gene_Symbol", "deltaBeta")])
  if (!nrow(x)) return(data.frame(Gene_Symbol = character(), n_CpGs = integer(),
                                   mean_deltaBeta = numeric(), mixed_direction = logical()))
  do.call(rbind, lapply(split(x, x$Gene_Symbol), function(d) {
    data.frame(Gene_Symbol = d$Gene_Symbol[1], n_CpGs = nrow(d),
               mean_deltaBeta = mean(d$deltaBeta),
               mixed_direction = any(d$deltaBeta > 0) && any(d$deltaBeta < 0))
  }))
}
significant_methylation <- function(x, config) {
  !is.na(x$adj.P.Val) & x$adj.P.Val < config$fdr & abs(x$deltaBeta) > config$min_delta_beta
}
significant_rna <- function(x, config) {
  !is.na(x$FDR) & x$FDR < config$fdr & abs(x$logFC) > config$min_abs_rna_logfc
}
clean_symbols <- function(x) sort(unique(trimws(x[!is.na(x) & nzchar(trimws(x))])))
