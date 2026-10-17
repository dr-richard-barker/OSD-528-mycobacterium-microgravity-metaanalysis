#!/usr/bin/env Rscript
# 02_differential_expression.R
# DESeq2 differential expression on the kallisto estimated counts written by
# 00_quantify_real_osdr_rnaseq.py.
#
# Two analyses, from config/samples.tsv:
#   primary      samples with primary == "yes" (RFPRPM4 excluded; RPM n = 2)
#   sensitivity  all nine samples, with the inferred RFPRPM4 library
# Contrasts: Clinostat vs Ground, RPM vs Ground, Clinostat vs RPM.
# A gene is differentially expressed (DEG) when the Benjamini-Hochberg adjusted
# Wald p-value is < 0.05. Fold changes are reported both as the maximum-likelihood
# estimate and shrunk with apeglm.

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(jsonlite)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
proj <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
processed <- file.path(proj, "data", "processed")

counts <- read.delim(file.path(processed, "osd528_counts_raw.tsv"), check.names = FALSE)
samples <- read.delim(file.path(proj, "config", "samples.tsv"))
counts <- counts[counts$pseudo == "no", ]
stopifnot(!anyDuplicated(counts$gene_id))
mat <- round(as.matrix(counts[, samples$sample_id]))
rownames(mat) <- counts$gene_id
annot <- counts[, c("gene_id", "gene_symbol", "protein")]

contrasts <- list(
  clinostat_vs_ground = c("Clinostat", "Ground", "deg_3dclinostat_vs_static1g.tsv"),
  rpm_vs_ground       = c("RPM", "Ground", "deg_rpm2_vs_static1g.tsv"),
  clinostat_vs_rpm    = c("Clinostat", "RPM", "deg_3dclinostat_vs_rpm2.tsv")
)
alpha <- 0.05

run_analysis <- function(sample_rows, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  s <- samples[sample_rows, ]
  m <- mat[, s$sample_id]
  keep <- rowSums(m >= 10) >= 2
  m <- m[keep, ]
  s$group <- factor(s$group, levels = c("Ground", "Clinostat", "RPM"))
  rownames(s) <- s$sample_id

  dds <- DESeq(DESeqDataSetFromMatrix(m, s, ~ group), quiet = TRUE)
  summary <- list(samples = s$sample_id, genes_input = nrow(mat), genes_tested = nrow(m),
                  size_factors = as.list(sizeFactors(dds)), contrasts = list())
  tables <- list()
  for (name in names(contrasts)) {
    num <- contrasts[[name]][1]; den <- contrasts[[name]][2]
    d <- dds
    if (levels(d$group)[1] != den) {
      d$group <- relevel(d$group, ref = den)
      d <- nbinomWaldTest(d, quiet = TRUE)
    }
    coef <- paste0("group_", num, "_vs_", den)
    res <- results(d, name = coef, alpha = alpha)
    shr <- lfcShrink(d, coef = coef, res = res, type = "apeglm", quiet = TRUE)
    tab <- data.frame(gene_id = rownames(res), baseMean = res$baseMean,
                      log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                      log2FoldChange_shrunk = shr$log2FoldChange,
                      pvalue = res$pvalue, padj = res$padj)
    tab <- merge(annot, tab, by = "gene_id")
    tab$significant <- ifelse(!is.na(tab$padj) & tab$padj < alpha, "YES", "NO")
    tab <- tab[order(tab$padj, tab$pvalue, na.last = TRUE), ]
    write.table(tab, file.path(out_dir, contrasts[[name]][3]), sep = "\t", quote = FALSE,
                row.names = FALSE)
    sig <- tab$significant == "YES"
    summary$contrasts[[name]] <- list(
      numerator = num, denominator = den,
      genes_with_padj = sum(!is.na(tab$padj)),
      deg = sum(sig), up = sum(sig & tab$log2FoldChange > 0), down = sum(sig & tab$log2FoldChange < 0),
      deg_abs_shrunk_lfc_ge_1 = sum(sig & abs(tab$log2FoldChange_shrunk) >= 1))
    tables[[name]] <- tab
  }

  # Shared simulated-microgravity response: DEG in both simulators, same direction.
  a <- tables$clinostat_vs_ground; b <- tables$rpm_vs_ground
  ab <- merge(a[, c("gene_id", "log2FoldChange", "significant")],
              b[, c("gene_id", "log2FoldChange", "significant")], by = "gene_id",
              suffixes = c("_clinostat", "_rpm"))
  both <- ab$significant_clinostat == "YES" & ab$significant_rpm == "YES"
  same <- both & sign(ab$log2FoldChange_clinostat) == sign(ab$log2FoldChange_rpm)
  either <- ab$significant_clinostat == "YES" | ab$significant_rpm == "YES"
  summary$simulator_overlap <- list(
    deg_both = sum(both), deg_both_same_direction = sum(same), deg_either = sum(either),
    jaccard_same_direction = if (sum(either) > 0) sum(same) / sum(either) else NA,
    spearman_lfc_all_tested = cor(ab$log2FoldChange_clinostat, ab$log2FoldChange_rpm,
                                  method = "spearman", use = "complete.obs"))
  shared <- ab[same, ]
  shared <- merge(annot, shared, by = "gene_id")
  write.table(shared, file.path(out_dir, "deg_shared_simulated_microgravity.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)

  # Variance-stabilised expression (blind to design) and PCA on the 500 most variable genes.
  v <- assay(vst(dds, blind = TRUE))
  write.table(data.frame(gene_id = rownames(v), v, check.names = FALSE),
              file.path(out_dir, "osd528_vst.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  norm <- counts(dds, normalized = TRUE)
  write.table(data.frame(gene_id = rownames(norm), round(norm, 3), check.names = FALSE),
              file.path(out_dir, "osd528_counts_normalized.tsv"), sep = "\t", quote = FALSE,
              row.names = FALSE)
  top <- head(order(apply(v, 1, var), decreasing = TRUE), 500)
  pca <- prcomp(t(v[top, ]))
  ve <- pca$sdev^2 / sum(pca$sdev^2)
  write.table(data.frame(sample_id = rownames(pca$x), group = s[rownames(pca$x), "group"],
                         PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3]),
              file.path(out_dir, "pca_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  summary$pca <- list(genes = length(top), variance_explained = ve[1:min(5, length(ve))])
  list(summary = summary, tables = tables)
}

primary <- run_analysis(samples$primary == "yes", processed)
sens <- run_analysis(rep(TRUE, nrow(samples)), file.path(processed, "sensitivity_with_rfprpm4"))

# How much the primary DEG calls move when the inferred RFPRPM4 library is added.
concord <- list()
for (name in names(contrasts)) {
  p <- primary$tables[[name]]; q <- sens$tables[[name]]
  pq <- merge(p[, c("gene_id", "log2FoldChange", "significant")],
              q[, c("gene_id", "log2FoldChange", "significant")], by = "gene_id",
              suffixes = c("_primary", "_sensitivity"))
  concord[[name]] <- list(
    deg_primary = sum(pq$significant_primary == "YES"),
    deg_sensitivity = sum(pq$significant_sensitivity == "YES"),
    deg_both = sum(pq$significant_primary == "YES" & pq$significant_sensitivity == "YES"),
    pearson_lfc = cor(pq$log2FoldChange_primary, pq$log2FoldChange_sensitivity,
                      use = "complete.obs"))
}

out <- list(method = list(deseq2 = as.character(packageVersion("DESeq2")),
                          apeglm = as.character(packageVersion("apeglm")),
                          R = R.version.string, alpha = alpha,
                          prefilter = "counts >= 10 in at least 2 samples",
                          deg_definition = "BH-adjusted Wald p < 0.05"),
            primary = primary$summary, sensitivity_with_rfprpm4 = sens$summary,
            primary_vs_sensitivity = concord)
write_json(out, file.path(processed, "de_summary.json"), auto_unbox = TRUE, pretty = TRUE,
           digits = NA)
cat(toJSON(lapply(out$primary$contrasts, function(x) x[c("deg", "up", "down")]), auto_unbox = TRUE), "\n")
cat("Simulator overlap:", toJSON(out$primary$simulator_overlap, auto_unbox = TRUE), "\n")
cat("PCA variance explained:", round(100 * out$primary$pca$variance_explained, 1), "\n")
