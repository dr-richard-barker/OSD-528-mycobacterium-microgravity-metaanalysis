#!/usr/bin/env Rscript
# 03_wgcna_coexpression_network.R
# Exploratory weighted gene co-expression network analysis (WGCNA R package) on the
# primary-sample DESeq2 variance-stabilised matrix from 02_differential_expression.R.
#
# Eight samples is very few for co-expression network inference: every pairwise
# correlation rests on eight points. Results here are exploratory and are reported as
# such. The soft-thresholding power is chosen from the measured scale-free fit; if no
# power reaches R^2 >= 0.8 the power with the best fit is used and the shortfall is
# recorded in wgcna_summary.json rather than hidden.

suppressPackageStartupMessages({
  library(WGCNA)
  library(jsonlite)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads(4)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
proj <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
processed <- file.path(proj, "data", "processed")

vst <- read.delim(file.path(processed, "osd528_vst.tsv"), check.names = FALSE)
annot <- read.delim(file.path(processed, "osd528_counts_raw.tsv"), check.names = FALSE)[
  , c("gene_id", "gene_symbol", "protein")]
samples <- read.delim(file.path(proj, "config", "samples.tsv"))
samples <- samples[samples$primary == "yes", ]
expr <- t(as.matrix(vst[, samples$sample_id]))
colnames(expr) <- vst$gene_id

gsg <- goodSamplesGenes(expr, verbose = 0)
expr <- expr[gsg$goodSamples, gsg$goodGenes]

powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(expr, powerVector = powers, networkType = "signed", verbose = 0)
fit <- sft$fitIndices
fit$scale_free_R2 <- -sign(fit$slope) * fit$SFT.R.sq
write.table(fit, file.path(processed, "wgcna_soft_threshold.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
reached <- fit$Power[fit$scale_free_R2 >= 0.8]
power <- if (length(reached)) min(reached) else fit$Power[which.max(fit$scale_free_R2)]
power_r2 <- fit$scale_free_R2[fit$Power == power]

net <- blockwiseModules(expr, power = power, networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25, maxBlockSize = 10000,
                        numericLabels = FALSE, pamRespectsDendro = FALSE, verbose = 0)
colors <- net$colors
MEs <- orderMEs(net$MEs)

adj <- adjacency(expr, power = power, type = "signed")
ic <- intramodularConnectivity(adj, colors)
kme <- signedKME(expr, MEs)
kme_own <- mapply(function(g, m) if (m == "grey") NA else kme[g, paste0("kME", m)],
                  colnames(expr), colors)
assign <- data.frame(gene_id = colnames(expr), module = colors, kME = kme_own,
                     kTotal = ic$kTotal, kWithin = ic$kWithin)
assign <- merge(annot, assign, by = "gene_id")
assign <- assign[order(assign$module, -assign$kWithin), ]
write.table(assign, file.path(processed, "wgcna_module_assignments.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(sample_id = rownames(expr), MEs, check.names = FALSE),
            file.path(processed, "wgcna_module_eigengenes.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

traits <- data.frame(
  simulated_microgravity = as.numeric(samples$microgravity == "yes"),
  clinostat = as.numeric(samples$group == "Clinostat"),
  rpm = as.numeric(samples$group == "RPM"))
r <- cor(MEs, traits)
p <- corPvalueStudent(r, nSamples = nrow(expr))
mt <- data.frame(module = rep(rownames(r), ncol(r)), trait = rep(colnames(r), each = nrow(r)),
                 pearson_r = as.vector(r), pvalue = as.vector(p))
mt$padj <- p.adjust(mt$pvalue, method = "BH")
write.table(mt, file.path(processed, "wgcna_module_trait_correlations.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

sizes <- table(colors)
summary <- list(
  wgcna_version = as.character(packageVersion("WGCNA")), samples = rownames(expr),
  genes = ncol(expr), network_type = "signed", power = power,
  scale_free_R2_at_power = power_r2, scale_free_fit_reached = length(reached) > 0,
  max_scale_free_R2 = max(fit$scale_free_R2), min_module_size = 30, merge_cut_height = 0.25,
  modules = as.list(sizes[names(sizes) != "grey"]), unassigned_grey = unname(sizes["grey"]),
  module_trait_tests = nrow(mt), module_trait_padj_lt_0.05 = sum(mt$padj < 0.05))
write_json(summary, file.path(processed, "wgcna_summary.json"), auto_unbox = TRUE,
           pretty = TRUE, digits = NA)
cat("Power", power, "R2", round(power_r2, 3), "fit reached:", length(reached) > 0, "\n")
print(sizes)
print(mt[mt$padj < 0.05, ])
