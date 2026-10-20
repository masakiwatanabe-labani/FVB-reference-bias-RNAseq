#!/usr/bin/env Rscript
# edgeR (3.40.2) per tissue x arm: filterByExpr (defaults) -> TMM -> estimateDisp -> glmQLFit -> glmQLFTest.
# Usage: run_edger.R <counts.tsv> <samples.tsv> <formula> <coef_name> <out.tsv> [sample_filter_expr]
#   counts.tsv: gene_id + one column per sample (raw counts; may be fractional for -M --fraction, rounded here)
#   samples.tsv: sample_id, strain, diet, ...   formula e.g. "~ diet + strain"   coef_name e.g. "strainFVB"
#   sample_filter_expr (R expression on the sample sheet), e.g. 'diet == "water"' or 'strain == "FVB"'
suppressPackageStartupMessages(library(edgeR))
a <- commandArgs(trailingOnly = TRUE)
cnt <- read.delim(a[1], check.names = FALSE, row.names = 1); ss <- read.delim(a[2], stringsAsFactors = FALSE)
ss$strain <- factor(ifelse(ss$strain == "FVB/NJ", "FVB", ifelse(ss$strain == "C57BL/6J", "B6", ss$strain)), levels = c("B6", "FVB"))
if ("diet" %in% names(ss)) ss$diet <- factor(ss$diet, levels = unique(c(intersect(c("water", "control", "standard chow"), ss$diet), setdiff(ss$diet, c("water", "control", "standard chow")))))
if ("sex" %in% names(ss)) ss$sex <- factor(ss$sex)
if (length(a) >= 6 && nzchar(a[6])) ss <- ss[with(ss, eval(parse(text = a[6]))), ]
ss <- ss[ss$sample_id %in% colnames(cnt), ]; cnt <- round(as.matrix(cnt[, ss$sample_id, drop = FALSE]))
ss$strain <- droplevels(ss$strain); if ("diet" %in% names(ss)) ss$diet <- droplevels(ss$diet); if ("sex" %in% names(ss)) ss$sex <- droplevels(ss$sex)
design <- model.matrix(as.formula(a[3]), data = ss)
y <- DGEList(counts = cnt, samples = ss)
keep <- filterByExpr(y, design = design); y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); y <- estimateDisp(y, design)
fit <- glmQLFit(y, design); coef <- which(colnames(design) == a[4])
if (length(coef) != 1) stop(paste("coefficient not found:", a[4], "in", paste(colnames(design), collapse = ",")))
res <- topTags(glmQLFTest(fit, coef = coef), n = Inf, sort.by = "none")$table
res$gene_id <- rownames(res)
write.table(res[, c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")], a[5], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: n=%d samples (%s), genes kept %d/%d, coef %s, FDR<0.05: %d\n", basename(a[5]), ncol(cnt), paste(table(ss$strain), collapse = "/"), sum(keep), length(keep), a[4], sum(res$FDR < 0.05)))
