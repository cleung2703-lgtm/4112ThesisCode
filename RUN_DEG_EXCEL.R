DATA_DIR <- "."
OUT_DIR  <- DATA_DIR
if (!file.exists(file.path(DATA_DIR, "vista_object_shared.rds"))) {
  if (file.exists("vista_object_shared.rds")) { DATA_DIR <- "."; OUT_DIR <- "."
  } else stop("vista_object_shared.rds not found in ", DATA_DIR,
              "\nRun RUN_VISTA_SELECTED_FIGURES.R first to build it.")
}

if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl", type = "binary")
suppressPackageStartupMessages({
  library(VISTA); library(SummarizedExperiment); library(writexl)
})

vista <- readRDS(file.path(DATA_DIR, "vista_object_shared.rds"))
comp  <- names(comparisons(vista))[1]
rd    <- as.data.frame(rowData(vista))
de    <- as.data.frame(comparisons(vista)[[comp]])
cat("VISTA object:", nrow(de), "genes |", comp, "\n")

raw  <- assay(vista, "counts")
norm <- assay(vista, "norm_counts")
grp  <- as.character(colData(vista)$group)
ko_c <- colnames(raw)[grp == "KO"]; wt_c <- colnames(raw)[grp == "WT"]
cat("samples:", length(ko_c), "KO (", paste(ko_c, collapse = ", "), ") |",
    length(wt_c), "WT (", paste(wt_c, collapse = ", "), ")\n")

i <- match(de$gene_id, rownames(raw))
rawm <- raw[i, , drop = FALSE]; normm <- norm[i, , drop = FALSE]
ko_n <- normm[, grp == "KO", drop = FALSE]; wt_n <- normm[, grp == "WT", drop = FALSE]

out <- data.frame(
  gene_id     = de$gene_id,
  gene_name   = unname(as.character(rd$SYMBOL)[match(de$gene_id, rownames(rd))]),
  description = unname(as.character(rd$GENENAME)[match(de$gene_id, rownames(rd))]),
  stringsAsFactors = FALSE
)
out <- cbind(out, as.data.frame(rawm))
names(out)[4:(3 + ncol(rawm))] <- paste0("raw_", colnames(rawm))
out <- cbind(out, round(as.data.frame(normm), 2))
names(out)[(4 + ncol(rawm)):(3 + 2 * ncol(rawm))] <- paste0("norm_", colnames(normm))

out$KO_avg <- round(rowMeans(ko_n), 2)
out$WT_avg <- round(rowMeans(wt_n), 2)
out$baseMean <- round(de$baseMean, 2)
out$log2FC <- round(de$log2fc, 4)
out$fold_change <- round(ifelse(de$log2fc > 0, 2^de$log2fc, -2^abs(de$log2fc)), 3)
out$lfcSE <- round(de$lfcSE, 4)
out$stat  <- round(de$stat, 4)
out$p_value <- de$pvalue
out$minuslog_p_value <- round(-log10(pmax(de$pvalue, .Machine$double.xmin)), 4)
out$FDR <- de$padj
out$minuslog_FDR <- round(-log10(pmax(de$padj, .Machine$double.xmin)), 4)
out$regulation <- de$regulation

out$separation <- ifelse(apply(ko_n, 1, min) > apply(wt_n, 1, max) |
                         apply(ko_n, 1, max) < apply(wt_n, 1, min),
                         "Clean", "Overlapping")
out$KO_cv <- round(apply(ko_n, 1, sd) / pmax(rowMeans(ko_n), 1e-9), 3)

out <- out[order(out$FDR, -abs(out$log2FC)), ]

degs <- out[out$regulation %in% c("Up", "Down"), ]
ups  <- out[out$regulation == "Up", ]
dns  <- out[out$regulation == "Down", ]

cat(sprintf("\nAll genes tested : %d\n", nrow(out)))
cat(sprintf("DEGs             : %d  (%d up, %d down)\n", nrow(degs), nrow(ups), nrow(dns)))
cat(sprintf("DEGs with clean sample separation: %d of %d\n",
            sum(degs$separation == "Clean"), nrow(degs)))

meanings <- data.frame(
  column = c("gene_id","gene_name","description",
             paste0("raw_", colnames(rawm)[1]), paste0("norm_", colnames(normm)[1]),
             "KO_avg","WT_avg","baseMean","log2FC","fold_change","lfcSE","stat",
             "p_value","minuslog_p_value","FDR","minuslog_FDR","regulation",
             "separation","KO_cv"),
  meaning = c(
    "Ensembl gene ID (the identifier the analysis is keyed on)",
    "Gene symbol from org.Mm.eg.db",
    "Full gene name from org.Mm.eg.db",
    "Raw read count for that sample (one column per sample)",
    "DESeq2 size-factor normalised count (one column per sample)",
    "Mean normalised count across the 3 KO samples",
    "Mean normalised count across the 2 WT samples",
    "DESeq2 baseMean: mean normalised count across all 5 samples",
    "log2 fold change, KO relative to WT",
    "Linear fold change; negative means fold DOWN in KO",
    "Standard error of the log2 fold change",
    "Wald test statistic",
    "Raw p-value from the Wald test",
    "-log10 of the raw p-value",
    "Benjamini-Hochberg adjusted p-value (equivalent to the FDR column in the lab file)",
    "-log10 of the adjusted p-value",
    "VISTA's call: Up or Down needs FDR < 0.05 AND |log2FC| >= 1; everything else is Other",
    "Clean = all 3 KO values fall outside the full WT range, so no single sample is driving it",
    "Coefficient of variation within the KO group; high values mean one KO sample differs"),
  stringsAsFactors = FALSE)

method <- data.frame(item = c(
  "Input", "Samples", "Genes in input", "Prefilter", "Genes tested",
  "DE method", "Package", "Significance rule", "Comparison", "Normalisation"),
  value = c(
    "MEF_CPD_KO_feature_count.txt (raw counts)",
    paste(length(ko_c), "KO vs", length(wt_c), "WT mouse embryonic fibroblasts"),
    "48,709", "min_counts = 10, min_replicates = 1",
    format(nrow(out), big.mark = ","),
    "DESeq2 (Wald test), run through VISTA",
    paste("VISTA", as.character(utils::packageVersion("VISTA"))),
    "FDR < 0.05 and |log2FC| >= 1", comp,
    "DESeq2 median-of-ratios size factors"),
  stringsAsFactors = FALSE)

write_xlsx(list(DEGs = degs, Upregulated = ups, Downregulated = dns,
                All_genes = out, Column_meanings = meanings, Method = method),
           file.path(OUT_DIR, "VISTA_DEG_results.xlsx"))
write.csv(degs, file.path(OUT_DIR, "VISTA_DEG_results.csv"), row.names = FALSE)

cat("\nWrote: VISTA_DEG_results.xlsx  (6 sheets)\n       VISTA_DEG_results.csv\n")
cat("\nTop 10 DEGs by FDR:\n")
print(degs[1:10, c("gene_name","log2FC","fold_change","p_value","FDR","regulation","separation")],
      row.names = FALSE)
