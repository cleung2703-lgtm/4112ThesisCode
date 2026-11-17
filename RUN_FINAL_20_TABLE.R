DATA_DIR <- "."
OUT_DIR  <- DATA_DIR
if (!file.exists(file.path(DATA_DIR, "vista_object_shared.rds"))) {
  if (file.exists("vista_object_shared.rds")) { DATA_DIR <- "."; OUT_DIR <- "."
  } else stop("vista_object_shared.rds not found in ", DATA_DIR,
              "\nRun RUN_VISTA_SELECTED_FIGURES.R first.")
}

CANDIDATES <- c("Cpd","Pitx2","Des","Mylk","Il1rl1","Spp1","Hmga1","Clu","Rasa1",
                "Hmga2","Rbp4","Ednrb","Serpinb6b","Fgf18","Smpd3","Tfap2b",
                "Hoxa7","Ablim1","Epas1","Hoxb6")

ROBUST_PATHWAYS <- c("MYOGENESIS", "IL2_STAT5_SIGNALING",
                     "EPITHELIAL_MESENCHYMAL_TRANSITION")
THESIS_A <- c("Smpd3","Pitx2","Epas1","Fgf18","Chd7","Rbp4","Snap25","Ednrb",
              "Hmga2","Oas1a","Oas1g","Dio2","Hmga1")

if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl", type = "binary")
suppressPackageStartupMessages({
  library(VISTA); library(SummarizedExperiment); library(writexl)
})

vista <- readRDS(file.path(DATA_DIR, "vista_object_shared.rds"))
comp  <- names(comparisons(vista))[1]
rd    <- as.data.frame(rowData(vista))
de    <- as.data.frame(comparisons(vista)[[comp]])
de$SYMBOL   <- unname(as.character(rd$SYMBOL)[match(de$gene_id, rownames(rd))])
de$GENENAME <- unname(as.character(rd$GENENAME)[match(de$gene_id, rownames(rd))])

miss <- setdiff(CANDIDATES, de$SYMBOL)
if (length(miss)) cat("NOT FOUND in the VISTA object:", paste(miss, collapse = ", "), fill = TRUE)

pick_row <- function(symbols, sym_vec, base_vec) {
  vapply(symbols, function(g) {
    k <- which(!is.na(sym_vec) & sym_vec == g)
    if (!length(k)) return(NA_integer_)
    if (length(k) > 1) {
      k <- k[order(base_vec[k], decreasing = TRUE)]
      cat(sprintf("  NOTE: %s maps to %d Ensembl IDs; using the highest-expressed"
                  , g, length(k)), fill = TRUE)
    }
    k[1]
  }, 1L)
}

i <- pick_row(CANDIDATES, de$SYMBOL, de$baseMean)
i <- i[!is.na(i)]
tbl <- data.frame(
  gene        = de$SYMBOL[i],
  gene_name   = de$GENENAME[i],
  ensembl     = de$gene_id[i],
  direction   = ifelse(de$log2fc[i] > 0, "Up in KO", "Down in KO"),
  log2FC      = round(de$log2fc[i], 3),
  fold_change = round(ifelse(de$log2fc[i] > 0, 2^de$log2fc[i], -2^abs(de$log2fc[i])), 2),
  padj        = signif(de$padj[i], 3),
  baseMean    = round(de$baseMean[i], 1),
  DEG_call    = de$regulation[i],
  stringsAsFactors = FALSE
)

nc  <- assay(vista, "norm_counts")[de$gene_id[i], , drop = FALSE]
grp <- as.character(colData(vista)$group)
ko  <- nc[, grp == "KO", drop = FALSE]; wt <- nc[, grp == "WT", drop = FALSE]
tbl$KO_mean <- round(rowMeans(ko), 1)
tbl$WT_mean <- round(rowMeans(wt), 1)
tbl$KO_cv   <- round(apply(ko, 1, sd) / pmax(rowMeans(ko), 1e-9), 2)
tbl$separation <- ifelse(apply(ko, 1, min) > apply(wt, 1, max) |
                         apply(ko, 1, max) < apply(wt, 1, min), "Clean", "Overlapping")
tbl <- cbind(tbl, round(as.data.frame(nc), 1))

gp <- file.path(DATA_DIR, "vista_tables", "go_bp_results.csv")
if (!file.exists(gp)) gp <- file.path(DATA_DIR, "03_go_bp_results.csv")
if (file.exists(gp)) {
  go <- read.csv(gp, stringsAsFactors = FALSE)
  go <- go[go$p.adjust < 0.05, ]
  hits <- strsplit(go$geneID, "/")
  tbl$n_GO_terms <- vapply(tbl$ensembl, function(e) sum(vapply(hits, function(h) e %in% h, TRUE)), 1L)
} else { tbl$n_GO_terms <- NA_integer_; cat("GO results not found\n") }

lp <- file.path(DATA_DIR, "vista_tables", "gsea_leading_edge_by_gene.csv")
if (!file.exists(lp)) lp <- file.path(DATA_DIR, "05_gsea_leading_edge_by_gene.csv")
if (file.exists(lp)) {
  le <- read.csv(lp, stringsAsFactors = FALSE)
  gcol <- if ("SYMBOL" %in% names(le)) "SYMBOL" else "gene"
  tbl$GSEA_pathways <- vapply(tbl$gene, function(g) {
    p <- unique(le$pathway[le[[gcol]] == g]); if (!length(p)) "" else paste(sort(p), collapse = "; ")
  }, "")
  tbl$n_GSEA_pathways <- vapply(strsplit(tbl$GSEA_pathways, "; "),
                                function(x) sum(nzchar(x)), 1L)
  tbl$in_robust_pathway <- grepl(paste(ROBUST_PATHWAYS, collapse = "|"), tbl$GSEA_pathways)
} else {
  tbl$GSEA_pathways <- ""; tbl$n_GSEA_pathways <- NA_integer_
  tbl$in_robust_pathway <- NA; cat("GSEA leading edge not found\n")
}

pp <- file.path(DATA_DIR, "Cpd_qPCR_primers.xlsx")
if (file.exists(pp)) {
  pr <- suppressMessages(readxl::read_excel(pp, col_names = FALSE))
  names(pr)[1:3] <- c("primer", "seq", "amplicon")
  pr <- pr[!is.na(pr$primer) & grepl("-5'", pr$primer) & !is.na(pr$amplicon), ]
  key <- toupper(gsub("^hm|-5'.*$", "", pr$primer))
  tbl$primer <- gsub("\\s+", " ", trimws(pr$amplicon))[match(toupper(tbl$gene), key)]
} else { tbl$primer <- NA_character_; cat("Primer file not found\n") }

tbl$thesis_A <- tbl$gene %in% THESIS_A

tbl$score <- (2 * (tbl$DEG_call %in% c("Up","Down"))) +
             (tbl$separation == "Clean") +
             (tbl$n_GSEA_pathways > 0) +
             (tbl$in_robust_pathway %in% TRUE) +
             (tbl$n_GO_terms >= 3) +
             (!is.na(tbl$primer)) +
             tbl$thesis_A
tbl <- tbl[order(-tbl$score, tbl$padj), ]

front <- c("gene","gene_name","direction","log2FC","fold_change","padj","DEG_call",
           "separation","KO_mean","WT_mean","KO_cv","n_GO_terms","n_GSEA_pathways",
           "in_robust_pathway","primer","thesis_A","score")
tbl <- tbl[, c(front, setdiff(names(tbl), front))]

cat("\n--- FINAL", nrow(tbl), "CANDIDATE GENES ---\n")
print(tbl[, c("gene","direction","log2FC","padj","separation","n_GO_terms",
              "n_GSEA_pathways","in_robust_pathway","primer","score")], row.names = FALSE)

cat(sprintf("\nAll DEGs: %s | all clean separation: %s | with primers: %d | Thesis A: %d\n",
            all(tbl$DEG_call %in% c("Up","Down")), all(tbl$separation == "Clean"),
            sum(!is.na(tbl$primer)), sum(tbl$thesis_A)))

write.csv(tbl, file.path(OUT_DIR, "FINAL_20_candidate_genes.csv"), row.names = FALSE)

meanings <- data.frame(column = front, meaning = c(
  "Gene symbol", "Full gene name", "Direction of change in KO vs WT",
  "log2 fold change from VISTA/DESeq2", "Linear fold change; negative = fold DOWN",
  "Benjamini-Hochberg adjusted p-value", "VISTA's call: Up/Down needs padj<0.05 AND |log2FC|>=1",
  "Clean = all 3 KO samples outside the full WT range (no outlier driving it)",
  "Mean normalised count, 3 KO samples", "Mean normalised count, 2 WT samples",
  "Coefficient of variation within KO; high = one sample differs",
  "Number of significant GO BP terms containing this gene",
  "Number of significant Hallmark GSEA leading edges containing it",
  "In MYOGENESIS, IL2_STAT5 or EMT (robust per the sensitivity analysis)",
  "Validated qPCR amplicon size, if one exists",
  "Appeared in the Thesis A shortlist or chord hub genes",
  "Sum of the evidence columns, max 8. Re-weight by editing the script."))
write_xlsx(list(candidates = tbl, column_meanings = meanings),
           file.path(OUT_DIR, "FINAL_20_candidate_genes.xlsx"))

cat("\nWrote: FINAL_20_candidate_genes.csv\n       FINAL_20_candidate_genes.xlsx\n")
