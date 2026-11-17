DATA_DIR <- "."
OUT_DIR  <- DATA_DIR

if (!file.exists(file.path(DATA_DIR, "MEF_CPD_KO_feature_count.txt"))) {
  if (file.exists("MEF_CPD_KO_feature_count.txt")) { DATA_DIR <- "."; OUT_DIR <- "."
  } else stop("Cannot find MEF_CPD_KO_feature_count.txt\n  Looked in: ", DATA_DIR,
              "\n  Fix: set DATA_DIR at the top, or run setwd to that folder.")
}

FIGDIR <- file.path(OUT_DIR, "vista_figures")
TABDIR <- file.path(OUT_DIR, "vista_tables")
dir.create(FIGDIR, showWarnings = FALSE); dir.create(TABDIR, showWarnings = FALSE)

PADJ_CUTOFF   <- 0.05
LOG2FC_CUTOFF <- 1
ONT           <- "BP"
TOP_N_ENRICH  <- 15
GSEA_SEED     <- 42
N_PERM        <- 100000
Y_CAP         <- 20
HEATMAP_N     <- 60
REBUILD       <- FALSE

COLORS <- c(Up = "#D73027", Down = "#2166AC", Other = "grey80")

CANDIDATES <- c("Cpd","Pitx2","Des","Mylk","Il1rl1","Spp1","Hmga1","Clu","Rasa1",
                "Hmga2","Rbp4","Ednrb","Serpinb6b","Fgf18","Smpd3","Tfap2b",
                "Hoxa7","Ablim1","Epas1","Hoxb6")

suppressPackageStartupMessages({
  library(VISTA); library(org.Mm.eg.db); library(SummarizedExperiment)
})

.status <- data.frame(figure = character(), ok = logical(), note = character())

fig <- function(name, fn, w = 10, h = 8, base_graphics = FALSE) {
  f <- file.path(FIGDIR, paste0(name, ".png"))
  r <- try({
    if (base_graphics) {
      grDevices::png(f, width = w, height = h, units = "in", res = 300)
      fn()
      grDevices::dev.off()
    } else {
      save_vista_plot(fn(), f, width = w, height = h, dpi = 300)
    }
  }, silent = TRUE)
  bad <- inherits(r, "try-error")
  if (bad && grDevices::dev.cur() > 1) grDevices::dev.off()
  msg <- if (bad) conditionMessage(attr(r, "condition")) else ""
  .status <<- rbind(.status, data.frame(figure = name, ok = !bad,
                                        note = substr(gsub("[ \t\r\n]+", " ", msg), 1, 80)))
  cat(sprintf("  %-42s %s\n", name,
              if (bad) paste("FAILED -", substr(msg, 1, 50)) else "ok"))
  invisible(!bad)
}

tab <- function(d, name) {
  utils::write.csv(d, file.path(TABDIR, paste0(name, ".csv")), row.names = FALSE)
  cat(sprintf("  %-42s %d rows\n", paste0(name, ".csv"), nrow(d)))
}

cat("\n== 1. VISTA object ==\n")
VISTA_RDS <- file.path(OUT_DIR, "vista_object_shared.rds")

if (!REBUILD && file.exists(VISTA_RDS)) {
  vista <- readRDS(VISTA_RDS); cat("Loaded cached object\n")
} else {
  cm <- as.matrix(read.delim(file.path(DATA_DIR, "MEF_CPD_KO_feature_count.txt"),
                             row.names = "gene_id", check.names = FALSE))
  sm <- data.frame(sample_names = colnames(cm),
                   group = ifelse(grepl("^CPD-KO", colnames(cm)), "KO", "WT"))
  stopifnot(sum(sm$group == "KO") == 3, sum(sm$group == "WT") == 2)
  mi <- match_vista_inputs(read_vista_counts(cm, format = "matrix"),
                           read_vista_metadata(sm))
  vista <- create_vista(
    counts = mi$counts, sample_info = mi$sample_info,
    column_geneid = mi$column_geneid, group_column = "group",
    group_numerator = "KO", group_denominator = "WT", method = "deseq2",
    log2fc_cutoff = LOG2FC_CUTOFF, pval_cutoff = PADJ_CUTOFF,
    p_value_type = "padj", min_counts = 10, min_replicates = 1)
  vista <- set_rowdata(vista, orgdb = org.Mm.eg.db,
                       columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                       keytype = "ENSEMBL")
  saveRDS(vista, VISTA_RDS); cat("Built and cached object\n")
}
comp <- names(comparisons(vista))[1]
cat("Comparison:", comp, "\n"); print(deg_summary(vista))

cat("\n== 2. QC ==\n")
fig("01_qc_pca",          function() get_pca_plot(vista), 8, 7)
fig("02_qc_mds",          function() get_mds_plot(vista), 8, 7)
fig("03_qc_corr_heatmap", function() get_corr_heatmap(vista), 8, 7)

cat("\n== 3. Differential expression ==\n")
fig("04_volcano", function()
  get_volcano_plot(vista, sample_comparison = comp, colors = COLORS,
                   display_id = "SYMBOL", label_size = 3, point_size = 1.2), 11, 9)

fig("05_volcano_zoomed", function()
  get_volcano_plot(vista, sample_comparison = comp, colors = COLORS,
                   display_id = "SYMBOL", label_size = 3, point_size = 1.2,
                   ylim = c(0, Y_CAP), xlim = c(-8, 8)), 11, 9)

fig("06_ma_plot", function()
  get_ma_plot(vista, sample_comparison = comp, colors = COLORS,
              display_id = "SYMBOL", label_n = 10, label_size = 3,
              point_size = 1.2, alpha = 0.6), 11, 8)

fig("07_deg_count_barplot", function() get_deg_count_barplot(vista), 7, 6)
fig("08_deg_count_pieplot", function() get_deg_count_pieplot(vista), 7, 6)
fig("09_deg_count_donut",   function() get_deg_count_donutplot(vista), 7, 6)

de <- as.data.frame(comparisons(vista)[[comp]])
rd <- as.data.frame(rowData(vista))
de$SYMBOL   <- unname(as.character(rd$SYMBOL)[match(de$gene_id, rownames(rd))])
de$GENENAME <- unname(as.character(rd$GENENAME)[match(de$gene_id, rownames(rd))])
keep <- c("gene_id","SYMBOL","GENENAME","baseMean","log2fc","lfcSE","stat",
          "pvalue","padj","regulation")
de <- de[order(de$padj), intersect(keep, colnames(de))]
tab(de, "de_results_all_genes")
tab(de[de$regulation %in% c("Up","Down"), ], "de_results_DEGs")

cat("\n== 4. Expression ==\n")
fig("10_expression_heatmap_top", function()
  get_expression_heatmap(vista, top_n = HEATMAP_N, display_id = "SYMBOL",
                         repair_genes = TRUE, summarise_replicates = FALSE,
                         show_row_names = TRUE, label_size = 6), 8, 13)

pick_row <- function(symbols, sym_vec, base_vec) {
  vapply(symbols, function(g) {
    k <- which(!is.na(sym_vec) & sym_vec == g)
    if (!length(k)) return(NA_integer_)
    if (length(k) > 1) k <- k[order(base_vec[k], decreasing = TRUE)]
    k[1]
  }, 1L)
}
CAND_ENS <- rownames(rd)[pick_row(CANDIDATES, as.character(rd$SYMBOL), rd$baseMean)]
if (anyNA(CAND_ENS)) {
  msg <- paste(CANDIDATES[is.na(CAND_ENS)], collapse = ", ")
  cat("  candidates with no Ensembl match:", msg, fill = TRUE)
}
CAND_ENS <- CAND_ENS[!is.na(CAND_ENS)]

fig("11_expression_heatmap_candidates", function()
  get_expression_heatmap(vista, genes = CAND_ENS, display_id = "SYMBOL",
                         repair_genes = TRUE, summarise_replicates = FALSE,
                         show_row_names = TRUE, label_size = 10), 8, 10)

fig("12_foldchange_boxplot_candidates", function()
  get_foldchange_boxplot(vista, genes = CANDIDATES, display_id = "SYMBOL"), 12, 8)

fig("13_foldchange_barplot_candidates", function()
  get_foldchange_barplot(vista, genes = CANDIDATES, display_id = "SYMBOL"), 11, 8)

cat("\n== 5. GO enrichment ==\n")
go_tabs <- list()
for (reg in c("Up", "Down", "Both")) {
  res <- try(get_go_enrichment(vista, sample_comparison = comp, regulation = reg,
                               ont = ONT, orgdb = org.Mm.eg.db,
                               species = "Mus musculus", from_type = "ENSEMBL"),
             silent = TRUE)
  if (inherits(res, "try-error")) { cat("  GO", reg, "FAILED\n"); next }
  obj <- if (!is.null(res[[reg]])) res[[reg]] else res[[1]]
  if (is.null(obj)) next
  d <- as.data.frame(obj@result)
  if (!nrow(d)) next
  d$regulation <- reg; go_tabs[[reg]] <- d
  nsig <- sum(d$p.adjust < PADJ_CUTOFF, na.rm = TRUE)
  cat(sprintf("  GO %-5s %d terms, %d significant\n", reg, nrow(d), nsig))
  if (nsig > 0) {
    local({
      o <- obj; r <- reg
      fig(sprintf("14_go_%s_%s_dotplot", tolower(ONT), tolower(r)), function()
        get_enrichment_plot(o, top_n = TOP_N_ENRICH,
          title = sprintf("GO %s enrichment - %s-regulated", ONT, r)), 11, 8)
    })
  }
}
if (length(go_tabs)) tab(do.call(rbind, go_tabs), paste0("go_", tolower(ONT), "_results"))

cat("\n== 6. MSigDB over-representation (ORA) ==\n")
ora_tabs <- list()
for (reg in c("Up", "Down", "Both")) {
  res <- try(get_msigdb_enrichment(vista, sample_comparison = comp,
                                   regulation = reg, from_type = "ENSEMBL",
                                   orgdb = org.Mm.eg.db, msigdb_category = "H",
                                   species = "Mus musculus"), silent = TRUE)
  if (inherits(res, "try-error")) { cat("  ORA", reg, "FAILED\n"); next }
  obj <- if (!is.null(res[[reg]])) res[[reg]] else res[[1]]
  if (is.null(obj)) next
  d <- as.data.frame(obj@result)
  if (!nrow(d)) next
  d$regulation <- reg; ora_tabs[[reg]] <- d
  nsig <- sum(d$p.adjust < PADJ_CUTOFF, na.rm = TRUE)
  cat(sprintf("  ORA %-5s %d sets, %d significant\n", reg, nrow(d), nsig))
  if (nsig > 0) {
    local({
      o <- obj; r <- reg
      fig(sprintf("15_msigdb_ora_%s_dotplot", tolower(r)), function()
        get_enrichment_plot(o, top_n = TOP_N_ENRICH,
          title = sprintf("Hallmark over-representation - %s", r)), 11, 8)
    })
  }
}
if (length(ora_tabs)) tab(do.call(rbind, ora_tabs), "msigdb_ora_results")

cat("\n== 7. GSEA (Hallmark) ==\n")
GSEA_RDS <- file.path(OUT_DIR, "vista_gsea_shared.rds")
if (!REBUILD && file.exists(GSEA_RDS)) {
  gsea <- readRDS(GSEA_RDS); cat("Loaded cached GSEA\n")
} else {
  set.seed(GSEA_SEED)
  gsea <- get_gsea(vista, sample_comparison = comp, set_type = "msigdb",
                   from_type = "ENSEMBL", msigdb_category = "H",
                   orgdb = org.Mm.eg.db, species = "Mus musculus",
                   pvalueCutoff = 1, nPermSimple = N_PERM, eps = 0)
  saveRDS(gsea, GSEA_RDS); cat("Built and cached GSEA (seed", GSEA_SEED, ")\n")
}
gt <- as.data.frame(gsea$enrich@result)
cat(sprintf("  %d Hallmark sets returned, %d significant at padj < %g\n",
            nrow(gt), sum(gt$p.adjust < PADJ_CUTOFF, na.rm = TRUE), PADJ_CUTOFF))
if (nrow(gt) < 50)
  cat("  NOTE:", 50 - nrow(gt), "Hallmark set(s) dropped by fgsea for unbalanced",
      "gene-level statistics.\n        Raising N_PERM can sometimes recover them.\n")
tab(gt[order(gt$p.adjust), ], "gsea_hallmark_results")

fig("16_gsea_dotplot", function()
  get_enrichment_plot(gsea$enrich, top_n = TOP_N_ENRICH,
                      title = "Hallmark GSEA - CPD KO vs WT"), 11, 8)

sig_paths <- gt$Description[!is.na(gt$p.adjust) & gt$p.adjust < PADJ_CUTOFF]
cat("  significant pathways:", length(sig_paths), "\n")

if (length(sig_paths)) {

  fig("17_pathway_heatmap", function()
    get_pathway_heatmap(vista, gsea, pathways = utils::head(sig_paths, 5),
                        max_genes = 80), 9, 11)

  pg <- as.data.frame(get_pathway_genes(gsea, pathways = sig_paths,
                                        return_type = "long"))
  gcol <- intersect(c("gene", "gene_id", "geneID"), colnames(pg))[1]
  sym  <- setNames(as.character(rd$SYMBOL), rownames(rd))
  pg$ensembl <- pg[[gcol]]
  pg$SYMBOL  <- unname(sym[pg$ensembl])
  pg$SYMBOL  <- ifelse(is.na(pg$SYMBOL) | pg$SYMBOL == "", pg$ensembl, pg$SYMBOL)
  m <- match(pg$ensembl, de$gene_id)
  pg$log2fc <- de$log2fc[m]; pg$padj <- de$padj[m]; pg$regulation <- de$regulation[m]
  pg$pathway <- sub("^HALLMARK_", "", pg$pathway)
  tab(pg[, c("pathway","SYMBOL","ensembl","log2fc","padj","regulation")],
      "gsea_leading_edge_by_gene")

  n_path <- aggregate(pathway ~ SYMBOL, data = pg,
                      FUN = function(v) length(unique(v)))
  names(n_path)[2] <- "n_pathways"
  paths <- aggregate(pathway ~ SYMBOL, data = pg,
                     FUN = function(v) paste(sort(unique(v)), collapse = "; "))
  hub <- merge(n_path, paths, by = "SYMBOL")
  m2 <- pick_row(hub$SYMBOL, de$SYMBOL, de$baseMean)
  hub$log2fc <- de$log2fc[m2]; hub$padj <- de$padj[m2]
  hub$regulation <- de$regulation[m2]
  tab(hub[order(-hub$n_pathways, hub$padj), ], "gsea_genes_by_pathway_count")

  up   <- utils::head(gt$Description[gt$NES > 0 & gt$p.adjust < PADJ_CUTOFF], 5)
  down <- utils::head(gt$Description[gt$NES < 0 & gt$p.adjust < PADJ_CUTOFF], 5)
  chord_paths <- c(up, down)
  cat("  chord pathways:", length(chord_paths), "\n")
  if (length(chord_paths) >= 2)
    fig("18_gsea_chord", function() {
      par(xpd = NA)
      get_enrichment_chord(gsea$enrich, vista = vista, sample_comparison = comp,
                           pathways = chord_paths, max_genes = 70,
                           gene_order_by = "foldchange", color_by = "regulation",
                           up_color = COLORS[["Up"]], down_color = COLORS[["Down"]],
                           display_id = "SYMBOL", label_cex = 0.45, gap_degree = 2,
                           title = "CPD KO vs WT: Hallmark pathway-gene network")
    }, 13, 13, base_graphics = TRUE)
}

cat("\n== Summary ==\n")
print(.status, row.names = FALSE)
cat(sprintf("\n%d of %d figures written to %s\n", sum(.status$ok), nrow(.status),
            normalizePath(FIGDIR, mustWork = FALSE)))
cat("Tables written to", normalizePath(TABDIR, mustWork = FALSE), "\n")
if (any(!.status$ok)) {
  cat("\nFailures to look at:\n")
  print(.status[!.status$ok, c("figure", "note")], row.names = FALSE)
}
cat("\nDone.\n")
