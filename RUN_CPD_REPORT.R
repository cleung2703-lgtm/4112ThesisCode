DATA_DIR <- "."
OUT_DIR  <- file.path(DATA_DIR, "CPD_report")
ASSETS   <- file.path(OUT_DIR, "CPD_report_assets")
PRIMERS  <- file.path(DATA_DIR, "Cpd_KO_qRT-PCR_primers_design.xlsx")

PADJ_CUTOFF  <- 0.05
LFC_CUTOFF   <- 1
TOP_N_ENRICH <- 15
GSEA_SEED    <- 42
N_PERM       <- 10000
MAX_GENES_SHOWN <- 10

LABEL_GENES <- "Cpd"

COLORS <- c(Up = "#D73027", Down = "#2166AC", Other = "grey80")

for (p in c("quarto","DT","plotly","htmlwidgets","readxl"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, type = "binary")
suppressPackageStartupMessages({
  library(VISTA); library(SummarizedExperiment); library(org.Mm.eg.db)
  library(ggplot2); library(plotly)
})

old_wd <- setwd(DATA_DIR)
dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)
if (is.null(quarto::quarto_path())) {
  q <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"
  if (file.exists(q)) Sys.setenv(QUARTO_PATH = q) else
    stop("Quarto CLI not found. Install RStudio or Quarto (https://quarto.org).")
}
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
`%||%` <- function(a, b) if (is.null(a)) b else a
cat("VISTA", as.character(packageVersion("VISTA")), "| Quarto",
    as.character(quarto::quarto_version()), "\n\n")

pf <- readxl::read_excel(PRIMERS, sheet = 1, col_names = FALSE, .name_repair = "minimal")
pv <- as.character(pf[[1]]); pv <- pv[!is.na(pv)]
PANEL <- unique(sub("-[35]'Q$", "", pv[grepl("-[35]'Q$", pv)]))
cat(sprintf("Gene panel from %s: %d genes\n", basename(PRIMERS), length(PANEL)))

if (!file.exists("vista_object_shared.rds"))
  stop("vista_object_shared.rds not found in ", DATA_DIR)
vista <- readRDS("vista_object_shared.rds")
rd    <- as.data.frame(rowData(vista))
SYM   <- function(ids) { s <- as.character(rd$SYMBOL)[match(ids, rownames(rd))]
                         ifelse(is.na(s) | s == "", ids, s) }
add_name <- function(df) { df$gene_name <- SYM(df$gene_id)
  df[, c("gene_id","gene_name", setdiff(names(df), c("gene_id","gene_name")))] }
md <- S4Vectors::metadata(vista)
md$de_results <- S4Vectors::SimpleList(lapply(md$de_results, add_name))
md$de_results_by_method <- lapply(md$de_results_by_method,
  function(b) S4Vectors::SimpleList(lapply(b, add_name)))
S4Vectors::metadata(vista) <- md
validate_vista(vista)

comp <- names(comparisons(vista))[1]
de   <- as.data.frame(comparisons(vista)[[comp]])
nc   <- as.data.frame(norm_counts(vista))
si   <- as.data.frame(sample_info(vista))
grp  <- as.character(si$group)
KO   <- colnames(nc)[grp == "KO"]; WT <- colnames(nc)[grp == "WT"]
n_up <- sum(de$regulation == "Up"); n_dn <- sum(de$regulation == "Down")
n_deg <- n_up + n_dn; n_tested <- nrow(de)
fmt <- function(x) format(x, big.mark = ",")
cat(sprintf("%s | %s genes tested | %d DEGs (%d up, %d down) | %d KO vs %d WT\n\n",
            comp, fmt(n_tested), n_deg, n_up, n_dn, length(KO), length(WT)))

pick_row <- function(symbols, tbl) {
  degs <- which(tbl$regulation %in% c("Up","Down"))
  vapply(symbols, function(g) {
    k <- degs[tbl$gene_name[degs] == g]
    if (!length(k)) k <- which(tbl$gene_name == g)
    if (!length(k)) return(NA_integer_)
    k[order(tbl$baseMean[k], decreasing = TRUE)][1]
  }, 1L)
}

gene_hover <- function(gene_id, lfc, padj, reg) {
  i <- match(gene_id, rownames(nc))
  j <- function(cols) apply(round(nc[i, cols, drop = FALSE], 1), 1,
         function(r) paste(formatC(r, format = "f", digits = 1), collapse = ", "))
  sprintf("<b>%s</b><br>log2FC: %s<br>FDR: %s<br>%s<br>KO: %s<br>WT: %s",
    SYM(gene_id), formatC(lfc, format = "f", digits = 2),
    ifelse(is.na(padj), "not tested", format(signif(padj, 3), scientific = TRUE)),
    reg, j(KO), j(WT))
}

plain_axes <- function(p, x, y) p + labs(x = x, y = y,
                                         title = NULL, subtitle = NULL, caption = NULL)

add_hover <- function(p, hover) {
  p$data$hover <- hover
  k <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), TRUE))
  if (!length(k)) stop("no point layer found")
  for (i in k) {
    m <- p$layers[[i]]$mapping
    p$layers[[i]]$mapping <- if (is.null(m)) ggplot2::aes(text = hover) else
                             modifyList(m, ggplot2::aes(text = hover))
  }
  p
}

repel_to_text <- function(p, size = 3.2) {
  k <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), TRUE))
  if (!length(k)) return(p)
  p$layers <- p$layers[-k]
  yv <- rlang::eval_tidy(p$mapping$y, p$data)
  p + ggplot2::geom_text(size = size, show.legend = FALSE,
                         nudge_y = 0.075 * diff(range(yv, na.rm = TRUE)))
}

tidy_legend <- function(w) {
  for (i in seq_along(w$x$data)) {
    n <- w$x$data[[i]]$name
    if (!is.null(n) && nzchar(n)) w$x$data[[i]]$name <- gsub("^\\(|,\\d+\\)$|\\)$", "", n)
  }
  w
}

interactive_plot <- function(p, tooltip = "text", height = 560) {
  w <- tryCatch(suppressWarnings(ggplotly(p, tooltip = tooltip, height = height)),
                error = function(e) { cat("    not convertible:", conditionMessage(e), "\n"); NULL })
  if (is.null(w)) return(NULL)
  w <- tidy_legend(w)
  w <- plotly::config(w, displaylogo = FALSE,
         modeBarButtonsToRemove = c("select2d","lasso2d","autoScale2d"))
  plotly::layout(w, hoverlabel = list(align = "left"),
                 margin = list(l = 60, r = 30, t = 30, b = 60))
}

save_png <- function(expr, file, width = 9, height = 7, res = 200, base = FALSE) {
  path <- file.path(ASSETS, file)
  png(path, width = width * res, height = height * res, res = res)
  on.exit(dev.off(), add = TRUE)
  if (base) force(expr) else print(expr)
  invisible(file.path("CPD_report_assets", file))
}
save_csv <- function(df, file) {
  write.csv(df, file.path(ASSETS, file), row.names = FALSE)
  invisible(file.path("CPD_report_assets", file))
}

shorten <- function(x, n = MAX_GENES_SHOWN) vapply(strsplit(x %||% "", ", "), function(g) {
  g <- g[nzchar(g)]
  if (!length(g)) "" else if (length(g) <= n) paste(g, collapse = ", ")
  else sprintf("%s ... (+%d more)", paste(g[seq_len(n)], collapse = ", "), length(g) - n)
}, "")

B <- list(plots = list(), widgets = list(), tables = list(), meta = list())

cat("== quality control ==\n")
p_pca <- get_pca_plot(vista, label = TRUE); p_mds <- get_mds_plot(vista, label = TRUE)
p_corr <- get_corr_heatmap(vista)
B$plots$pca  <- save_png(p_pca,  "qc_pca.png", 9, 6)
B$plots$mds  <- save_png(p_mds,  "qc_mds.png", 9, 6)
B$plots$corr <- save_png(p_corr, "qc_corr_heatmap.png", 8, 7)
B$widgets$pca  <- interactive_plot(repel_to_text(p_pca), tooltip = c("x","y","colour","label"))
B$widgets$mds  <- interactive_plot(repel_to_text(p_mds), tooltip = c("x","y","colour","label"))
B$widgets$corr <- interactive_plot(p_corr, tooltip = c("x","y","fill"))

cat("== differential expression ==\n")
B$plots$pie <- save_png(get_deg_count_pieplot(vista, colors = COLORS[c("Up","Down")]),
                        "de_pie.png", 7, 6)

p_vol <- suppressWarnings(get_volcano_plot(vista, sample_comparison = comp,
           colors = COLORS, display_id = "SYMBOL", label_genes = LABEL_GENES,
           label_size = 3.5, point_size = 1.2))
p_vol <- p_vol + labs(colour = "Regulation")
B$plots$volcano <- save_png(p_vol, "de_volcano.png", 10, 8)
pv <- plain_axes(p_vol, "log2 fold change", "-log10 adjusted p-value")
pv <- add_hover(pv, gene_hover(pv$data$gene_id, pv$data$log2fc, pv$data$padj, pv$data$regulation))
B$widgets$volcano <- interactive_plot(pv, height = 620)

p_ma <- get_ma_plot(vista, sample_comparison = comp, colors = COLORS,
         display_id = "SYMBOL", label_n = 0, point_size = 1.2, alpha = 0.6)
p_ma <- p_ma + labs(colour = "Regulation")
ma_lab <- p_ma$data[p_ma$data$display_gene %in% LABEL_GENES, ]
if (nrow(ma_lab))
  p_ma <- p_ma + geom_text(data = ma_lab, inherit.aes = FALSE, size = 3.5, vjust = -1.2,
            aes(x = log10_baseMean, y = log2fc, label = display_gene))
B$plots$ma <- save_png(p_ma, "de_ma.png", 10, 7)
pm <- plain_axes(p_ma, "log10(mean expression + 1)", "log2 fold change")
pm <- add_hover(pm, gene_hover(pm$data$gene_id, pm$data$log2fc, pm$data$padj, pm$data$regulation))
B$widgets$ma <- interactive_plot(pm, height = 620)

ord <- function(d) d[order(d$padj, -abs(d$log2fc)), ]
de_s <- ord(de)
round_de <- function(d) { d$baseMean <- round(d$baseMean,1); d$log2fc <- round(d$log2fc,3)
  d$lfcSE <- round(d$lfcSE,3); d$stat <- round(d$stat,3)
  d$pvalue <- signif(d$pvalue,3); d$padj <- signif(d$padj,3); d }
B$tables$de_top  <- round_de(head(de_s, 200))
B$tables$de_up   <- round_de(ord(de[de$regulation == "Up", ]))
B$tables$de_down <- round_de(ord(de[de$regulation == "Down", ]))
save_csv(round_de(de_s), "de_results_all_genes.csv")
save_csv(B$tables$de_up, "de_upregulated.csv")
save_csv(B$tables$de_down, "de_downregulated.csv")

cat("== panel heatmap ==\n")
idx  <- pick_row(PANEL, de)
miss <- PANEL[is.na(idx)]
if (length(miss)) cat("  WARNING not found:", paste(miss, collapse = ", "), "\n")
pt <- de[idx[!is.na(idx)], ]
B$plots$heatmap <- save_png(
  get_expression_heatmap(vista, genes = pt$gene_id, display_id = "SYMBOL",
    value_transform = "zscore", show_row_names = TRUE, label_size = 9,
    annotate_columns = FALSE, summarise_replicates = FALSE,
    cluster_columns = TRUE, return_type = "plot"),
  "expr_heatmap_panel.png", 7, 10)
cat(sprintf("  %d of %d panel genes\n", nrow(pt), length(PANEL)))

B$tables$panel <- data.frame(
  Gene = pt$gene_name, `Ensembl ID` = pt$gene_id, baseMean = round(pt$baseMean, 1),
  log2FC = round(pt$log2fc, 3), FDR = signif(pt$padj, 3), Regulation = pt$regulation,
  check.names = FALSE)[order(pt$padj), ]
save_csv(B$tables$panel, "gene_panel.csv")

cat("== functional enrichment ==\n")
cached <- function(file, fn) if (file.exists(file)) readRDS(file) else { o <- fn(); saveRDS(o, file); o }
pick <- function(res) if (!is.null(res$Both)) res$Both else res[[1]]

go   <- cached("vista_go_both_shared.rds", function()
         get_go_enrichment(vista, sample_comparison = comp, regulation = "Both",
           ont = "BP", orgdb = org.Mm.eg.db, species = "Mus musculus", from_type = "ENSEMBL"))
ora  <- cached("vista_ora_both_shared.rds", function()
         get_msigdb_enrichment(vista, sample_comparison = comp, regulation = "Both",
           from_type = "ENSEMBL", orgdb = org.Mm.eg.db, msigdb_category = "H",
           species = "Mus musculus"))
gsea <- cached("vista_gsea_shared.rds", function() { set.seed(GSEA_SEED)
         get_gsea(vista, sample_comparison = comp, set_type = "msigdb",
           from_type = "ENSEMBL", msigdb_category = "H", orgdb = org.Mm.eg.db,
           species = "Mus musculus", pvalueCutoff = 1, nPermSimple = N_PERM, eps = 0) })
go_e <- pick(go); ora_e <- pick(ora); gsea_e <- gsea$enrich

enr_table <- function(enr) {
  d <- as.data.frame(enr@result, stringsAsFactors = FALSE)
  if (!nrow(d)) return(d)
  gcol <- intersect(c("geneID","core_enrichment"), names(d))[1]
  if (!is.na(gcol)) d$genes <- vapply(strsplit(d[[gcol]], "/"), function(g)
      paste(unique(SYM(g)), collapse = ", "), "")
  keep <- intersect(c("ID","Description","GeneRatio","setSize","NES","Count",
                      "pvalue","p.adjust","genes"), names(d))
  d <- d[, keep, drop = FALSE]
  for (k in intersect(c("pvalue","p.adjust"), names(d))) d[[k]] <- signif(d[[k]], 3)
  if ("NES" %in% names(d)) d$NES <- round(d$NES, 3)
  d
}
B$tables$go <- enr_table(go_e); B$tables$ora <- enr_table(ora_e); B$tables$gsea <- enr_table(gsea_e)
for (k in c("go","ora","gsea")) save_csv(B$tables[[k]], paste0("enrichment_", k, ".csv"))

B$tables$go_view   <- transform(B$tables$go,   genes = shorten(B$tables$go$genes))
B$tables$ora_view  <- transform(B$tables$ora,  genes = shorten(B$tables$ora$genes))
B$tables$gsea_view <- transform(B$tables$gsea, genes = shorten(B$tables$gsea$genes))

nsig <- function(d) sum(d$p.adjust < PADJ_CUTOFF, na.rm = TRUE)
cat(sprintf("  GO %d/%d sig | Hallmark ORA %d/%d | GSEA %d/%d\n",
  nsig(B$tables$go), nrow(B$tables$go), nsig(B$tables$ora), nrow(B$tables$ora),
  nsig(B$tables$gsea), nrow(B$tables$gsea)))

B$plots$go   <- save_png(get_enrichment_plot(go_e,   top_n = TOP_N_ENRICH,
                  title = paste("GO biological process -", comp)), "enrichment_go.png", 10, 7)
B$plots$ora  <- save_png(get_enrichment_plot(ora_e,  top_n = TOP_N_ENRICH,
                  title = paste("Hallmark over-representation -", comp)), "enrichment_ora.png", 10, 7)
B$plots$gsea <- save_png(get_enrichment_plot(gsea_e, top_n = TOP_N_ENRICH,
                  title = "Hallmark GSEA - CPD KO vs WT"), "enrichment_gsea.png", 10, 7)

sig_sets <- rbind(
  cbind(B$tables$ora[B$tables$ora$p.adjust < PADJ_CUTOFF, c("ID","p.adjust","genes")],
        Test = "Over-representation"),
  cbind(B$tables$gsea[B$tables$gsea$p.adjust < PADJ_CUTOFF, c("ID","p.adjust","genes")],
        Test = "GSEA"))
if (nrow(sig_sets)) {
  gl <- strsplit(sig_sets$genes, ", ")
  pg <- data.frame(
    Pathway = sub("^HALLMARK_", "", sig_sets$ID), Test = sig_sets$Test,
    FDR = sig_sets$p.adjust, Genes = lengths(gl),
    `Panel genes` = vapply(gl, function(g) paste(intersect(PANEL, g), collapse = ", "), ""),
    `Other genes` = shorten(vapply(gl, function(g) paste(setdiff(g, PANEL), collapse = ", "), "")),
    check.names = FALSE, stringsAsFactors = FALSE)
  pg <- pg[order(pg$Test, pg$FDR), ]
  B$tables$pathway_genes <- pg
  save_csv(data.frame(pg[, c("Pathway","Test","FDR","Genes","Panel genes")],
                      `All genes` = sig_sets$genes[order(sig_sets$Test, sig_sets$p.adjust)],
                      check.names = FALSE), "pathway_gene_membership.csv")
  cat(sprintf("  pathway membership: %d significant sets (%d panel genes involved)\n",
              nrow(pg), length(unique(unlist(strsplit(pg$`Panel genes`, ", "))))))
}

cat("== chord diagrams ==\n")
B$plots$chord_go <- save_png({ par(xpd = NA)
  get_enrichment_chord(go_e, vista = vista, sample_comparison = comp,
    top_n = 12, min_pathways = 3, max_genes = 70, gene_order_by = "foldchange",
    color_by = "regulation", up_color = COLORS[["Up"]], down_color = COLORS[["Down"]],
    display_id = "SYMBOL", label_cex = 0.5, gap_degree = 2,
    title = "Hub genes bridging enriched GO pathways (CPD KO vs WT)") },
  "chord_go_hub_genes.png", 12, 12, base = TRUE)

top12 <- head(B$tables$go[order(B$tables$go$p.adjust), ], 12)
hc <- sort(table(unlist(lapply(strsplit(top12$genes, ", "), unique))), decreasing = TRUE)
hub <- data.frame(Gene = names(hc), Pathways = as.integer(hc),
                  check.names = FALSE, stringsAsFactors = FALSE)
hub <- hub[hub$Pathways >= 3, ]
mm <- de[pick_row(hub$Gene, de), ]
hub$log2FC <- round(mm$log2fc, 3); hub$FDR <- signif(mm$padj, 3)
hub$Regulation <- mm$regulation
hub$`In panel` <- ifelse(hub$Gene %in% PANEL, "yes", "")
B$tables$hub <- hub; save_csv(hub, "chord_go_hub_genes.csv")
cat(sprintf("  GO chord: top 12 pathways, %d hub genes in >= 3\n", nrow(hub)))

gt <- B$tables$gsea; sg <- gt[!is.na(gt$p.adjust) & gt$p.adjust < PADJ_CUTOFF, ]
chord_paths <- c(head(sg$Description[sg$NES > 0], 5), head(sg$Description[sg$NES < 0], 5))
if (length(chord_paths) >= 2) {
  B$plots$chord_gsea <- save_png({ par(xpd = NA)
    get_enrichment_chord(gsea_e, vista = vista, sample_comparison = comp,
      pathways = chord_paths, max_genes = 70, gene_order_by = "foldchange",
      color_by = "regulation", up_color = COLORS[["Up"]], down_color = COLORS[["Down"]],
      display_id = "SYMBOL", label_cex = 0.45, gap_degree = 2,
      title = "CPD KO vs WT: Hallmark pathway-gene network") },
    "chord_hallmark_gsea.png", 12, 12, base = TRUE)
  cat(sprintf("  Hallmark chord: %d pathways\n", length(chord_paths)))
}

cat("== normalised counts ==\n")
i <- match(de$gene_id, rownames(nc))
cts <- data.frame(gene_id = de$gene_id, gene_name = de$gene_name,
  round(nc[i, c(KO, WT)], 2),
  KO_mean = round(rowMeans(nc[i, KO, drop = FALSE]), 2),
  WT_mean = round(rowMeans(nc[i, WT, drop = FALSE]), 2),
  log2FC = round(de$log2fc, 3), FDR = signif(de$padj, 3),
  regulation = de$regulation, check.names = FALSE, stringsAsFactors = FALSE)
B$tables$counts <- cts[order(cts$FDR, -abs(cts$log2FC)), ]
save_csv(B$tables$counts, "normalised_counts_per_replicate.csv")
cat(sprintf("  %s genes x %d replicates\n", fmt(nrow(cts)), length(KO) + length(WT)))

B$meta$summary <- data.frame(
  Item = c("Comparison","Samples","Genes tested","Differentially expressed",
           "Significance rule","Method","Gene panel"),
  Value = c(comp,
    sprintf("%d KO (%s) vs %d WT (%s)", length(KO), paste(KO, collapse = ", "),
            length(WT), paste(WT, collapse = ", ")),
    fmt(n_tested), sprintf("%d (%d up, %d down)", n_deg, n_up, n_dn),
    sprintf("FDR < %s and |log2FC| >= %s", PADJ_CUTOFF, LFC_CUTOFF),
    sprintf("DESeq2 Wald test via VISTA %s", as.character(packageVersion("VISTA"))),
    sprintf("%d genes, from %s", length(PANEL), basename(PRIMERS))),
  stringsAsFactors = FALSE)
B$meta$counts <- list(n_tested = n_tested, n_deg = n_deg, n_up = n_up, n_dn = n_dn,
                      n_panel = nrow(pt), comp = comp,
                      vista = as.character(packageVersion("VISTA")))

BUNDLE <- file.path(tempdir(), "cpd_report_bundle.rds")
saveRDS(B, BUNDLE)
QMD <- file.path(OUT_DIR, "CPD_KO_vs_WT_report.qmd")
file.copy(file.path(DATA_DIR, "cpd_report_template.qmd"), QMD, overwrite = TRUE)
file.copy(file.path(DATA_DIR, "report.css"), file.path(OUT_DIR, "report.css"), overwrite = TRUE)
quarto::quarto_render(input = QMD, output_file = "CPD_KO_vs_WT_report.html",
                      execute_params = list(bundle_rds = normalizePath(BUNDLE, winslash = "/")))
html <- file.path(OUT_DIR, "CPD_KO_vs_WT_report.html")
cat(sprintf("\nDone. %s (%.1f MB)\n", html, file.size(html)/1024^2))
setwd(old_wd)
