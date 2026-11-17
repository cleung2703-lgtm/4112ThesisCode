DATA_DIR <- "."
OUT_DIR  <- DATA_DIR

if (!file.exists(file.path(DATA_DIR, "MEF_CPD_KO_feature_count.txt"))) {
  if (file.exists("MEF_CPD_KO_feature_count.txt")) {
    DATA_DIR <- "."; OUT_DIR <- "."
    message("Using current working directory: ", normalizePath("."))
  } else stop("Cannot find MEF_CPD_KO_feature_count.txt\n  Looked in: ", DATA_DIR,
              "\n  and in   : ", normalizePath(".", mustWork = FALSE),
              "\n\nFix: set DATA_DIR at the top, or run setwd(\"<folder>\").")
}
cat("Data folder:", normalizePath(DATA_DIR), "\n")

TOP_N        <- 15
PADJ_GSEA    <- 0.05
GSEA_SEED    <- 42
N_PERM       <- 100000

PADJ_CUTOFF <- 0.05; LOG2FC_CUTOFF <- 1
VISTA_RDS <- file.path(OUT_DIR, "vista_object_shared.rds")
GSEA_RDS  <- file.path(OUT_DIR, "vista_gsea_shared.rds")
REBUILD_VISTA <- FALSE; REBUILD_GSEA <- FALSE
COL_UP <- "#D73027"; COL_DOWN <- "#2166AC"

suppressPackageStartupMessages({
  library(VISTA); library(org.Mm.eg.db); library(msigdbr); library(SummarizedExperiment)
  library(dplyr); library(tibble); library(ggplot2)
  library(plotly); library(htmlwidgets)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; slice_head <- dplyr::slice_head; count <- dplyr::count
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

if (Sys.which("pandoc") == "" && !nzchar(Sys.getenv("RSTUDIO_PANDOC"))) {
  .c <- c("C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
          "C:/Program Files/RStudio/bin/quarto/bin/tools",
          "C:/Program Files/RStudio/bin/pandoc")
  .h <- .c[file.exists(file.path(.c, "pandoc.exe"))]
  if (length(.h)) Sys.setenv(RSTUDIO_PANDOC = .h[1])
}
SELFCONTAINED <- nzchar(Sys.which("pandoc")) || nzchar(Sys.getenv("RSTUDIO_PANDOC"))

.label_text <- function(l, fallback = NULL) {
  if (is.null(l)) return(fallback)
  if (is.character(l) && length(l) == 1) return(l)
  if (is.expression(l) && length(l) >= 1) l <- l[[1]]
  txt <- paste(deparse(l), collapse = " ")
  txt <- gsub('"', "", txt); txt <- gsub("^~|~", " ", txt)
  txt <- gsub("\\[|\\]", "", txt); txt <- gsub("\\s+", " ", trimws(txt))
  if (!nzchar(txt)) fallback else txt
}
.scale_title <- function(p, aes_names) {
  for (s in p$scales$scales)
    if (any(s$aesthetics %in% aes_names)) return(.label_text(s$name))
  NULL
}
.plain_labels <- function(p) {
  for (n in names(p$labels)) {
    l <- p$labels[[n]]
    if (is.null(l) || (is.character(l) && length(l) == 1)) next
    p$labels[[n]] <- .label_text(l)
  }
  p
}
.add_size_legend <- function(pl, size_values, title = "Gene count") {
  if (is.null(size_values)) return(pl)
  px <- pl$x$data[[1]]$marker$size
  if (is.null(px) || length(px) != length(size_values)) return(pl)
  df <- unique(data.frame(val = as.numeric(size_values), px = as.numeric(px)))
  df <- df[is.finite(df$val) & is.finite(df$px), ]
  if (!nrow(df)) return(pl)
  df <- df[order(df$val), ]
  idx <- unique(round(seq(1, nrow(df), length.out = min(4, nrow(df)))))
  for (i in idx) {
    pl$x$data[[length(pl$x$data) + 1]] <- list(
      x = list(NULL), y = list(NULL),
      type = "scatter", mode = "markers",
      name = as.character(round(df$val[i])),
      legendgroup = "sizekey", showlegend = TRUE, hoverinfo = "skip",
      marker = list(size = df$px[i], color = "grey45",
                    line = list(width = 0.5, color = "grey25")),
      xaxis = "x", yaxis = "y")
  }
  plotly::layout(pl, showlegend = TRUE,
                 legend = list(title = list(text = title), itemsizing = "constant",
                               x = 1.02, y = 0.35, xanchor = "left",
                               bgcolor = "rgba(255,255,255,0)"))
}

save_widget_html <- function(p, file, title, size_title = NULL) {
  if (inherits(p, "ggplot")) {
    cb <- .scale_title(p, c("colour", "color", "fill")) %||%
          .label_text(p$labels$colour %||% p$labels$fill)
    if (is.null(size_title)) size_title <- .scale_title(p, "size") %||% "Gene count"
    sv <- if (".plot_size" %in% colnames(p$data)) p$data$.plot_size else
          if ("Count" %in% colnames(p$data)) p$data$Count else NULL
    p <- ggplotly(.plain_labels(p))
    if (!is.null(cb) && nzchar(cb))
      for (i in seq_along(p$x$data)) {
        mk <- p$x$data[[i]]$marker
        if (!is.null(mk) && !is.null(mk$colorbar))
          p$x$data[[i]]$marker$colorbar$title <- list(text = cb, side = "right")
      }
    p <- .add_size_legend(p, sv, size_title)
  }
  out <- file.path(OUT_DIR, file)
  saveWidget(plotly::config(p, displaylogo = FALSE), out,
             selfcontained = SELFCONTAINED, title = title)
  cat("Wrote:", out, "\n"); invisible(out)
}

if (!REBUILD_VISTA && file.exists(VISTA_RDS)) {
  vista <- readRDS(VISTA_RDS); cat("Loaded cached VISTA object\n")
} else {
  cm <- read.delim(file.path(DATA_DIR, "MEF_CPD_KO_feature_count.txt"),
                   row.names = "gene_id", check.names = FALSE) |> as.matrix()
  sm <- data.frame(sample_names = colnames(cm),
                   group = ifelse(grepl("^CPD-KO", colnames(cm)), "KO", "WT"))
  stopifnot(sum(sm$group == "KO") == 3, sum(sm$group == "WT") == 2)
  mi <- match_vista_inputs(read_vista_counts(cm, format = "matrix"),
                           read_vista_metadata(sm))
  vista <- create_vista(counts = mi$counts, sample_info = mi$sample_info,
    column_geneid = mi$column_geneid, group_column = "group",
    group_numerator = "KO", group_denominator = "WT", method = "deseq2",
    log2fc_cutoff = LOG2FC_CUTOFF, pval_cutoff = PADJ_CUTOFF,
    p_value_type = "padj", min_counts = 10, min_replicates = 1)
  vista <- set_rowdata(vista, orgdb = org.Mm.eg.db,
                       columns = c("SYMBOL", "GENENAME", "ENTREZID"),
                       keytype = "ENSEMBL")
  saveRDS(vista, VISTA_RDS); cat("Built and cached VISTA object\n")
}
comp <- names(comparisons(vista))[1]
cat("Comparison:", comp, "\n")

if (!REBUILD_GSEA && file.exists(GSEA_RDS)) {
  gsea_results <- readRDS(GSEA_RDS); cat("Loaded cached GSEA\n")
} else {
  set.seed(GSEA_SEED)
  gsea_results <- get_gsea(
    vista, sample_comparison = comp, set_type = "msigdb",
    from_type = "ENSEMBL",
    msigdb_category = "H",
    orgdb = org.Mm.eg.db, species = "Mus musculus",
    pvalueCutoff = 1, nPermSimple = N_PERM, eps = 0
  )
  saveRDS(gsea_results, GSEA_RDS)
  cat("Built and cached GSEA (seed", GSEA_SEED, ")\n")
}

gsea_table <- as.data.frame(gsea_results$enrich@result) |> as_tibble() |>
  select(any_of(c("ID", "Description", "NES", "pvalue", "p.adjust", "setSize"))) |>
  arrange(p.adjust)

write.csv(gsea_table, file.path(OUT_DIR, "05_gsea_results_table.csv"), row.names = FALSE)
cat(sprintf("\n%d Hallmark sets returned | %d significant at padj < %g\n",
            nrow(gsea_table), sum(gsea_table$p.adjust < PADJ_GSEA, na.rm = TRUE),
            PADJ_GSEA))

hm_all <- suppressMessages(msigdbr(species = "Mus musculus", collection = "H"))
missing <- setdiff(sort(unique(hm_all$gs_name)), gsea_table$ID)
if (length(missing)) {
  cat("\n!! ", length(missing), " Hallmark set(s) returned NO result and were dropped:\n", sep = "")
  for (m in missing) cat("     ", m, "\n")
  cat("   Cause: fgsea could not compute a p-value (unbalanced gene-level\n",
      "   statistics). Raising N_PERM above ", N_PERM, " may not help.\n",
      "   Report these as 'not evaluable under VISTA's GSEA', not as 'not significant'.\n", sep = "")
} else cat("\nAll 50 Hallmark sets evaluated.\n")

cat("\nTop pathways:\n")
print(gsea_table |> slice_head(n = 20) |>
        mutate(Description = gsub("^HALLMARK_", "", Description),
               NES = round(NES, 3), p.adjust = signif(p.adjust, 3)) |>
        select(Description, NES, p.adjust, setSize) |> as.data.frame(), row.names = FALSE)

p_enr <- get_enrichment_plot(gsea_results$enrich, top_n = TOP_N,
                             title = "Hallmark GSEA: CPD KO vs WT MEFs")
ggsave(file.path(OUT_DIR, "05_gsea_dotplot.png"), p_enr, width = 11, height = 8, dpi = 300)
save_widget_html(p_enr, "05_gsea_dotplot.html", "Hallmark GSEA - CPD KO vs WT")

sig <- gsea_table |> filter(p.adjust < PADJ_GSEA)
if (nrow(sig)) {
  bar <- sig |> arrange(desc(abs(NES))) |> slice_head(n = TOP_N) |>
    mutate(label = gsub("_", " ", gsub("^HALLMARK_", "", Description)),
           Direction = ifelse(NES > 0, "Enriched in KO", "Enriched in WT")) |>
    arrange(NES) |> mutate(label = factor(label, levels = label))

  p_bar <- ggplot(bar, aes(NES, label, fill = Direction)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("Enriched in KO" = COL_UP,
                                 "Enriched in WT" = COL_DOWN)) +
    labs(title = "Hallmark GSEA: CPD KO vs WT MEFs",
         subtitle = sprintf("%d pathways at padj < %g, ranked by |NES|",
                            nrow(bar), PADJ_GSEA),
         x = "Normalised enrichment score (NES)", y = NULL) +
    theme_minimal(base_size = 12)

  ggsave(file.path(OUT_DIR, "05_gsea_bar.png"), p_bar, width = 9, height = 6, dpi = 300)
  save_widget_html(p_bar, "05_gsea_bar.html", "Hallmark GSEA NES - CPD KO vs WT")
} else cat("\nNothing significant at padj <", PADJ_GSEA, "- skipping the bar chart.\n")

pg <- get_pathway_genes(gsea_results, pathways = sig$Description,
                        return_type = "long") |> as_tibble()
cat("\nLeading-edge genes per significant pathway:\n")
print(as.data.frame(count(pg, pathway)), row.names = FALSE)

rd  <- as.data.frame(SummarizedExperiment::rowData(vista))
sym <- setNames(as.character(rd$SYMBOL), rownames(rd))

de <- as.data.frame(comparisons(vista)[[comp]]) |>
  select(gene_id, log2fc, padj, regulation)

gene_col <- intersect(c("gene", "gene_id", "geneID"), colnames(pg))[1]
pg <- pg |>
  mutate(ensembl = .data[[gene_col]],
         symbol  = unname(sym[ensembl]),
         symbol  = ifelse(is.na(symbol) | symbol == "", ensembl, symbol)) |>
  left_join(de, by = c("ensembl" = "gene_id"))

n_unmapped <- sum(pg$symbol == pg$ensembl)
cat(sprintf("Symbol mapping: %d of %d leading-edge entries had no symbol (kept as Ensembl)\n",
            n_unmapped, nrow(pg)))

pg |>
  mutate(pathway = gsub("^HALLMARK_", "", pathway)) |>
  select(pathway, gene = symbol, ensembl, log2FC = log2fc, padj, regulation) |>
  arrange(pathway, desc(abs(log2FC))) |>
  write.csv(file.path(OUT_DIR, "05_gsea_leading_edge_by_gene.csv"), row.names = FALSE)

pg |> group_by(pathway) |>
  summarise(n_genes = dplyr::n(),
            genes = paste(sort(unique(symbol)), collapse = ", "), .groups = "drop") |>
  left_join(gsea_table |> select(Description, NES, p.adjust),
            by = c("pathway" = "Description")) |>
  arrange(p.adjust) |>
  write.csv(file.path(OUT_DIR, "05_gsea_leading_edge_genes.csv"), row.names = FALSE)

hub <- pg |>
  group_by(symbol) |>
  summarise(n_pathways = dplyr::n_distinct(pathway),
            pathways = paste(sort(gsub("^HALLMARK_", "", unique(pathway))), collapse = "; "),
            log2FC = dplyr::first(log2fc), padj = dplyr::first(padj),
            regulation = dplyr::first(regulation), .groups = "drop") |>
  arrange(desc(n_pathways), padj)
write.csv(hub, file.path(OUT_DIR, "05_gsea_genes_by_pathway_count.csv"), row.names = FALSE)

cat("Wrote: 05_gsea_leading_edge_by_gene.csv (one row per gene per pathway)\n")
cat("Wrote: 05_gsea_leading_edge_genes.csv   (one row per pathway)\n")
cat("Wrote: 05_gsea_genes_by_pathway_count.csv (one row per gene, hub view)\n")

cat("\nLeading-edge genes in the most pathways:\n")
print(hub |> slice_head(n = 15) |>
        mutate(log2FC = round(log2FC, 2), padj = signif(padj, 3)) |>
        select(symbol, n_pathways, log2FC, padj, regulation) |>
        as.data.frame(), row.names = FALSE)

cat("\nNOTE ON METHOD: VISTA's get_gsea() ranks genes by log2FC (verified: the\n",
    "geneList range matches the log2fc column exactly). It does not expose a\n",
    "ranking-metric argument, so the Wald-statistic arm of the sensitivity\n",
    "analysis cannot be run through VISTA. The laboratory record therefore notes\n",
    "fgsea-direct as the method for that comparison. EMT surfaces only under\n",
    "Wald ranking, so do not expect it in this table.\n", sep = "")

cat("\nDone.\n")
