DATA_DIR <- "."
OUT_DIR  <- "."

PADJ_CUTOFF   <- 0.05
LOG2FC_CUTOFF <- 1
COL_UP <- "#D73027"; COL_DOWN <- "#2166AC"; COL_NS <- "#BDBDBD"

VISTA_RDS     <- file.path(OUT_DIR, "vista_object_shared.rds")
GSEA_RDS      <- file.path(OUT_DIR, "vista_gsea_shared.rds")
GSEA_SEED     <- 42
REBUILD_VISTA <- FALSE
REBUILD_GSEA  <- FALSE

suppressPackageStartupMessages({
  library(VISTA); library(org.Mm.eg.db)
  library(dplyr); library(tibble); library(ggplot2)
  library(plotly); library(htmlwidgets); library(base64enc)
})

select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
mutate <- dplyr::mutate; arrange <- dplyr::arrange; count <- dplyr::count
group_by <- dplyr::group_by; summarise <- dplyr::summarise
pull <- dplyr::pull; slice_head <- dplyr::slice_head
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

if (Sys.which("pandoc") == "" && !nzchar(Sys.getenv("RSTUDIO_PANDOC"))) {
  .c <- c("C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
          "C:/Program Files/RStudio/bin/quarto/bin/tools",
          "C:/Program Files/RStudio/bin/pandoc")
  .h <- .c[file.exists(file.path(.c, "pandoc.exe"))]
  if (length(.h)) Sys.setenv(RSTUDIO_PANDOC = .h[1])
}
SELFCONTAINED <- nzchar(Sys.which("pandoc")) || nzchar(Sys.getenv("RSTUDIO_PANDOC"))

.plain_labels <- function(p) {
  if (!inherits(p, "ggplot")) return(p)
  for (n in names(p$labels)) {
    l <- p$labels[[n]]
    if (is.null(l) || (is.character(l) && length(l) == 1)) next
    txt <- paste(deparse(l), collapse = " ")
    txt <- gsub('"', "", txt)
    txt <- gsub("^~|~", " ", txt)
    txt <- gsub("\\[|\\]", "", txt)
    txt <- gsub("\\s+", " ", trimws(txt))
    p$labels[[n]] <- txt
  }
  p
}

save_widget_html <- function(p, file, title, tooltip = "text") {
  if (inherits(p, "ggplot")) p <- ggplotly(.plain_labels(p), tooltip = tooltip)
  out <- file.path(OUT_DIR, file)
  saveWidget(plotly::config(p, displaylogo = FALSE), out,
             selfcontained = SELFCONTAINED, title = title)
  cat("Wrote:", out, "\n"); invisible(out)
}

save_static_html <- function(draw_fn, file, title, subtitle = NULL,
                             width = 10, height = 12, res = 200) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = width, height = height, units = "in", res = res)
  on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
  draw_fn(); dev.off()
  out <- file.path(OUT_DIR, file)
  writeLines(c(
    "<!doctype html><html><head><meta charset='utf-8'>",
    sprintf("<title>%s</title>", title),
    "<style>body{margin:0;padding:24px;background:#fff;color:#111;",
    "font-family:-apple-system,Segoe UI,Roboto,Arial,sans-serif}",
    "h1{font-size:19px;margin:0 0 4px}p.sub{color:#666;font-size:13px;margin:0 0 18px}",
    "img{max-width:100%;height:auto;display:block;border:1px solid #e5e5e5;border-radius:4px}",
    "@media(prefers-color-scheme:dark){body{background:#111;color:#eee}",
    "p.sub{color:#aaa}img{background:#fff;border-color:#333}}</style></head><body>",
    sprintf("<h1>%s</h1>", title),
    if (!is.null(subtitle)) sprintf("<p class='sub'>%s</p>", subtitle) else "",
    sprintf("<img src='data:image/png;base64,%s' alt='%s'>",
            base64enc::base64encode(tmp), title),
    "</body></html>"), out)
  cat("Wrote:", out, "\n"); invisible(out)
}

build_vista_shared <- function() {
  cm <- read.delim(file.path(DATA_DIR, "MEF_CPD_KO_feature_count.txt"),
                   row.names = "gene_id", check.names = FALSE) %>% as.matrix()
  sm <- data.frame(sample_names = colnames(cm),
                   group = ifelse(grepl("^CPD-KO", colnames(cm)), "KO", "WT"))
  stopifnot(sum(sm$group == "KO") == 3, sum(sm$group == "WT") == 2)

  mi <- match_vista_inputs(read_vista_counts(cm, format = "matrix"),
                           read_vista_metadata(sm))
  v <- create_vista(
    counts = mi$counts, sample_info = mi$sample_info,
    column_geneid = mi$column_geneid, group_column = "group",
    group_numerator = "KO", group_denominator = "WT", method = "deseq2",
    log2fc_cutoff = LOG2FC_CUTOFF, pval_cutoff = PADJ_CUTOFF,
    p_value_type = "padj",
    min_counts = 10, min_replicates = 1
  )
  set_rowdata(v, orgdb = org.Mm.eg.db,
              columns = c("SYMBOL", "GENENAME", "ENTREZID"), keytype = "ENSEMBL")
}

if (!REBUILD_VISTA && file.exists(VISTA_RDS)) {
  vista <- readRDS(VISTA_RDS); cat("Loaded cached VISTA object\n")
} else {
  vista <- build_vista_shared(); saveRDS(vista, VISTA_RDS)
  cat("Built and cached VISTA object\n")
}

comp <- names(comparisons(vista))[1]
cat("Comparison:", comp, "\n")
print(deg_summary(vista))
cat("\nNOTE: expect ~84 up / 41 down, not 97/59. The 97/59 figure is\n",
    "the LAB's call from their pre-filtered 13,080-gene file; this re-runs DESeq2\n",
    "on all 48,709 rows, so a different number of genes enters the BH correction.\n")

get_gsea_cached <- function() {
  if (!REBUILD_GSEA && file.exists(GSEA_RDS)) {
    cat("Loaded cached GSEA\n"); return(readRDS(GSEA_RDS))
  }
  set.seed(GSEA_SEED)
  g <- get_gsea(vista, sample_comparison = comp, set_type = "msigdb",
                from_type = "ENSEMBL", msigdb_category = "H",
                orgdb = org.Mm.eg.db, species = "Mus musculus")
  saveRDS(g, GSEA_RDS); cat("Built and cached GSEA (seed", GSEA_SEED, ")\n")
  g
}

get_de_table <- function() {
  d <- as.data.frame(comparisons(vista)[[comp]])
  cat("VISTA DE table columns:", paste(colnames(d), collapse = ", "), "\n")
  pick <- function(pats) {
    for (p in pats) {
      h <- grep(p, colnames(d), ignore.case = TRUE, value = TRUE)
      if (length(h)) return(h[1])
    }
    NA_character_
  }
  ci <- pick(c("^gene_id$", "^gene$", "^ensembl"))
  cl <- pick(c("^log2fc$", "^log2foldchange$", "log2.*fold", "lfc"))
  cp <- pick(c("^padj$", "adj.*p", "p.*adj", "fdr"))
  cs <- pick(c("^stat$", "wald", "statistic"))
  cy <- pick(c("^symbol$", "gene_name"))
  if (any(is.na(c(ci, cl, cp))))
    stop("Could not find gene_id / log2FC / padj above; map them by hand.")
  tibble(gene_id = as.character(d[[ci]]),
         log2fc  = as.numeric(d[[cl]]),
         padj    = as.numeric(d[[cp]]),
         stat    = if (!is.na(cs)) as.numeric(d[[cs]]) else NA_real_,
         symbol  = if (!is.na(cy)) as.character(d[[cy]]) else NA_character_)
}
