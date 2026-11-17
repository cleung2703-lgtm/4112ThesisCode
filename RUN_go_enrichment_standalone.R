DATA_DIR <- "."
OUT_DIR  <- DATA_DIR

if (!file.exists(file.path(DATA_DIR, "MEF_CPD_KO_feature_count.txt"))) {
  if (file.exists("MEF_CPD_KO_feature_count.txt")) {
    DATA_DIR <- "."; OUT_DIR <- "."
    message("Using the current working directory: ", normalizePath("."))
  } else {
    stop("Cannot find MEF_CPD_KO_feature_count.txt\n",
         "  Looked in : ", DATA_DIR, "\n",
         "  and in    : ", normalizePath(".", mustWork = FALSE), "\n\n",
         "Fix: set DATA_DIR at the top of this script, or run setwd(\"<folder>\").")
  }
}
cat("Data folder:", normalizePath(DATA_DIR), "\n")

ONT   <- "BP"
TOP_N <- 15
REGS  <- c("Up", "Down", "Both")

PADJ_CUTOFF   <- 0.05
LOG2FC_CUTOFF <- 1
VISTA_RDS     <- file.path(OUT_DIR, "vista_object_shared.rds")
REBUILD_VISTA <- FALSE

suppressPackageStartupMessages({
  library(VISTA); library(org.Mm.eg.db)
  library(dplyr); library(tibble); library(ggplot2)
  library(plotly); library(htmlwidgets)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; slice_head <- dplyr::slice_head
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

.plain_labels <- function(p) {
  if (!inherits(p, "ggplot")) return(p)
  for (n in names(p$labels)) {
    l <- p$labels[[n]]
    if (is.null(l) || (is.character(l) && length(l) == 1)) next
    p$labels[[n]] <- .label_text(l)
  }
  p
}

.scale_title <- function(p, aes_names) {
  for (s in p$scales$scales)
    if (any(s$aesthetics %in% aes_names)) return(.label_text(s$name))
  NULL
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
  cat("Wrote:", out, "
"); invisible(out)
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
  vista <- create_vista(
    counts = mi$counts, sample_info = mi$sample_info,
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
print(deg_summary(vista))

go_objects <- list()
go_tables  <- list()

for (reg in REGS) {
  cat("\n=== GO", ONT, "-", reg, "===\n")

  res <- try(get_go_enrichment(
    x = vista, sample_comparison = comp, regulation = reg, ont = ONT,
    orgdb = org.Mm.eg.db, species = "Mus musculus", from_type = "ENSEMBL"
  ), silent = TRUE)

  if (inherits(res, "try-error")) {
    cat("FAILED:", conditionMessage(attr(res, "condition")), "\n"); next
  }

  obj <- if (!is.null(res[[reg]])) res[[reg]] else res[[1]]
  if (is.null(obj)) { cat("no result object returned\n"); next }

  d <- try(as.data.frame(obj@result), silent = TRUE)
  if (inherits(d, "try-error") || !nrow(d)) { cat("empty result\n"); next }

  n_sig <- sum(d$p.adjust < PADJ_CUTOFF, na.rm = TRUE)
  cat(sprintf("%d GO terms, %d significant at p.adjust < %g\n",
              nrow(d), n_sig, PADJ_CUTOFF))
  if (n_sig == 0) {
    cat("Nothing significant, so the dotplot for", reg, "would be empty. Skipping.\n")
    next
  }
  print(utils::head(d[order(d$p.adjust), c("Description", "GeneRatio", "p.adjust", "Count")], 5))

  go_objects[[reg]] <- obj
  d$regulation <- reg
  go_tables[[reg]] <- d

  p <- try(get_enrichment_plot(x = obj, top_n = TOP_N,
                               title = sprintf("GO %s enrichment - %s-regulated", ONT, reg)),
           silent = TRUE)
  if (inherits(p, "try-error")) {
    cat("get_enrichment_plot failed:", conditionMessage(attr(p, "condition")), "\n"); next
  }

  ggsave(file.path(OUT_DIR, sprintf("03_go_%s_%s_dotplot.png", tolower(ONT), tolower(reg))),
         p, width = 11, height = 8, dpi = 300)
  save_widget_html(p, sprintf("03_go_%s_%s_dotplot.html", tolower(ONT), tolower(reg)),
                   sprintf("GO %s enrichment - %s", ONT, reg),
                   size_title = "Gene count")
}

if (!length(go_objects))
  stop("No GO enrichment produced a usable result. Read the messages above.")

full <- bind_rows(go_tables)
write.csv(full, file.path(OUT_DIR, sprintf("03_go_%s_results.csv", tolower(ONT))),
          row.names = FALSE)
cat("\nWrote:", file.path(OUT_DIR, sprintf("03_go_%s_results.csv", tolower(ONT))), "\n")

cat("\nSummary of significant terms per direction:\n")
print(as.data.frame(dplyr::count(dplyr::filter(full, p.adjust < PADJ_CUTOFF),
                                 regulation)))

cat("\nDone.\n")
