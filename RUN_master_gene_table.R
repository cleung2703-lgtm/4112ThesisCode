DATA_DIR <- "."
OUT_DIR  <- DATA_DIR

if (!file.exists(file.path(DATA_DIR, "MEF_CPD_KO_feature_count.txt"))) {
  if (file.exists("MEF_CPD_KO_feature_count.txt")) { DATA_DIR <- "."; OUT_DIR <- "."
  } else stop("Cannot find MEF_CPD_KO_feature_count.txt\n  Looked in: ", DATA_DIR,
              "\n  Fix: set DATA_DIR at the top, or setwd(\"<folder>\").")
}
cat("Data folder:", normalizePath(DATA_DIR), "\n")

PADJ_CUTOFF <- 0.05; LOG2FC_CUTOFF <- 1
N_SHORTLIST <- 30
GO_PADJ     <- 0.05

ROBUST_PATHWAYS <- c("MYOGENESIS", "IL2_STAT5_SIGNALING",
                     "EPITHELIAL_MESENCHYMAL_TRANSITION")
THESIS_A <- c("Smpd3","Pitx2","Epas1","Fgf18","Chd7","Rbp4","Snap25","Ednrb",
              "Hmga2","Oas1a","Oas1g","Dio2","Hmga1")

for (p in c("DT", "writexl")) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, type = "binary")
suppressPackageStartupMessages({
  library(VISTA); library(SummarizedExperiment); library(org.Mm.eg.db)
  library(dplyr); library(tibble); library(readxl); library(writexl)
  library(DT); library(htmlwidgets)
})
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; count <- dplyr::count; slice_head <- dplyr::slice_head

if (Sys.which("pandoc") == "" && !nzchar(Sys.getenv("RSTUDIO_PANDOC"))) {
  .c <- c("C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
          "C:/Program Files/RStudio/bin/quarto/bin/tools")
  .h <- .c[file.exists(file.path(.c, "pandoc.exe"))]
  if (length(.h)) Sys.setenv(RSTUDIO_PANDOC = .h[1])
}
SC <- nzchar(Sys.which("pandoc")) || nzchar(Sys.getenv("RSTUDIO_PANDOC"))

VISTA_RDS <- file.path(OUT_DIR, "vista_object_shared.rds")
if (!file.exists(VISTA_RDS))
  stop("vista_object_shared.rds not found. Run RUN_volcano_standalone.R or\n",
       "RUN_gsea_standalone.R first to build and cache the VISTA object.")
vista <- readRDS(VISTA_RDS)
comp  <- names(comparisons(vista))[1]
cat("Loaded VISTA object |", comp, "\n")

de <- as.data.frame(comparisons(vista)[[comp]]) |>
  select(ensembl = gene_id, baseMean, log2FC = log2fc, padj, regulation)

rd  <- as.data.frame(rowData(vista))
de$gene     <- unname(as.character(rd$SYMBOL)[match(de$ensembl, rownames(rd))])
de$gene_name<- unname(as.character(rd$GENENAME)[match(de$ensembl, rownames(rd))])
de$gene     <- ifelse(is.na(de$gene) | de$gene == "", de$ensembl, de$gene)

nc  <- assay(vista, "norm_counts")
grp <- as.character(colData(vista)$group)
ko  <- nc[, grp == "KO", drop = FALSE]; wt <- nc[, grp == "WT", drop = FALSE]
cat(sprintf("Per-sample norm_counts: %d KO, %d WT\n", ncol(ko), ncol(wt)))

samp <- data.frame(
  ensembl = rownames(nc),
  KO_mean = rowMeans(ko), WT_mean = rowMeans(wt),
  KO_min = apply(ko, 1, min), KO_max = apply(ko, 1, max),
  WT_min = apply(wt, 1, min), WT_max = apply(wt, 1, max),
  KO_cv  = apply(ko, 1, sd) / pmax(rowMeans(ko), 1e-9)
) |>
  mutate(separation = ifelse(KO_min > WT_max | KO_max < WT_min,
                             "Clean", "Overlapping"))
samp <- cbind(samp, as.data.frame(nc))

go_path <- file.path(DATA_DIR, "03_go_bp_results.csv")
if (file.exists(go_path)) {
  go <- read.csv(go_path, stringsAsFactors = FALSE) |> filter(p.adjust < GO_PADJ)
  go_long <- go |>
    select(Description, regulation_go = regulation, geneID) |>
    tidyr::separate_rows(geneID, sep = "/") |>
    rename(ensembl = geneID) |>
    filter(ensembl != "")
  go_sum <- go_long |> group_by(ensembl) |>
    summarise(n_go_terms = dplyr::n_distinct(Description),
              go_terms = paste(utils::head(sort(unique(Description)), 8), collapse = "; "),
              .groups = "drop")
  cat(sprintf("GO: %d significant terms, %d genes with membership\n",
              dplyr::n_distinct(go$Description), nrow(go_sum)))
} else {
  go_sum <- tibble(ensembl = character(), n_go_terms = integer(), go_terms = character())
  cat("GO results not found - run RUN_go_enrichment_standalone.R for those columns\n")
}

gs_path <- file.path(DATA_DIR, "05_gsea_leading_edge_by_gene.csv")
if (file.exists(gs_path)) {
  gs <- read.csv(gs_path, stringsAsFactors = FALSE)
  gs_sum <- gs |> group_by(ensembl) |>
    summarise(n_gsea_pathways = dplyr::n_distinct(pathway),
              gsea_pathways = paste(sort(unique(pathway)), collapse = "; "),
              .groups = "drop") |>
    mutate(robust_pathway = grepl(paste(ROBUST_PATHWAYS, collapse = "|"), gsea_pathways))
  cat(sprintf("GSEA: %d genes in a significant leading edge\n", nrow(gs_sum)))
} else {
  gs_sum <- tibble(ensembl = character(), n_gsea_pathways = integer(),
                   gsea_pathways = character(), robust_pathway = logical())
  cat("GSEA leading edge not found - run RUN_gsea_standalone.R for those columns\n")
}

pr_path <- file.path(DATA_DIR, "Cpd_qPCR_primers.xlsx")
if (file.exists(pr_path)) {
  pr <- suppressMessages(read_excel(pr_path, col_names = FALSE))
  names(pr)[1:3] <- c("primer", "seq", "amplicon")
  primers <- pr |>
    filter(!is.na(primer), grepl("-5'", primer), !is.na(amplicon)) |>
    transmute(gene_upper = toupper(gsub("^hm|-5'.*$", "", primer)),
              primer_amplicon = gsub("\\s+", " ", trimws(amplicon)))
  cat(sprintf("Primers: %d validated genes\n", nrow(primers)))
} else {
  primers <- tibble(gene_upper = character(), primer_amplicon = character())
  cat("Primer file not found\n")
}

master <- de |>
  left_join(samp,   by = "ensembl") |>
  left_join(go_sum, by = "ensembl") |>
  left_join(gs_sum, by = "ensembl") |>
  mutate(gene_upper = toupper(gene)) |>
  left_join(primers, by = "gene_upper") |>
  select(-gene_upper) |>
  mutate(
    n_go_terms      = tidyr::replace_na(n_go_terms, 0L),
    n_gsea_pathways = tidyr::replace_na(n_gsea_pathways, 0L),
    robust_pathway  = tidyr::replace_na(robust_pathway, FALSE),
    is_DEG    = regulation %in% c("Up", "Down"),
    thesis_a  = gene %in% THESIS_A,
    has_primer = !is.na(primer_amplicon),
    score = (2 * is_DEG) +
            (separation == "Clean" & !is.na(separation)) +
            (n_gsea_pathways > 0) + robust_pathway +
            (n_go_terms >= 3) + has_primer + thesis_a
  ) |>
  arrange(desc(score), padj)

cat(sprintf("\nMaster table: %d genes | %d DEGs | %d with clean separation\n",
            nrow(master), sum(master$is_DEG),
            sum(master$is_DEG & master$separation == "Clean", na.rm = TRUE)))

shortlist <- master |> filter(is_DEG) |> arrange(desc(score), padj) |>
  slice_head(n = N_SHORTLIST)

cat(sprintf("\n--- Top %d DEGs by evidence score ---\n", N_SHORTLIST))
print(shortlist |>
        mutate(log2FC = round(log2FC, 2), padj = signif(padj, 3)) |>
        select(gene, log2FC, padj, regulation, separation,
               n_go_terms, n_gsea_pathways, robust_pathway,
               primer = primer_amplicon, score) |>
        as.data.frame(), row.names = FALSE)

cols <- c("gene","gene_name","ensembl","regulation","log2FC","padj","baseMean",
          "is_DEG","separation","KO_mean","WT_mean","KO_cv",
          "n_go_terms","n_gsea_pathways","robust_pathway","thesis_a",
          "primer_amplicon","score","go_terms","gsea_pathways",
          colnames(nc))

tidy <- function(d) d |> select(any_of(cols)) |>
  mutate(across(any_of(c("log2FC","KO_mean","WT_mean","KO_cv","baseMean")), ~round(.x, 3)),
         across(any_of("padj"), ~signif(.x, 3)),
         across(any_of(colnames(nc)), ~round(.x, 1)))

write.csv(tidy(shortlist), file.path(OUT_DIR, "09_candidate_shortlist.csv"), row.names = FALSE)
write_xlsx(list(
  shortlist   = tidy(shortlist),
  all_DEGs    = tidy(master |> filter(is_DEG)),
  all_genes   = tidy(master),
  gsea_genes  = if (exists("gs")) gs else tibble(),
  go_terms    = if (exists("go")) go |> select(any_of(c("Description","regulation",
                    "GeneRatio","p.adjust","Count"))) else tibble()
), file.path(OUT_DIR, "09_master_gene_table.xlsx"))
cat("\nWrote: 09_master_gene_table.xlsx\n       09_candidate_shortlist.csv\n")

html_rows <- master |> filter(is_DEG | n_gsea_pathways > 0 | n_go_terms > 0 | has_primer)
dt <- datatable(tidy(html_rows), filter = "top", rownames = FALSE,
  extensions = "Buttons", caption = paste0(
    "CPD KO vs WT - ", nrow(html_rows), " genes that are DEGs, in a GO term, ",
    "in a GSEA leading edge, or have a validated primer. ",
    "All statistics from the shared VISTA object."),
  options = list(pageLength = 30, scrollX = TRUE, dom = "Bfrtip",
                 buttons = c("copy","csv","excel")))
saveWidget(dt, file.path(OUT_DIR, "09_master_gene_table.html"),
           selfcontained = SC, title = "CPD KO vs WT - master gene table")
cat("Wrote: 09_master_gene_table.html (", nrow(html_rows), "rows )\n")

cat("\nScore components: 2 x DEG, +1 clean sample separation, +1 in a GSEA\n",
    "leading edge, +1 in a dependable pathway, +1 in >=3 GO terms, +1 validated\n",
    "primer, +1 Thesis A gene. Max 8. Every component is its own column, so\n",
    "the total can be re-weighted or filtered instead.\n", sep = "")
cat("\nDone.\n")
