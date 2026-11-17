# Carboxypeptidase D knockout MEF analysis code

R code for the Bachelor of Engineering (Honours) thesis *Transcriptional Consequences of Carboxypeptidase D Loss in Mouse Embryonic Fibroblasts* (School of Biomedical Engineering, 2026). The same scripts are printed in Appendix D of the thesis, with the section numbers given below.

The analysis compares three *Cpd* knockout and two wild-type mouse embryonic fibroblast lines by RNA-seq (differential expression, Gene Ontology and Hallmark enrichment, candidate gene scoring) and follows up a gene panel by qRT-PCR.

## Scripts

| Appendix | Script | What it does |
|---|---|---|
| D.3 | [`00_vista_common.R`](00_vista_common.R) | Shared configuration: file paths, the significance rule, the colour scheme and helper functions used by every analysis script. |
| D.4 | [`RUN_CPD_REPORT.R`](RUN_CPD_REPORT.R) | Produces the VISTA figures used in Chapter 4: quality control, MA plot, Gene Ontology and Hallmark enrichment, the chord diagrams and the panel heatmap. |
| D.5 | [`RUN_ALL_VISTA_FIGURES.R`](RUN_ALL_VISTA_FIGURES.R) | Generates the full VISTA figure set from the cached analysis object, including the figures not carried into the thesis. |
| D.6 | [`VOLCANO_BROKEN_AXIS.R`](VOLCANO_BROKEN_AXIS.R) | Builds the broken-axis volcano plot of Figure 4.2 by calling the VISTA volcano function twice over different vertical ranges and stacking the two panels. |
| D.7 | [`RUN_go_enrichment_standalone.R`](RUN_go_enrichment_standalone.R) | Gene Ontology over-representation for the combined, upregulated and downregulated gene lists. |
| D.8 | [`RUN_gsea_standalone.R`](RUN_gsea_standalone.R) | Rank-based gene set enrichment with fgsea, including the alternative ranking metrics and gene set versions compared in Section 4.4.3. |
| D.9 | [`RUN_master_gene_table.R`](RUN_master_gene_table.R) | Assembles the per-gene evidence table combining differential expression, Gene Ontology connectivity and enrichment membership. |
| D.10 | [`RUN_FINAL_20_TABLE.R`](RUN_FINAL_20_TABLE.R) | Applies the seven scoring criteria of Section 3.2.7 and produces the ranked candidate table (Table F.2). |
| D.11 | [`RUN_PROLIFERATION_SCORE.R`](RUN_PROLIFERATION_SCORE.R) | Computes the per-sample Hallmark scores of Figure 4.5 from the normalised counts, using gene set membership taken from MSigDB. |
| D.12 | [`RUN_GSEA_SENSITIVITY_FIGURE.R`](RUN_GSEA_SENSITIVITY_FIGURE.R) | Regenerates Figure 4.6 from the adjusted p-values produced by the sensitivity analysis of Section 3.2.6, comparing each Hallmark set under the two gene ranking metrics. |
| D.13 | [`RUN_DEG_EXCEL.R`](RUN_DEG_EXCEL.R) | Exports the differential expression results used for the pipeline comparison in Table F.1. |
| D.14 | [`qPCR/RUN_QPCR_ANALYSIS.R`](qPCR/RUN_QPCR_ANALYSIS.R) | The complete qRT-PCR pipeline: reads every instrument export, assigns samples by replicate position, applies the quality criteria of Section 3.7.2 and produces Table 4.6, Figure 4.8, Table E.5, Figure E.5 and Figures G.4 and G.5. |
| D.15 | [`RUN_QPCR_POWER.R`](RUN_QPCR_POWER.R) | Estimates the smallest change the qRT-PCR design could have detected (Figure 4.9). |
| D.16 | [`RUN_EXTERNAL_VALIDATION.R`](RUN_EXTERNAL_VALIDATION.R) | Tests the knockout transcriptome against pre-specified MSigDB perturbation signatures of IGF1R, PI3K, AKT and mTOR signalling (Figure 4.10 and Table E.6). |

## Which script produced each thesis figure and table

| Thesis item | Script |
|---|---|
| Table 3.1, library size and gene detection | RUN_DEG_EXCEL.R |
| Figure 4.1 and Figures G.1 and G.2, quality control | RUN_CPD_REPORT.R |
| Figure 4.2, broken-axis volcano | VOLCANO_BROKEN_AXIS.R |
| Table 4.1, most significant genes | RUN_DEG_EXCEL.R |
| Table 4.2, functional groups | prepared manually from the output of RUN_DEG_EXCEL.R |
| Figure 4.3 and Table 4.3, Gene Ontology enrichment and themes | RUN_CPD_REPORT.R and RUN_go_enrichment_standalone.R |
| Table 4.4 and Table 4.5, Hallmark results and confidence categories | RUN_gsea_standalone.R |
| Figure 4.4 and Figure G.3, Hallmark chord and enrichment | RUN_CPD_REPORT.R |
| Figure 4.5, per-sample proliferation scores | RUN_PROLIFERATION_SCORE.R |
| Figure 4.6, ranking metric sensitivity | RUN_GSEA_SENSITIVITY_FIGURE.R |
| Table 4.6 and Figure 4.8, qRT-PCR results | qPCR/RUN_QPCR_ANALYSIS.R |
| Figure 4.9, qRT-PCR detectability | RUN_QPCR_POWER.R |
| Figure 4.10 and Table E.6, independent signature test | RUN_EXTERNAL_VALIDATION.R |
| Table F.1, pipeline comparison | RUN_DEG_EXCEL.R |
| Table F.2, candidate scores | RUN_FINAL_20_TABLE.R, from RUN_master_gene_table.R |
| Table F.3 and Figure F.1, hub genes | RUN_go_enrichment_standalone.R and RUN_CPD_REPORT.R |
| Table E.5 and Figure E.5, failed assays and well viability | qPCR/RUN_QPCR_ANALYSIS.R |
| Figures G.4 and G.5, qRT-PCR supporting figures | qPCR/RUN_QPCR_ANALYSIS.R |
| Figures G.6 and G.7, supporting enrichment figures | RUN_ALL_VISTA_FIGURES.R |
| Figures B.1 and B.2, direction-specific enrichment | RUN_go_enrichment_standalone.R |

## Input data (not included)

The data belong to the host laboratory and are not in this repository. The scripts expect these files in the working directory:

- `MEF_CPD_KO_feature_count.txt`: raw count matrix
- `MEF_CPD_KO_de_genes.xlsx`: the laboratory DESeq2 results, used for the pipeline comparison
- `vista_object_shared.rds`: the cached VISTA analysis object, so every figure comes from the same run
- `qPCR/raw/`: the Mx3005P instrument exports

A few scripts also read small tables (for example the primer list or normalised counts) named at the top of each script.

## Running the scripts

Run every script from the folder that holds the data, so relative paths such as `qPCR/raw` resolve. Every script takes its data from `DATA_DIR`, set to the working directory on the first line. Run them from the folder holding the data, or edit that line to point elsewhere.

```bash
Rscript RUN_ALL_VISTA_FIGURES.R
Rscript qPCR/RUN_QPCR_ANALYSIS.R
```

Some scripts read files written by others:

- `RUN_ALL_VISTA_FIGURES.R` writes `vista_tables/de_results_all_genes.csv`, which `RUN_QPCR_POWER.R` and `RUN_EXTERNAL_VALIDATION.R` read.
- `RUN_go_enrichment_standalone.R` and `RUN_gsea_standalone.R` write the GO and GSEA tables that `RUN_master_gene_table.R` and `RUN_FINAL_20_TABLE.R` read.
- `qPCR/RUN_QPCR_ANALYSIS.R` writes `qPCR/results/qPCR_results.xlsx`, which `RUN_QPCR_POWER.R` reads.

## Software

The thesis results were generated with R 4.6.1, VISTA 1.1.5 (with DESeq2 1.52.0), clusterProfiler 4.20.0, org.Mm.eg.db 3.18.0, msigdbr 26.1.0, fgsea 1.38.0, ggplot2 4.0.3, EnhancedVolcano 1.30.0, ComplexHeatmap 2.28.0 and circlize 0.4.18.

```r
install.packages(c("BiocManager", "dplyr", "tidyr", "tibble", "readr", "readxl", "writexl",
                   "ggplot2", "ggrepel", "patchwork", "plotly", "htmlwidgets", "DT",
                   "base64enc", "msigdbr", "circlize"))
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Mm.eg.db", "fgsea",
                       "SummarizedExperiment", "S4Vectors", "EnhancedVolcano", "ComplexHeatmap"))
```

VISTA is the host laboratory's transcriptomics package; see its documentation at <https://cparsania.github.io/VISTA/> for installation.
