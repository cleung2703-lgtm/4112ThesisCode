suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(msigdbr); library(ggrepel)
})

COUNTS  <- "CPD_KO_normalised_counts.xlsx"
OUT_PNG <- "vista_figures/proliferation_score.png"
OUT_CSV <- "proliferation_scores.csv"
COL_UP  <- "#D73027"; COL_DOWN <- "#2166AC"

SETS <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
          "HALLMARK_E2F_TARGETS",    "HALLMARK_G2M_CHECKPOINT",
          "HALLMARK_MITOTIC_SPINDLE","HALLMARK_MTORC1_SIGNALING")

hm <- tryCatch(msigdbr(species = "Mus musculus", collection = "H"),
               error = function(e) msigdbr(species = "Mus musculus", category = "H"))
sym_col <- intersect(c("gene_symbol", "gene_symbol_mouse"), names(hm))[1]
set_col <- intersect(c("gs_name", "gs_id"), names(hm))[1]
hm <- hm %>% filter(.data[[set_col]] %in% SETS) %>%
  transmute(set = .data[[set_col]], symbol = .data[[sym_col]]) %>% distinct()

cnt <- read_excel(COUNTS)
samples <- setdiff(names(cnt), "gene_id")
cnt$symbol <- sub(".*:", "", cnt$gene_id)
m <- as.matrix(cnt[, samples]); rownames(m) <- cnt$symbol
m <- log2(m + 1)

z <- t(scale(t(m)))
z <- z[is.finite(rowSums(z)), , drop = FALSE]

score <- hm %>% inner_join(
  as.data.frame(z) %>% tibble::rownames_to_column("symbol"), by = "symbol") %>%
  pivot_longer(all_of(samples), names_to = "sample", values_to = "z") %>%
  group_by(set, sample) %>%
  summarise(n_genes = n(), score = mean(z), .groups = "drop") %>%
  mutate(group = ifelse(grepl("KO", sample), "Knockout", "Wild type"),
         label = tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", set)))))

readr_ok <- requireNamespace("readr", quietly = TRUE)
write.csv(score, OUT_CSV, row.names = FALSE)

cat("genes matched per set:\n")
print(score %>% distinct(label, n_genes) %>% as.data.frame(), row.names = FALSE)
cat("\nmean score by genotype:\n")
print(score %>% group_by(label, group) %>% summarise(mean = round(mean(score), 3), .groups = "drop") %>%
      pivot_wider(names_from = group, values_from = mean) %>%
      mutate(difference = round(Knockout - `Wild type`, 3)) %>% as.data.frame(), row.names = FALSE)
cat("\nsamples below zero in every set:\n")
print(score %>% group_by(sample) %>% summarise(all_neg = all(score < 0), .groups = "drop") %>% as.data.frame(),
      row.names = FALSE)

p <- ggplot(score, aes(group, score, colour = group)) +
  geom_hline(yintercept = 0, colour = "grey75") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, colour = "grey30", linewidth = 0.3) +
  geom_point(size = 2.8) +
  geom_text_repel(aes(label = sub("CPD-", "", sample)), size = 2.4, colour = "grey35",
                  nudge_x = 0.3, direction = "y", hjust = 0, segment.size = 0.2,
                  segment.colour = "grey70", min.segment.length = 0.1, seed = 1) +
  scale_colour_manual(values = c(`Wild type` = COL_DOWN, Knockout = COL_UP), guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  facet_wrap(~ label, nrow = 2) +
  labs(x = NULL, y = "Mean z-score of set members (per sample)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5))

dir.create("vista_figures", showWarnings = FALSE)
ggsave(OUT_PNG, p, width = 8.4, height = 5.6, dpi = 300, bg = "white")
cat("\nwrote", OUT_PNG, "and", OUT_CSV, "\n")
