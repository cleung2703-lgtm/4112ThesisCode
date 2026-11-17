suppressPackageStartupMessages({
  library(msigdbr); library(fgsea); library(dplyr); library(ggplot2); library(writexl)
})

set.seed(42)
out_dir <- "vista_tables"

predicted <- tribble(
  ~gs_name,                                        ~sign, ~tier,     ~basis,
  "MULLIGAN_NTF3_SIGNALING_VIA_INSR_AND_IGF1R_UP",   -1,  "primary", "induced by signalling through INSR and IGF1R",
  "MULLIGAN_NTF3_SIGNALING_VIA_INSR_AND_IGF1R_DN",   +1,  "primary", "repressed by signalling through INSR and IGF1R",
  "PACHER_TARGETS_OF_IGF1_AND_IGF2_UP",              -1,  "primary", "induced by IGF1 and IGF2",
  "CREIGHTON_AKT1_SIGNALING_VIA_MTOR_UP",            -1,  "primary", "induced by AKT1 signalling through mTOR",
  "CREIGHTON_AKT1_SIGNALING_VIA_MTOR_DN",            +1,  "primary", "repressed by AKT1 signalling through mTOR",
  "PARENT_MTOR_SIGNALING_UP",                        -1,  "primary", "induced by mTOR signalling",
  "PARENT_MTOR_SIGNALING_DN",                        +1,  "primary", "repressed by mTOR signalling",
  "XU_AKT1_TARGETS_6HR",                            -1,  "primary", "induced by AKT1 at 6 h",
  "XU_AKT1_TARGETS_48HR",                           -1,  "primary", "induced by AKT1 at 48 h",
  "REICHERT_G1S_REGULATORS_AS_PI3K_TARGETS",         -1,  "primary", "G1/S regulators downstream of PI3K",
  "BILANGES_RAPAMYCIN_SENSITIVE_GENES",             -1,  "primary", "expression requires mTOR activity",
  "SCHLOSSER_SERUM_RESPONSE_UP",                    -1,  "secondary", "induced by serum growth factors",
  "SCHLOSSER_SERUM_RESPONSE_DN",                    +1,  "secondary", "repressed by serum growth factors",
  "GRAESSMANN_APOPTOSIS_BY_SERUM_DEPRIVATION_UP",   +1,  "secondary", "induced by serum withdrawal",
  "GRAESSMANN_APOPTOSIS_BY_SERUM_DEPRIVATION_DN",   -1,  "secondary", "repressed by serum withdrawal")

de <- read.csv("vista_tables/de_results_all_genes.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(stat), !is.na(SYMBOL), SYMBOL != "") %>%
  group_by(SYMBOL) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()

ranks <- setNames(de$stat, de$SYMBOL)
ranks <- sort(ranks, decreasing = TRUE)
cat(sprintf("\nRanked list: %d genes, statistic from %.1f to %.1f\n",
            length(ranks), min(ranks), max(ranks)))

cgp <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CGP")
sets <- split(cgp$gene_symbol, cgp$gs_name)
sets <- lapply(sets, function(g) unique(intersect(g, names(ranks))))
sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
cat(sprintf("Testing %d curated perturbation sets\n", length(sets)))

res <- fgsea(pathways = sets, stats = ranks, minSize = 10, maxSize = 500,
             nPermSimple = 10000) %>%
  as_tibble() %>%
  mutate(rank_nes = rank(NES))

missing <- setdiff(predicted$gs_name, res$pathway)
if (length(missing)) cat("\nNot available after filtering:", paste(missing, collapse = ", "), "\n")

tested <- predicted %>%
  inner_join(select(res, gs_name = pathway, size, NES, pval, padj, rank_nes),
             by = "gs_name") %>%
  mutate(as_predicted = sign(NES) == sign,
         percentile = 100 * rank_nes / nrow(res),
         competitive = ifelse(sign < 0, percentile < 50, percentile > 50)) %>%
  arrange(tier, desc(as_predicted), NES)

cat("\n--- Pre-specified sets, predicted direction against observed ---\n")
tested %>%
  transmute(set = substr(gs_name, 1, 44), tier, n = size,
            predicted = ifelse(sign < 0, "down", "up"),
            NES = round(NES, 2), padj = signif(padj, 2),
            matched = ifelse(as_predicted, "yes", "no"),
            pctile = round(percentile)) %>%
  as.data.frame() %>% print(right = FALSE)

global_med <- median(res$NES, na.rm = TRUE)
cat(sprintf("
Median NES across all %d sets: %.2f, so the collection is globally depleted.
",
            nrow(res), global_med))

naive <- filter(tested, tier == "primary")
k <- sum(naive$as_predicted)
cat(sprintf("Naive sign test against zero: %d of %d primary sets matched, p = %.3f.
",
            k, nrow(naive), binom.test(k, nrow(naive), 0.5, alternative = "greater")$p.value))
cat("That test is reported for completeness but the global skew makes it uninformative.
")

down <- filter(tested, tier == "primary", sign < 0)
up   <- filter(tested, tier == "primary", sign > 0)
bg   <- filter(res, !pathway %in% predicted$gs_name)

w1 <- wilcox.test(down$NES, bg$NES, alternative = "less")
cat(sprintf("
Sets predicted down (n = %d) versus the %d background sets: p = %.4f
",
            nrow(down), nrow(bg), w1$p.value))
cat(sprintf("  median NES %.2f against %.2f, median percentile %.0f
",
            median(down$NES), median(bg$NES), median(down$percentile)))

comp <- filter(tested, tier == "primary")
kc <- sum(comp$competitive)
cat(sprintf("
Competitive sign test against the collection median: %d of %d primary sets
",
            kc, nrow(comp)))
cat(sprintf("fell on the predicted side, p = %.3f.
",
            binom.test(kc, nrow(comp), 0.5, alternative = "greater")$p.value))

if (nrow(up) >= 2) {
  w2 <- wilcox.test(down$NES, up$NES, alternative = "less")
  cat(sprintf("
Sets predicted down versus the %d predicted up: p = %.3f, medians %.2f and %.2f.
",
              nrow(up), w2$p.value, median(down$NES), median(up$NES)))
}

sig <- tested %>% filter(padj < 0.05) %>% arrange(NES)
cat(sprintf("
%d pre-specified sets reached padj < 0.05:
", nrow(sig)))
sig %>% transmute(gs_name, NES = round(NES, 2), padj = signif(padj, 2),
                  matched = ifelse(as_predicted, "as predicted", "against prediction"),
                  basis) %>%
  as.data.frame() %>% print(right = FALSE)

COL_UP <- "#D73027"; COL_DOWN <- "#2166AC"
plot_df <- tested %>%
  mutate(label = gsub("_", " ", tolower(gs_name)),
         label = paste0(label, "  (n = ", size, ")"),
         label = factor(label, levels = label[order(NES)]),
         fill = ifelse(competitive, "as predicted", "against prediction"))

p <- ggplot(plot_df, aes(NES, label, fill = fill)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_vline(xintercept = global_med, colour = "grey45", linetype = 2, linewidth = 0.4) +
  geom_col(width = 0.68) +
  geom_text(aes(label = ifelse(padj < 0.05, "*", ""),
                hjust = ifelse(NES < 0, 1.6, -0.6)), size = 4, colour = "grey20") +
  annotate("text", x = global_med, y = 0.3, label = "collection median",
           angle = 90, vjust = -0.5, hjust = 0, size = 2.6, colour = "grey45") +
  scale_fill_manual(values = c("as predicted" = COL_DOWN, "against prediction" = COL_UP),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  facet_grid(tier ~ ., scales = "free_y", space = "free_y") +
  labs(x = "normalised enrichment score in Cpd KO versus WT", y = NULL,
       caption = paste("* padj < 0.05. The direction each set should move was fixed before testing:",
                       "a set induced by the pathway
should be depleted if the pathway is attenuated.",
                       "Because most sets in the collection are depleted, the dashed
line rather than zero is the",
                       "reference against which a set counts as moved.")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        plot.caption = element_text(size = 7, colour = "grey35", hjust = 0),
        plot.caption.position = "plot",
        strip.background = element_rect(fill = "grey92", colour = "grey30"))
ggsave("external_validation_igf1r.png", p, width = 7.6, height = 5.4, dpi = 300)

write_xlsx(list(prespecified = tested,
                all_cgp = arrange(res, NES) %>%
                  mutate(leadingEdge = vapply(leadingEdge, paste, "", collapse = " "))),
           file.path(out_dir, "external_validation_cgp.xlsx"))
cat("\nWrote external_validation_igf1r.png and", file.path(out_dir, "external_validation_cgp.xlsx"), "\n")
