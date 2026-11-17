suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(ggplot2); library(writexl)
})

ALPHA <- 0.05
POWER <- 0.80
ATTEN <- 0.6

out_dir <- "qPCR/results"

ps <- read_excel(file.path(out_dir, "qPCR_results.xlsx"), sheet = "per_sample") %>%
  mutate(across(c(pass, viable, balanced), as.logical)) %>%
  filter(pass, !is.na(dct))

de <- read.csv("vista_tables/de_results_all_genes.csv", stringsAsFactors = FALSE) %>%
  select(gene = SYMBOL, rnaseq_log2fc = log2fc, rnaseq_padj = padj)

var_parts <- ps %>%
  group_by(run, gene, group) %>%
  summarise(n = n(), v = if (n() > 1) var(dct) else NA_real_, .groups = "drop") %>%
  filter(!is.na(v))

df_total  <- sum(var_parts$n - 1)
sd_pooled <- sqrt(sum((var_parts$n - 1) * var_parts$v) / df_total)

cat(sprintf("\nPooled within-group SD of dCt: %.3f cycles (%d df from %d group estimates)\n",
            sd_pooled, df_total, nrow(var_parts)))
cat(sprintf("Range of the individual SDs: %.2f to %.2f cycles\n",
            min(sqrt(var_parts$v)), max(sqrt(var_parts$v))))

keep <- ps %>% distinct(run, sample)
dups <- read_excel(file.path(out_dir, "qPCR_results.xlsx"), sheet = "well_qc") %>%
  mutate(ct_diff = as.numeric(ct_diff)) %>%
  filter(!is.na(ct_diff)) %>%
  semi_join(keep, by = c("run", "sample"))

sd_tech <- sqrt(mean(dups$ct_diff^2) / 2)
sd_sample <- sqrt(max(sd_pooled^2 - sd_tech^2, 0))

cat(sprintf("
Technical SD from %d duplicate pairs: %.3f cycles, which carries to dCt unchanged
",
            nrow(dups), sd_tech))
cat(sprintf("Sample-to-sample SD after removing it: %.3f cycles (%.0f%% of the total variance)
",
            sd_sample, 100 * sd_sample^2 / sd_pooled^2))

power_at <- function(delta, n1, n2, sd, alpha = ALPHA) {
  if (n1 < 2 || n2 < 2) return(NA_real_)
  df  <- n1 + n2 - 2
  ncp <- delta / (sd * sqrt(1 / n1 + 1 / n2))
  tc  <- qt(1 - alpha / 2, df)
  pt(-tc, df, ncp) + (1 - pt(tc, df, ncp))
}

mdd <- function(n1, n2, sd, power = POWER) {
  if (n1 < 2 || n2 < 2) return(NA_real_)
  f <- function(d) power_at(d, n1, n2, sd) - power
  if (f(20) < 0) return(NA_real_)
  uniroot(f, c(1e-6, 20), tol = 1e-6)$root
}

grid <- expand.grid(n_wt = 1:3, n_ko = 1:3) %>%
  rowwise() %>%
  mutate(mdd_dct  = mdd(n_wt, n_ko, sd_pooled),
         mdd_fold = 2^mdd_dct) %>%
  ungroup() %>%
  arrange(desc(n_wt + n_ko))

cat("\nSmallest detectable change at 80% power, by group size:\n")
grid %>% filter(!is.na(mdd_dct)) %>%
  transmute(n_wt, n_ko,
            `min detectable log2FC` = round(mdd_dct, 2),
            `as fold change` = round(mdd_fold, 1)) %>%
  as.data.frame() %>% print(right = FALSE)

achieved <- ps %>%
  group_by(gene) %>%
  summarise(n_wt = sum(group == "WT"), n_ko = sum(group == "KO"), .groups = "drop") %>%
  left_join(de, by = "gene") %>%
  rowwise() %>%
  mutate(
    mdd_dct   = mdd(n_wt, n_ko, sd_pooled),
    expected  = abs(rnaseq_log2fc) * ATTEN,
    power_seq = power_at(abs(rnaseq_log2fc), n_wt, n_ko, sd_pooled),
    power_exp = power_at(expected, n_wt, n_ko, sd_pooled)) %>%
  ungroup() %>%
  mutate(reach = case_when(
    is.na(mdd_dct)            ~ "not testable (one group had a single sample)",
    expected >= mdd_dct       ~ "within reach",
    TRUE                      ~ "below the detection limit")) %>%
  arrange(desc(abs(rnaseq_log2fc)))

cat("\nPer gene, against the change the sequencing predicted:\n")
achieved %>%
  transmute(gene, n_wt, n_ko,
            seq_log2FC = round(rnaseq_log2fc, 2),
            expected_log2FC = round(expected, 2),
            min_detectable = round(mdd_dct, 2),
            power = round(power_exp, 2), reach) %>%
  as.data.frame() %>% print(right = FALSE)

n_test <- sum(!is.na(achieved$mdd_dct))
n_in   <- sum(achieved$reach == "within reach")
cat(sprintf("\n%d of the %d genes with a usable comparison had a predicted change\n", n_in, n_test))
cat(sprintf("large enough to be detectable; median power across them was %.2f.\n",
            median(achieved$power_exp, na.rm = TRUE)))

med_exp <- median(achieved$expected, na.rm = TRUE)
ns <- 2:400
need <- sapply(ns, function(n) power_at(med_exp, n, n, sd_pooled))
n_need <- ns[which(need >= POWER)[1]]
cat(sprintf("\nFor the median expected change (log2FC %.2f, %.1f-fold), reaching 80%% power\n",
            med_exp, 2^med_exp))
cat(sprintf("would require n = %d per group at the observed variability.\n", n_need))

COL_UP <- "#D73027"; COL_DOWN <- "#2166AC"
curve_df <- do.call(rbind, lapply(c("2 v 3" = 2, "3 v 3" = 3), function(n) {
  d <- seq(0, 4, length.out = 400)
  data.frame(delta = d, power = sapply(d, power_at, n1 = n, n2 = 3, sd = sd_pooled),
             design = paste0("n = ", n, " vs 3"))
}))

pts <- achieved %>% filter(!is.na(mdd_dct)) %>%
  mutate(design = paste0("n = ", pmin(n_wt, n_ko), " vs ", pmax(n_wt, n_ko)))

p <- ggplot(curve_df, aes(delta, power, linetype = design)) +
  geom_hline(yintercept = POWER, colour = COL_UP, linewidth = 0.4) +
  annotate("text", x = 3.9, y = POWER + 0.04, label = "80% power",
           hjust = 1, size = 2.9, colour = COL_UP) +
  geom_line(linewidth = 0.6, colour = "grey25") +
  geom_rug(data = achieved, aes(x = expected), inherit.aes = FALSE,
           sides = "b", colour = "grey45", length = unit(0.035, "npc")) +
  annotate("segment", x = min(achieved$expected, na.rm = TRUE),
           xend = max(achieved$expected, na.rm = TRUE),
           y = 0.145, yend = 0.145, colour = "grey45", linewidth = 0.4) +
  annotate("segment", x = 0.42, xend = 0.42, y = 0.155, yend = 0.56,
           colour = "grey60", linewidth = 0.3) +
  annotate("text", x = 0.30, y = 0.58, hjust = 0, vjust = 0, size = 2.8, colour = "grey30",
           label = "every change the sequencing
predicted for the 17 assayed genes
falls in this range") +
  ggrepel::geom_text_repel(data = pts, aes(x = expected, y = power_exp, label = gene),
                           inherit.aes = FALSE, size = 2.8, colour = "grey25",
                           nudge_y = -0.05, direction = "x", min.segment.length = 0, seed = 1) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
  scale_linetype_manual(values = c("n = 2 vs 3" = 2, "n = 3 vs 3" = 1), name = NULL) +
  labs(x = "expected difference in dCt (log2 fold change after attenuation)",
       y = "probability of reaching p < 0.05") +
  theme_bw(base_size = 10) +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5))
ggsave(file.path(out_dir, "fig5_detectability.png"), p, width = 6.2, height = 4.4, dpi = 300)

write_xlsx(list(per_gene = achieved, group_sizes = grid, variance_parts = var_parts),
           file.path(out_dir, "qPCR_power.xlsx"))

cat("\nWrote", file.path(out_dir, "qPCR_power.xlsx"), "and fig5_detectability.png\n")
