suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(writexl)
})

raw_dir <- "qPCR/raw"
out_dir <- "qPCR/results"
dir.create(out_dir, showWarnings = FALSE)

TECH_SD_MAX   <- 0.5
NORT_GAP_MIN  <- 5
NTC_GAP_MIN   <- 5
REF_CT_MAX    <- 30

samples <- c("WT1", "WT2", "WT3", "KO1", "KO2", "KO3")

name_fix <- c(edrnb = "Ednrb", rpb4 = "Rbp4", lgfbp2 = "Igfbp2", papp2 = "Pappa2")

title_case <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

read_run <- function(path) {
  d <- read_excel(path, sheet = "Sheet2", col_types = "text")
  if (!"Well" %in% names(d)) d$Well <- NA_character_
  d <- d %>% transmute(
    well      = Well,
    well_name = `Well Name`,
    well_type = `Well Type`,
    rep       = as.integer(Replicate),
    ct        = suppressWarnings(as.numeric(`Ct (dR)`)))

  amp <- tolower(d$well_name)
  amp <- gsub("\\b(ntc)\\b|\\b(wt|ko)\\s*[0-9]+\\b", "", amp, perl = TRUE)
  amp <- gsub("[^a-z0-9]", "", amp)
  amp <- ifelse(amp %in% names(name_fix), unname(name_fix[amp]), title_case(amp))

  d %>% mutate(
    file     = basename(path),
    run      = sub("\\.xlsx$", "", basename(path)),
    plate_row = substr(well, 1, 1),
    plate_col = as.integer(substring(well, 2)),
    amplicon = amp,
    kind     = case_when(well_type == "NTC"   ~ "NTC",
                         well_type == "No RT" ~ "noRT",
                         TRUE                 ~ "RT"),
    sample   = ifelse(well_type == "NTC", "NTC",
                      samples[((rep - 1) %% 12) %/% 2 + 1]))
}

wells <- bind_rows(lapply(list.files(raw_dir, "\\.xlsx$", full.names = TRUE), read_run))
wells$sample <- factor(wells$sample, levels = c(samples, "NTC"))

targets <- wells %>% filter(amplicon != "B2m") %>% distinct(run, amplicon)

label_check <- wells %>%
  filter(kind != "NTC") %>%
  mutate(typed = toupper(gsub("[^a-z0-9]", "", sub(".*?(wt|ko) ?([0-9]).*", "\\1\\2",
                                                     tolower(well_name), perl = TRUE)))) %>%
  filter(typed != as.character(sample)) %>%
  distinct(file, rep, well_name, layout_sample = sample, amplicon)

reps <- wells %>%
  group_by(run, amplicon, sample, kind) %>%
  summarise(n_wells = n(), n_ct = sum(!is.na(ct)),
            ct_mean = if (n_ct > 0) mean(ct, na.rm = TRUE) else NA_real_,
            ct_diff = if (n_ct == 2) abs(diff(ct)) else NA_real_,
            cts     = paste(ifelse(is.na(ct), "NoCt", sprintf("%.2f", ct)), collapse = " / "),
            cols    = paste(sort(unique(plate_col)), collapse = ","),
            .groups = "drop")

ntc  <- reps %>% filter(kind == "NTC")  %>% select(run, amplicon, ntc_ct = ct_mean, ntc_cts = cts)
nort <- reps %>% filter(kind == "noRT") %>% select(run, amplicon, sample, nort_ct = ct_mean)

qc <- reps %>% filter(kind == "RT") %>%
  left_join(nort, by = c("run", "amplicon", "sample")) %>%
  left_join(ntc,  by = c("run", "amplicon")) %>%
  mutate(
    flag_undetected = n_ct == 0,
    flag_single     = n_ct == 1,
    flag_discordant = !is.na(ct_diff) & ct_diff > 2 * TECH_SD_MAX,
    flag_gDNA       = !is.na(ct_mean) & !is.na(nort_ct) & (nort_ct - ct_mean) < NORT_GAP_MIN,
    flag_near_NTC   = !is.na(ct_mean) & !is.na(ntc_ct) & (ntc_ct - ct_mean) < NTC_GAP_MIN,
    flag_ref_late   = amplicon == "B2m" & (is.na(ct_mean) | ct_mean > REF_CT_MAX)
  )

flag_text <- function(q) {
  f <- c(undetected = "not detected", single = "1 of 2 wells", discordant = "duplicates >1 Ct apart",
         gDNA = "noRT <5 Ct behind", near_NTC = "<5 Ct from NTC", ref_late = "B2m >30 Ct")
  apply(q[paste0("flag_", names(f))], 1, function(r) paste(f[as.logical(r)], collapse = "; "))
}
qc$flags <- flag_text(qc)

REF_EXCESS_MAX <- 5

viability <- qc %>% filter(amplicon == "B2m") %>%
  transmute(run, sample, plate_col = cols, n_ct,
            ref_ct = ct_mean, ref_cts = cts, ref_ntc = ntc_ct, ref_nort = nort_ct,
            ref_ntc_gap = ntc_ct - ct_mean,
            ref_discordant = flag_discordant,
            viable = n_ct == 2 & !is.na(ct_mean) & ct_mean <= REF_CT_MAX &
                     !is.na(ntc_ct) & (ntc_ct - ct_mean) >= NTC_GAP_MIN &
                     !flag_discordant,
            why = case_when(
              n_ct == 0                               ~ "B2m not detected",
              n_ct == 1                               ~ "B2m in only 1 of 2 wells",
              flag_discordant                         ~ "B2m duplicates disagree",
              ct_mean > REF_CT_MAX                    ~ sprintf("B2m Ct %.1f (>%d)", ct_mean, REF_CT_MAX),
              is.na(ntc_ct)                           ~ "no B2m NTC on plate",
              (ntc_ct - ct_mean) < NTC_GAP_MIN        ~ sprintf("B2m only %.1f Ct before its NTC",
                                                                ntc_ct - ct_mean),
              TRUE                                    ~ "ok")) %>%
  group_by(run) %>%
  mutate(ref_best   = suppressWarnings(min(ref_ct[viable], na.rm = TRUE)),
         ref_excess = ref_ct - ref_best,
         balanced   = viable & !is.na(ref_excess) & ref_excess <= REF_EXCESS_MAX) %>%
  ungroup() %>%
  mutate(ref_best = ifelse(is.finite(ref_best), ref_best, NA_real_))

tgt <- qc %>% filter(amplicon != "B2m")
clean_ntc <- ntc %>% filter(amplicon != "B2m") %>%
  transmute(run, gene = amplicon, target_ntc_clean = is.na(ntc_ct) | ntc_ct >= 37)

ref <- viability %>% select(run, sample, ref_ct, ref_cts, viable, balanced,
                            ref_why = why, ref_ntc_gap, ref_excess)

per_sample <- tgt %>%
  select(run, gene = amplicon, sample, tgt_ct = ct_mean, tgt_cts = cts, tgt_flags = flags,
         flag_undetected, flag_single, flag_discordant, flag_near_NTC, flag_gDNA) %>%
  left_join(ref, by = c("run", "sample")) %>%
  left_join(clean_ntc, by = c("run", "gene")) %>%
  mutate(group = factor(substr(as.character(sample), 1, 2), levels = c("WT", "KO")),
         dct = tgt_ct - ref_ct,
         pass = viable & !flag_undetected & !flag_single & !flag_discordant &
                !flag_near_NTC & !is.na(dct),
         censored = balanced & flag_undetected & target_ntc_clean & !is.na(ref_ct),
         dct_cens = ifelse(censored, 40 - ref_ct, dct))

summarise_run <- function(ps, set_name) {
  ps %>% filter(!is.na(dct)) %>%
    group_by(run, gene) %>%
    group_modify(function(d, k) {
      wt <- d$dct[d$group == "WT"]; ko <- d$dct[d$group == "KO"]
      if (length(wt) < 1 || length(ko) < 1)
        return(tibble(n_wt = length(wt), n_ko = length(ko)))
      ddct <- mean(ko) - mean(wt)
      p <- if (length(wt) >= 2 && length(ko) >= 2) t.test(ko, wt)$p.value else NA_real_
      tibble(n_wt = length(wt), n_ko = length(ko),
             dct_wt = mean(wt), dct_ko = mean(ko),
             sd_wt = if (length(wt) > 1) sd(wt) else NA, sd_ko = if (length(ko) > 1) sd(ko) else NA,
             log2fc_qpcr = -ddct, fold_change = 2^-ddct, welch_p = p,
             wt_samples = paste(d$sample[d$group == "WT"], collapse = ","),
             ko_samples = paste(d$sample[d$group == "KO"], collapse = ","))
    }) %>% ungroup() %>% mutate(sample_set = set_name, .before = 1)
}

summary_all  <- summarise_run(per_sample, "all samples with a Ct")
summary_pass <- summarise_run(filter(per_sample, pass), "QC-passed samples only")
summary_cens <- summarise_run(per_sample %>% filter(pass | censored) %>% mutate(dct = dct_cens),
                              "QC-passed + non-detects at Ct 40") %>%
  filter(paste(run, gene) %in% with(filter(per_sample, censored), paste(run, gene)))
summary_bal  <- summarise_run(filter(per_sample, pass & balanced), "QC-passed, balanced input only")

de <- read_excel("MEF_CPD_KO_de_genes.xlsx") %>%
  mutate(gene = sub(".*:", "", gene_id)) %>%
  select(gene, rnaseq_log2fc = log2FoldChange, rnaseq_padj = padj, baseMean)

summary_tbl <- bind_rows(summary_all, summary_pass, summary_bal, summary_cens) %>%
  left_join(de, by = "gene") %>%
  mutate(same_direction = sign(log2fc_qpcr) == sign(rnaseq_log2fc),
         verdict = case_when(
           is.na(log2fc_qpcr) ~ "not assessable",
           n_wt < 2 | n_ko < 2 ~ "too few samples",
           same_direction & !is.na(welch_p) & welch_p < 0.05 ~ "validated",
           same_direction ~ "same direction, not significant",
           TRUE ~ "opposite direction"))

cat("\nRuns read:", length(unique(wells$run)),
    " targets:", paste(sort(unique(targets$amplicon)), collapse = ", "), "\n\n")
cat("Wells that held template (viable), by run:\n")
viability %>%
  group_by(run) %>%
  summarise(held = paste(sample[viable], collapse = " "),
            depleted = paste(sample[viable & !balanced], collapse = " "),
            n_wt = sum(balanced & grepl("^WT", sample)),
            n_ko = sum(balanced & grepl("^KO", sample)), .groups = "drop") %>%
  as.data.frame() %>% print(right = FALSE)
cat("\nGenes with a usable KO vs WT comparison:\n")
summary_tbl %>% filter(sample_set == "QC-passed samples only", !is.na(log2fc_qpcr)) %>%
  select(gene, n_wt, n_ko, log2fc_qpcr, fold_change, welch_p, rnaseq_log2fc, verdict) %>%
  as.data.frame() %>% print(digits = 3, right = FALSE)

write_xlsx(list(
  summary      = summary_tbl,
  per_sample   = per_sample,
  viability    = viability,
  well_qc      = qc %>% select(run, amplicon, sample, cols, cts, ct_mean, ct_diff,
                               nort_ct, ntc_cts, flags),
  label_check  = label_check,
  raw_wells    = wells),
  file.path(out_dir, "qPCR_results.xlsx"))

COL_UP <- "#D73027"; COL_DOWN <- "#2166AC"; COL_NS <- "#BDBDBD"
pal <- c(WT = COL_DOWN, KO = COL_UP)

run_label <- c(
  "20261001cl.cpd.qpcr (1)"         = "1   Cpd",
  "20261002cl.cpd.pcr"              = "2   Cpd",
  "20261006cl.cpd.pcr"              = "3   Cpd",
  "20261008cl.lgfbp2.qpcr"          = "4   Igfbp2",
  "20261008cl.mylk.qpcr"            = "5   Mylk",
  "20261009cl.hmga2.qpcr"           = "6   Hmga2",
  "20261011cl.fgf18.qpcr"           = "7   Fgf18",
  "20261013cl.msx1.qpcr"            = "8   Msx1",
  "20261017cl.ednrb.rasa1.qpcr"     = "9   Ednrb + Rasa1",
  "202161017cl.isl2.rbp4.qpcr"      = "10  Isl2 + Rbp4",
  "202611089cl.ablim1.shroom3.qpcr" = "11  Ablim1 + Shroom3",
  "20261108cl.il1rl1.epas1.qpcr"    = "12  Il1rl1 + Epas1",
  "20261108cl.smpd3.hoxa7.qpcr"     = "13  Smpd3 + Hoxa7",
  "20261108cl.papp2.qpcr"           = "14  Pappa2")
short_run <- function(x) ifelse(x %in% names(run_label), unname(run_label[x]), x)
run_levels <- rev(unname(run_label))

usable <- summary_tbl %>% filter(sample_set == "QC-passed samples only",
                                 !is.na(log2fc_qpcr)) %>% pull(gene)

ps_plot <- per_sample %>% filter(pass, gene %in% usable) %>%
  group_by(gene) %>%
  mutate(rel = 2^-(dct - mean(dct[group == "WT"], na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(gene = factor(gene, levels = usable[order(-summary_tbl$rnaseq_log2fc[
           match(usable, summary_tbl$gene)])]))

p1 <- ggplot(ps_plot, aes(group, log2(rel), colour = group)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, colour = "grey30", linewidth = 0.3) +
  geom_point(size = 2.8) +
  geom_text(aes(label = sample), size = 2.5, nudge_x = 0.3, colour = "grey35") +
  geom_hline(data = distinct(select(ps_plot, gene)) %>%
               left_join(distinct(select(summary_tbl, gene, rnaseq_log2fc)), by = "gene"),
             aes(yintercept = rnaseq_log2fc), linetype = 2, colour = "grey35", linewidth = 0.4) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  facet_wrap(~ gene, nrow = 2) +
  labs(x = NULL, y = "log2 expression relative to the wild-type mean (B2m-normalised)") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5))
ggsave(file.path(out_dir, "fig1_per_sample_expression.png"), p1, width = 8.5, height = 5.4, dpi = 200)

cmp <- summary_tbl %>%
  filter(sample_set %in% c("QC-passed samples only", "QC-passed, balanced input only"),
         !is.na(log2fc_qpcr))
if (nrow(cmp)) {
  p2 <- ggplot(cmp, aes(rnaseq_log2fc, log2fc_qpcr, colour = sample_set)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "grey95") +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey95") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(size = 2.8) +
    geom_text(aes(label = gene), size = 2.9, vjust = -1, show.legend = FALSE) +
    scale_colour_manual(values = c("QC-passed samples only" = COL_DOWN,
                                   "QC-passed, balanced input only" = "#7A9A3B")) +
    coord_cartesian(xlim = c(-2, 3.4), ylim = c(-2, 4.2)) +
    labs(x = "RNA-seq log2 fold change (DESeq2, knockout versus wild type)",
         y = "qPCR log2 fold change (2^-ddCt)", colour = NULL,
         caption = "shaded quadrants indicate disagreement in direction; dashed line is exact agreement") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", legend.direction = "vertical",
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
          plot.caption = element_text(size = 7.5, colour = "grey40"))
  ggsave(file.path(out_dir, "fig2_qpcr_vs_rnaseq.png"), p2, width = 6, height = 6.4, dpi = 200)
}

b2m <- viability %>% mutate(run_s = factor(short_run(run), levels = rev(run_levels)))
p3 <- ggplot(b2m, aes(sample, ref_ct)) +
  geom_hline(aes(yintercept = ref_ntc), colour = COL_UP, linetype = 2, linewidth = 0.4) +
  geom_point(aes(colour = viable), size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = COL_DOWN, `FALSE` = COL_NS),
                      labels = c(`TRUE` = "held template", `FALSE` = "failed"), name = NULL) +
  facet_wrap(~ run_s, ncol = 4) +
  labs(x = NULL, y = "B2m Ct (dashed line = that plate's B2m NTC)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "fig3_b2m_reference_ct.png"), p3, width = 10, height = 6, dpi = 200)

p4 <- viability %>%
  mutate(run_s  = factor(short_run(run), levels = run_levels),
         state  = case_when(!viable   ~ "no usable template",
                            !balanced ~ "template, but depleted",
                            TRUE      ~ "held template")) %>%
  ggplot(aes(sample, run_s, fill = state)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = ifelse(is.na(plate_col) | plate_col == "" | plate_col == "NA",
                              "", paste0("col ", plate_col))),
            size = 2.3, colour = "white") +
  scale_fill_manual(values = c("held template"          = COL_DOWN,
                               "template, but depleted" = "#92B4D4",
                               "no usable template"     = COL_NS), name = NULL) +
  labs(x = NULL, y = NULL,
       title = "Wells whose B2m amplified clear of that plate's own no-template control") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        plot.title = element_text(size = 9, face = "plain"))
ggsave(file.path(out_dir, "fig4_well_viability.png"), p4, width = 8.5, height = 4.2, dpi = 200)

cat("\nWrote", file.path(out_dir, "qPCR_results.xlsx"), "and 4 figures.\n")
