suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

IN_CSV  <- "gsea_sensitivity_padj_comparison.csv"
OUT_PNG <- "vista_figures/gsea_sensitivity_heatmap.png"
ALPHA   <- 0.05

d <- read.csv(IN_CSV, stringsAsFactors = FALSE, check.names = FALSE)

tidy <- d %>%
  pivot_longer(c(log2FC_full, stat_full), names_to = "metric", values_to = "padj") %>%
  mutate(
    set    = gsub("_", " ", sub("^HALLMARK_", "", pathway)),
    metric = factor(metric,
                    levels = c("log2FC_full", "stat_full"),
                    labels = c("Ranked by\nlog2 fold change", "Ranked by\nWald statistic")),
    sig      = !is.na(padj) & padj < ALPHA,
    neglog10 = ifelse(is.na(padj), 0, -log10(pmax(padj, 1e-40))))

keep <- tidy %>% group_by(set) %>% filter(any(sig)) %>% ungroup()
dropped <- setdiff(unique(tidy$set), unique(keep$set))

ord <- keep %>% group_by(set) %>% summarise(m = max(neglog10), .groups = "drop") %>% arrange(m)
keep$set <- factor(keep$set, levels = ord$set)

n_both <- keep %>% group_by(set) %>% summarise(n = sum(sig), .groups = "drop") %>% count(n)
cat(sprintf("%d of %d Hallmark sets significant under at least one ranking\n",
            length(unique(keep$set)), length(unique(tidy$set))))
print(n_both)
cat("significant under neither:\n  ", paste(dropped, collapse = ", "), "\n")

p <- ggplot(keep, aes(metric, set, fill = neglog10)) +
  geom_tile(colour = "white", linewidth = 0.9) +
  geom_text(aes(label = ifelse(sig, "*", ""),
                colour = neglog10 > 0.45 * max(keep$neglog10)),
            vjust = 0.72, size = 6, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey25")) +
  scale_fill_gradient(low = "#F7F3F2", high = "#D73027", name = expression(-log[10] ~ p[adj])) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       caption = sprintf("* marks an adjusted p-value below %.2f", ALPHA)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 8.5),
        axis.text.x = element_text(size = 9, face = "bold"),
        legend.position = "right",
        legend.key.height = unit(1.1, "cm"),
        plot.caption = element_text(hjust = 0, size = 8, colour = "grey35"))

dir.create("vista_figures", showWarnings = FALSE)
ggsave(OUT_PNG, p, width = 7.2, height = 6.2, dpi = 300, bg = "white")
cat("wrote", OUT_PNG, "\n")
