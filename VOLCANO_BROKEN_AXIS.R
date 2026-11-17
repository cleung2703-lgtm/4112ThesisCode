suppressPackageStartupMessages({
  library(VISTA); library(ggplot2); library(patchwork)
})
COLORS <- c(Up = "#D73027", Down = "#2166AC", Other = "grey80")
PADJ <- 0.05; LFC <- 1
vista <- readRDS("vista_object_shared.rds")
comp  <- names(comparisons(vista))[1]

p0 <- suppressWarnings(get_volcano_plot(
  vista, sample_comparison = comp, colors = COLORS,
  display_id = "SYMBOL", label_genes = NULL, point_size = 1.2))
d <- p0$data
d$y <- -log10(d$padj)
d <- d[is.finite(d$y) & is.finite(d$xvals), ]

BRK_LO <- 15; BRK_HI <- 65
stopifnot(sum(d$y > BRK_LO & d$y < BRK_HI) == 0)
cat(sprintf("genes plotted %d | in the break interval %d\n",
            nrow(d), sum(d$y > BRK_LO & d$y < BRK_HI)))

top_genes <- d$lab[d$y >= BRK_HI]
cat("upper panel:", paste(top_genes, collapse = ", "), "\n")

lower <- d[d$y < BRK_HI & d$regulation != "Other", ]
lower <- lower[order(lower$padj), ]
lab_lo <- c(head(lower$lab[lower$regulation == "Up"],   11),
            head(lower$lab[lower$regulation == "Down"], 11))

XR <- range(d$xvals)
XL <- c(XR[1] - diff(XR) * 0.055, XR[2] + diff(XR) * 0.01)

slash <- function(yc, h) {
  w <- diff(XL) * 0.007
  lapply(c(XL[1] + diff(XL) * 0.013, XL[1] + diff(XL) * 0.032), function(xc)
    annotate("segment", x = xc - w, xend = xc + w, y = yc - h, yend = yc + h,
             colour = "grey20", linewidth = 0.7))
}

vista_panel <- function(ylim, labels) suppressWarnings(get_volcano_plot(
  vista, sample_comparison = comp, colors = COLORS, display_id = "SYMBOL",
  label_genes = labels, label_size = 3.1, point_size = 1.2,
  pval_cutoff = PADJ, log2fc_cutoff = LFC,
  ylim = ylim, xlim = XL,
  legendPosition = "right", legendLabSize = 10, legendIconSize = 3,
  drawConnectors = TRUE, widthConnectors = 0.3, colConnectors = "grey55",
  min.segment.length = 0.1, max.overlaps = 20,
  gridlines.major = FALSE, gridlines.minor = FALSE,
  border = "full", borderWidth = 0.5, borderColour = "grey30",
  axisLabSize = 11))

p_top <- vista_panel(c(64.4, 71.2), top_genes) +
  labs(title = NULL, subtitle = NULL, caption = NULL, x = NULL, y = NULL) +
  scale_y_continuous(limits = c(64.4, 71.2), breaks = c(65, 68, 71)) +
  slash(64.95, 0.45) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none", plot.margin = margin(4, 5, 0, 5))

p_bot <- vista_panel(c(0, 15.8), lab_lo) +
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  scale_y_continuous(limits = c(0, 15.8), breaks = seq(0, 15, 5)) +
  slash(15.45, 0.25) +
  theme(plot.margin = margin(0, 5, 4, 5))

fig <- (p_top / p_bot) +
  plot_layout(heights = c(1, 5.0), guides = "collect") &
  theme(legend.position = "right")

dir.create("vista_figures", showWarnings = FALSE)
ggsave("vista_figures/volcano_broken_axis.png", fig,
       width = 9, height = 7.6, dpi = 300, bg = "white")
cat("wrote vista_figures/volcano_broken_axis.png\n")
