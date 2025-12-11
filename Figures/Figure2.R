load("Figure2.rda")

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

# NEJM colors
nejm.colors <- function (n = 8, alpha = 1){
  pal <- c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
           "#6F99AD", "#FFDC91", "#EE4C97")
  acode <- substr(rgb(0, 0, 0, alpha = alpha), 8, 9)
  paste0(pal, acode)[1:n]
}

col_map <- c("COLOC"="#20854E", "SuSiE"="#0072B5", "SharePro"="#7876B1", "HDL-C 0"="#BC3C29","HDL-C 0.5" = "#BC3C29")
linemap <- c(
  "HDL-C 0"   = "solid",
  "HDL-C 0.5" = "dashed",
  "COLOC"     = "dashed",
  "SuSiE"     = "dashed",
  "SharePro"  = "dashed"
)

p2real = ggplot(
  Figure2,
  aes(
    x = topN,
    y = rediscoveryRate,
    color = MethodOnly,      # Different color per method
    #shape = MethodOnly,      # Different point shape per method
    linetype = MethodOnly    # Different line type per method
  )
) +
  geom_point(size = 4) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~Direction, ncol = 2) +      # Stack panels: one for M->F, one for F->M
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "",
    x = "Top N in training set",
    y = "Rediscovery rate in test set",
    color = "Method",
    shape = "Method",
    linetype = "Method"
  ) +
  theme_classic() +
  # Use your NEJM color palette for the three method labels
  scale_color_manual(values = col_map) +
  scale_linetype_manual(values = linemap) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.line    = element_line(color = "black"),
    strip.text   = element_text(size = 12, color = "black"),
    text         = element_text( size = 12)
  )

p2real

