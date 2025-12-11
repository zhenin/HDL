load("SFigure2.rda")

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

col_map <- c("COLOC"="#20854E", "SuSiE"="#0072B5", "SharePro"="#7876B1", "HDL-C (0)"="#BC3C29","HDL-C (0.5)" = "#E18727")

tproc = ggplot(SFigure2, aes(x = Method, y = TPR, fill = Method)) +
  geom_col(width = 0.8) +
  facet_grid(cau_lab ~ h11.0 + target, labeller = labeller(
    h11_lab = function(x) paste0(x),
    cau_lab = function(x) paste0(x))) +
  coord_cartesian(ylim = c(0,1)) +
  scale_fill_manual(values = col_map) +
  labs(x = "", y = "True positive rate") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 0.9))

tproc
