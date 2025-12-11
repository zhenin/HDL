load("SupFigure3.rda")


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

col_map <- c("COLOC"="#20854E", "SuSiE"="#0072B5", "SharePro"="#7876B1", "HDL-C"= "#BC3C29")
SupFigure3$method = factor(SupFigure3$method, levels = c("HDL-C", "COLOC","SuSiE","SharePro"))

p_overall <- ggplot(SupFigure3, aes(x = method, y = elapsed, fill = method)) +
  geom_boxplot(width = 0.65, outlier.shape = 46) +
  scale_y_log10(name = "Execution time (seconds, log scale)",
                breaks = c(0.01, 0.1, 1, 10),
                labels = c("0.01", "0.1", "1", "10")) +
  coord_cartesian(ylim = c(NA, 20)) +  
  scale_fill_manual(values = col_map) +
  labs(x = NULL, title = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
p_overall
