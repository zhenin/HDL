load("Figure1.rda")

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


p1 <- ggplot(Figure1a %>% filter(rg_threshold>0), aes(x = rg_threshold, y = AUC, group = interaction(h11, Method),
                                                         color = Method, shape = h11)) +
  geom_line(linetype = "dashed", size = 0.6) +
  geom_point(
    aes(shape = factor(h11), fill = Method),
    #color = "black",    # outline color
    size = 2.5
    # outline thickness
    # alpha = 0.8
  ) +
  #facet_wrap(~ h2, labeller = labeller(h2 = function(x) paste("heritability = ",x)), scales = "free_x") +
  labs(title = "10% causal SNPs",
       x = "Genetic correlation threshold for defining colocalization",
       y = "AUC") +
  theme_classic() +
  theme(legend.position = "none",
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 12, color = "black"),  # Remove x-axis title
        axis.title.y = element_text(size = 12, color = "black"),  # Change Y-axis title size and color
        axis.text.x = element_text(size = 12, color = "black"),  # Change X-axis text size and color
        axis.text.y = element_text(size = 12, color = "black"),  # Change Y-axis text size and color
        axis.line = element_line(color = "black"),
        strip.text = element_text(size =14,face = "bold", color = "black")) +
  scale_shape_manual(
    name = "Heritability for trait 1 in corresponding cis-pQTL region",
    values = c("0.001"= 21,"0.01"= 22,"0.1"= 23)  # For however many levels h11 has
  ) +
  scale_fill_manual(
    values =col_map
  ) +
  scale_color_manual(values = col_map)+
  guides(#color = guide_legend(title = "Method"),
    shape = guide_legend(title = "Heritability for trait 1 in corresponding cis-pQTL region"))+
  scale_x_continuous(
    breaks = c(0, 0.1, 0.3, 0.5, 0.7, 0.9),  # Specify x-axis breaks
    labels = c("0", ".1", ".3", ".5", ".7",".9")  # Specify x-axis labels
  )+
  scale_y_continuous(
    breaks = c(0.4,0.5, 0.7, 0.9,1),  # Specify x-axis breaks
    labels = c(".4",".5",".7", ".9","1"), # Specify x-axis labels,
    limits = c(0.45,1)
  )

p1

df_base <- Figure1bcd %>%
  filter(cau == 1) %>%
  mutate(h11_fac = factor(h11, levels = c("0.001", "0.01", "0.1")),
         x0 = as.numeric(h11_fac))

p2 <- ggplot() +
  # Lines: centered on the tick (no nudge)
  geom_line(data = df_base,
            aes(x = x0, y = AUC, color = Method, group = Method),
            linewidth = 1, alpha = 0.6, linetype = "dashed") +
  # Points for all methods EXCEPT HDL-C (0) (no nudge)
  geom_point(data = subset(df_base, Method != "HDL-C (0)"),
             aes(x = x0, y = AUC, color = Method),
             shape = 17, size = 3) +
  # Points for HDL-C (0) ONLY: nudge horizontally (tune ±0.08 as needed)
  geom_point(data = subset(df_base, Method == "HDL-C (0)"),
             aes(x = x0 - 0.02, y = AUC, color = Method),  # try -0.08 or +0.08
             shape = 17, size = 3) +
  scale_color_manual(values = col_map) +
  scale_y_continuous(breaks = c(0.7, 0.8, 0.9, 1),
                     labels = c(".7", ".8", ".9", "1"),
                     limits = c(0.65, 1)) +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c(".001", ".01", ".1")) +
  labs(title = "One causal SNP",
       x = expression(Heritability~"for"~trait~"1"~"in"~corresponding~"cis-pQTL"~regions),
       y = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.line    = element_line(color = "black"),
    plot.title   = element_text(hjust = 0.5, size = 12)
  )

p2


p3 <- ggplot(Figure1bcd%>%filter(cau==3)%>% mutate(x0 = as.numeric(h11)), aes(x = x0, y = AUC, color = Method, group = Method)) +
  geom_point(size = 3, shape =17) +
  geom_line(linewidth = 1,alpha=0.6, linetype = "dashed") +
  labs(
    title = "Three causal SNPs",
    x = expression(Heritability~"for"~trait~"1"~"in"~corresponding~"cis-pQTL"~regions),
    y = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.line    = element_line(color = "black"),
    strip.text   = element_text(size =12, face = "bold", color = "black")
  ) +
  scale_color_manual(values = col_map)+
  scale_y_continuous(
    breaks = c(0.5, 0.6,0.7,0.8),  # Specify x-axis breaks
    labels = c(".5",".6", ".7",".8"), # Specify x-axis labels,
    limits = c(0.45,0.75)
  )+
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c(".001", ".01", ".1")) 

p3


p4 <- ggplot(Figure1bcd%>%filter(cau==5)%>% mutate(x0 = as.numeric(h11)), aes(x = x0, y = AUC, color = Method, group = Method)) +
  geom_point(size = 3, shape =17) +
  geom_line(linewidth = 1,alpha=0.6, linetype = "dashed") +
  labs(
    title = "Five causal SNPs",
    x = expression(Heritability~"for"~trait~"1"~"in"~corresponding~"cis-pQTL"~regions),
    y = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black"),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black"),
    axis.line    = element_line(color = "black"),
    strip.text   = element_text(size =12, face = "bold", color = "black")
  ) +
  scale_color_manual(values = col_map)+
  scale_y_continuous(
    breaks = c(0.5, 0.6,0.7,0.8),  # Specify x-axis breaks
    labels = c(".5",".6", ".7",".8"), # Specify x-axis labels,
    limits = c(0.45,0.8)
  )+
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c(".001", ".01", ".1")) 

p4


p_combined <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
p_combined
