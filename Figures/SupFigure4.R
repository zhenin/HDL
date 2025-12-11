load("SupFigure4.rda")

library(dplyr)
library(ggplot2)
library(patchwork)  # for combining plots


cat_levels <- c("Druggable", "Repurposing", "Matched+", "Matched−", "New")
cat_colors <- c(
  "Druggable"   = "#20854E",
  "Repurposing" = "#0072B5",
  "Matched+"    = "#7876B1",
  "Matched−"    = "#BC3C29",
  "New"         = "#E18727"
)


## --- Panel A: Composition across Top N ---
p_stack <- ggplot(summary_topN, aes(x = factor(TopN), y = Proportion, fill = Category)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = cat_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "",
    x = "Top N",
    y = "Proportion",
    fill = "Category"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    axis.title   = element_text(size = 12, color = "black"),
    axis.text    = element_text(size = 11, color = "black"),
    plot.title   = element_text(hjust = 0.0, size = 13, face = "bold")
  )

## --- Panel B: Matched+ proportion with exact binomial CIs ---
matched_df <- summary_topN %>%
  filter(Category == "Matched+") %>%
  transmute(
    TopN,
    Count,
    Total,
    Prop = Count / Total
  ) %>%
  rowwise() %>%
  mutate(
    lwr = binom.test(Count, Total)$conf.int[1],
    upr = binom.test(Count, Total)$conf.int[2]
  ) %>%
  ungroup()

p_matched <- ggplot(matched_df, aes(x = TopN, y = Prop)) +
  geom_line(linetype = "dashed") +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50)) +
  labs(
    title = "",
    x = "Top N",
    y = "Matched+ (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.0, size = 13, face = "bold"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text  = element_text(size = 11, color = "black")
  )

## --- Combine panels vertically (A on top, B below) ---
combined_fig <- p_stack / p_matched + plot_layout(heights = c(2, 1))

# Print
combined_fig
