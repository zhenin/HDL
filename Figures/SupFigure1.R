load("SupFigure1.rda")

library(ggplot2)
library(grid)
library(dplyr)
library(cowplot)
#Draw a plot for checking the distribution of herita
# Define the custom color palette function
npgc <- function(n = 15, alpha = 1) {
  pal <- c("#000000","#9c3386", "#b97dac", "#d3bad6","#ffb82a","#eeaf65","#007bbf","#5db8e1", "#009675","#00c69b","#c65848", "#f55d3b","#f38750", "#8e816d","#dadcdb")
  acode <- substr(rgb(0, 0, 0, alpha = alpha), 8, 9)
  return(paste0(pal, acode)[1:n])
}

# Kolmogorov-Smirnov test
set.seed(2025)
lab = sample(1:nrow(Topsnph2), 300)
ks_test <- ks.test(Topsnph2$h2, Topsnph2$h2[lab])
p_value <- ks_test$p.value
p_value_text <- expression(italic(P) * "-value of K-S test: 0.2404")

# Common theme settings for axis
common_theme <- theme_classic() +
  theme(
    axis.line = element_line(color = "black"),    # Set axis lines to black
    axis.ticks = element_line(color = "black"),   # Set axis ticks to black
    axis.text = element_text(color = "black", size =12),    # Set axis text to black
    axis.title = element_text(color = "black", size =14)    # Set axis titles to black
  )

# Histogram of Original Data
p1 <- ggplot(Topsnph2, aes(x = h2)) +
  geom_histogram(binwidth = diff(range(Topsnph2$h2))/30, fill = npgc()[7], color = "black") +
  xlab("Heritability") +
  ylab("Frequency") +
  annotate("text", x = Inf, y = Inf, label = p_value_text, hjust = 1.1, vjust = 2, size = 5, color = "black") +
  scale_x_continuous(expand = c(0, 0)) +  # Start x-axis from 0
  scale_y_continuous(expand = c(0, 0)) +  # Start y-axis from 0
  common_theme 

# Histogram of Sample Data
p2 <- ggplot(Topsnph2[lab, ], aes(x = h2)) +
  geom_histogram(binwidth = diff(range(Topsnph2$h2[lab]))/30, fill = npgc()[11], color = "black") +
  xlab("Heritability") +
  ylab("Frequency") +
  scale_x_continuous(expand = c(0, 0)) +  # Start x-axis from 0
  scale_y_continuous(expand = c(0, 0)) +  # Start y-axis from 0
  common_theme

Topsnph2$Group <- "Original"
Topsnph2_sample <- Topsnph2[lab, ]
Topsnph2_sample$Group <- "Sample"

# Combine the original and sample data
combined_data <- rbind(Topsnph2, Topsnph2_sample)

# Plot with legend
p3 <- ggplot(combined_data, aes(x = h2, color = Group)) +
  geom_density(size = 1, adjust = 1) +
  xlab("Heritability") +
  ylab("Density") +
  scale_x_continuous(expand = c(0, 0)) +  # Start x-axis from 0
  scale_y_continuous(expand = c(0, 0)) +  # Start y-axis from 0
  scale_color_manual(values = c("Original" = npgc()[7], "Sample" = npgc()[11])) +  # Manually set colors
  common_theme+
  theme(legend.position = c(.85, .92), # Place the legend at the top
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

# Plotting all together using gridExtra
p <- plot_grid(p1, p2, p3, labels = c("a", "b", "c"), label_size = 16, ncol = 3)
p
