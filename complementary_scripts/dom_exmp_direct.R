library(ggplot2)

# 1. Dataset with explicit coordinates for 4 scenarios
df_clean <- data.frame(
  scenario = rep(c("βA(+) βD(+)\nEffect Allele Dominant", "βA(+) βD(-)\nNon-effect Allele Dominant", 
                   "βA(-) βD(+)\nNon-effect Allele Dominant", "βA(-) βD(-)\nEffect Allele Dominant"), each = 3),
  # Using 1, 2, 3 for x-axis to ensure math works for segments
  x_idx     = rep(c(1, 2, 3), 4), 
  geno_lab  = rep(c("AA", "AT", "TT"), 4),
  obs_y     = c(-1, 0.5, 1,   # Scenario 1
                -1, -0.5, 1,  # Scenario 2
                1, 0.5, -1,  # Scenario 3
                1, -0.5, -1), # Scenario 4
  additive_y = c(-1, 0, 1,    # Additive nulls
                 -1, 0, 1,
                 1, 0, -1,
                 1, 0, -1)
)

# Fix factor levels for correct X-axis order
df_clean$geno_lab <- factor(df_clean$geno_lab, levels = c("AA", "AT", "TT"))
bracket_data <- subset(df_clean, geno_lab == "AT")

# 2. Plotting 2x2
p <- ggplot(df_clean, aes(x = geno_lab, group = 1)) +
  geom_hline(yintercept = 0, color = "grey92", linewidth = 0.3) +
  # LAYER 1: The Additive Baseline (Dashed)
  geom_line(aes(y = additive_y), linetype = "dashed", color = "grey60", size = 0.5) +
  
  # LAYER 2: The Observed Relationship (Solid)
  geom_line(aes(y = obs_y), color = "black", size = 0.8) +
  
  geom_segment(data = subset(df_clean, geno_lab == "AT"),
               aes(x = 2, xend = 2, y = additive_y, yend = obs_y),
               color = "firebrick3", linetype = "dashed", linewidth = 0.8) +
  
  geom_text(data = subset(df_clean, geno_lab == "AT"),
            aes(x = 2.5, y = ((additive_y + obs_y) / 2)+0.05, label = "}"),
            color = "black", size = 12, family = "sans") +
  
  # Adjusted the "d" label to sit nicely outside the new bracket
  geom_text(data = subset(df_clean, geno_lab == "AT"),
            aes(x = 2.80, y = (additive_y + obs_y) / 2, label = "d"),
            color = "firebrick3", fontface = "italic", size = 4.5) +
  
  # LAYER 4: Genotype Points
  geom_point(aes(y = obs_y, fill = geno_lab), shape = 21, size = 4.5, stroke = 1) +
  
  # Color palette: White (Major), Red (Het), Black (Minor)
  scale_fill_manual(values = c("AA" = "white", "AT" = "firebrick3", "TT" = "black")) +
  
  # 2x2 FACETING
  facet_wrap(~scenario, nrow = 2) + 
  
  # Nature Genetics Styling
  labs(title = "The Direction of Dominance Deviations",
       y = "Phenotype Value (Standardized)", x = "Genotype") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0),
    strip.background = element_blank(),
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 10), hjust = 0.5),
    axis.title = element_text(face = "bold", size = 12, color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.spacing = unit(2, "lines") # Extra space between plots
  )


p

ggsave(filename = "/Users/sezgi/Documents/overlapped_SNPs/plots/dom_scenarios.png", plot = p,
  width = 170, 
  height = 150, 
  units = "mm",
  dpi = 600
)
