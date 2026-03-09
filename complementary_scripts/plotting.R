library(data.table)
library(tidyverse)
library(ggrepel)
library(viridis)
library(scales)

rm(list = setdiff(ls(all.names = TRUE), c("all_sig_df", "clumped_df")))

#  Load the data 
all_sig_path <- "/Users/sezgi/Documents/overlapped_SNPs/significant_SNPs/all_sig_SNPs.tsv.gz"
all_sig_df <- fread(all_sig_path, select = c("variant", "rsid", "chr", "pos", "add_sig_total", "dom_sig_total"),
                    colClasses = c(
                      variant = "character",
                      rsid    = "character",
                      chr     = "numeric",
                      pos     = "numeric",
                      add_sig_total = "numeric",
                      dom_sig_total = "numeric"
                    ))

clumped_path <- "/Users/sezgi/Documents/overlapped_SNPs/significant_SNPs/clumped_SNPs.tsv.gz"
clumped_df <- fread(clumped_path, select = c("variant", "rsid", "chr", "pos", "add_sig_total", 
                                          "dom_sig_total", "sig_traits", "trait_count_dom",  
                                          "traits_dom_pval", "avg_dom_pval", "is_subsumed",
                                          "subsumed_by"),
                    colClasses = c(
                      variant = "character",
                      rsid    = "character",
                      chr     = "numeric",
                      pos     = "numeric",
                      add_sig_total = "numeric",
                      dom_sig_total = "numeric",
                      avg_dom_pval = "numeric"
                    ))


# Replace rsID with chr:pos:ref:alt
nrow(all_sig_df[is.na(rsid)])
all_sig_df[is.na(rsid) | rsid == "", rsid := variant]
setorder(all_sig_df, chr, pos)

# Adding the subsumed information
all_sig_df[clumped_df, on = .(variant), is_subsumed := i.is_subsumed]
all_sig_df[clumped_df, on = .(variant), sig_traits := i.sig_traits]

# Prepare the dataset
man_df <- all_sig_df %>% 
  
  # chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # cumulative positions of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(all_sig_df, ., by=c("chr"="chr")) %>%
  
  # cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate(BPcum=pos+tot)


# Center chr positions
axisdf <- man_df %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


# Manhattan Plot
ggplot(man_df, aes(x=BPcum, y=dom_sig_total)) +
  
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0.5, 5)), 
                     limits = c(0, NA)) +  
  
  # Custom the theme
  theme_classic(base_size = 14) +
  theme( 
    legend.position = "none",
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(color = "black", size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 15, face = "bold"),
    
    plot.tag.position = c(1, 1), 
    plot.tag = element_text(size = 12, face = "bold", hjust = 1, vjust = 2.85),
    
    # Add a faint dashed line strictly for the Y-axis to help the reader gauge peak height
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Pleiotropic Variants",
    subtitle = "Dominance",
    tag = bquote(bold(italic(p) < 4.72 %*% 10^{-11})),
    x = "Chromosome",
    y = "Number of Significant Traits",
  )



# 1. Identify the top 4 SNPs based on the number of significant traits
zoom_df <- man_df %>% filter(chr %in% c(1, 6, 17, 19))

top_4_snps <- zoom_df %>%
  filter(is_subsumed == FALSE) %>%
  group_by(chr) %>%
  slice_max(order_by = dom_sig_total, n = 1, with_ties = FALSE) %>%
  ungroup()

# 2. Define the window size (250 Kb = 250,000 bp)
window_kb <- 1000000

# 3. Filter man_df for regions around these 4 SNPs
# We create a Locus label for each to ensure clean faceting
top_regions_df <- top_4_snps %>%
  split(.$variant) %>%
  map_dfr(function(lead) {
    man_df %>%
      filter(chr == lead$chr, 
             pos >= (lead$pos - window_kb), 
             pos <= (lead$pos + window_kb)) %>%
      mutate(Locus = paste0("Chromose ", lead$chr))
  })

# 4. Plot the 4 Loci
ggplot(top_regions_df, aes(x = pos / 1e6, y = dom_sig_total)) + 
  
  # Base points colored by chromosome (alternating colors logic)
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 3) +
  
  # Highlight the specific Lead SNP in each facet
  geom_point(data = subset(top_regions_df, variant %in% top_4_snps$variant), 
             shape = 21, color = "firebrick", size = 4, stroke = 1.2) +
  
  geom_label_repel( data=subset(top_regions_df, is_subsumed=="FALSE"), aes(label=rsid), size=3) +
  
  # Facet by Locus (independent x-axes for each 500Kb window)
  facet_wrap(~ Locus, scales = "free_x", ncol = 2) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0.5, 5)), 
                     limits = c(0, NA)) +  
  
  # Use 3-4 breaks to keep the Mb x-axis clean in small windows
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  
  theme_classic(base_size = 14) +
  theme( 
    legend.position = "none",
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold"),
    
    # Facet header styling
    strip.background = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    
    plot.title = element_text(size = 18, face = "bold"),
    plot.tag.position = c(1, 1), 
    plot.tag = element_text(size = 12, face = "bold", hjust = 1, vjust = 2.85),
    
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed")
  ) +
  labs(
    title = "Top 4 Pleiotropic Loci",
    subtitle = "Dominance (\u00B1 250 Kb Zoom)",
    x = "Genomic Position (Mb)", 
    y = "Number of Significant Traits"
  )

















# Manhattan plot
ggplot(man_df, aes(x=BPcum, y=dom_sig_total)) +
  
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0.5, 5)), 
                     limits = c(0, NA)) +  
  
  # Independent SNPs
  geom_point(data = subset(man_df, is_subsumed == "FALSE"), 
             shape = 21,          
             fill = "skyblue",   
             color = "firebrick",
             size = 3,
             stroke = 1.1) +
  # Annotation
  geom_label_repel( data=subset(man_df, is_subsumed == "FALSE"), aes(label=rsid), size=2) +
  
  # Custom the theme
  theme_classic(base_size = 14) +
  theme( 
    legend.position = "none",
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(color = "black", size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 15, face = "bold"),
    
    plot.tag.position = c(1, 1), 
    plot.tag = element_text(size = 12, face = "bold", hjust = 1, vjust = 2.85),
    
    # Add a faint dashed line strictly for the Y-axis to help the reader gauge peak height
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Pleiotropic Variants",
    subtitle = "Dominance",
    tag = bquote(bold(italic(p) < 4.72 %*% 10^{-11})),
    x = "Chromosome",
    y = "Number of Significant Traits",
  )




only_clumped <- clumped_df[is_subsumed == FALSE]


