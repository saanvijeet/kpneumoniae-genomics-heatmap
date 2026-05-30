# ============================================================
# Consolidated Resistance Gene Family Prevalence across High-Risk STs
# Optimized for Poster Readability (Compressed Row Count)
# ============================================================

# ---- 1. Setup ----
library(tidyverse)
library(viridis)

# ---- 2. Load data ----
kleborate <- read_csv("pathogenwatch-kpneumoniae-i5WstZ5MxAKsB8n4G1AXHQ-kleborate.csv")

# ---- 3. Define gene categories ----
gene_columns <- c(
  "Bla_Carb_acquired",       
  "Bla_ESBL_acquired",       
  "Bla_acquired",            
  "AGly_acquired",           
  "Flq_acquired",            
  "Sul_acquired",            
  "Tmt_acquired"             
)

# ---- 4. Parse gene calls into long format ----
genes_long <- kleborate %>%
  select(`Genome ID`, ST, all_of(gene_columns)) %>%
  pivot_longer(
    cols = all_of(gene_columns),
    names_to = "gene_class",
    values_to = "genes"
  ) %>%
  filter(genes != "-", !is.na(genes)) %>%
  separate_rows(genes, sep = ";") %>%
  mutate(
    gene = str_trim(genes),
    gene_clean = str_remove_all(gene, "\\.v\\d+"),
    gene_clean = str_remove_all(gene_clean, "[\\*\\?\\^]"),
    gene_clean = str_trim(gene_clean)
  ) %>%
  filter(gene_clean != "")

# ---- 5. Group individual variants into readable Gene Families ----
gene_classes <- genes_long %>%
  mutate(
    drug_class = case_when(
      gene_class == "Bla_Carb_acquired" ~ "Carbapenems",
      gene_class == "Bla_ESBL_acquired" ~ "Cephalosporins (ESBL)",
      gene_class == "Bla_acquired"      ~ "Other β-lactams",
      gene_class == "AGly_acquired"     ~ "Aminoglycosides",
      gene_class == "Flq_acquired"      ~ "Fluoroquinolones",
      gene_class == "Sul_acquired"      ~ "Sulphonamides",
      gene_class == "Tmt_acquired"      ~ "Trimethoprim",
      TRUE ~ "Other"
    ),
    # Collapse specific allele variants into broader diagnostic groups
    gene_family = case_when(
      str_detect(gene_clean, "^NDM") ~ "blaNDM-type",
      str_detect(gene_clean, "^KPC") ~ "blaKPC-type",
      str_detect(gene_clean, "^OXA-181|^OXA-48") ~ "blaOXA-48-like",
      str_detect(gene_clean, "^CTX-M") ~ "blaCTX-M-type",
      str_detect(gene_clean, "^OXA-1|^OXA-9") ~ "blaOXA-1/9-like",
      str_detect(gene_clean, "^strA|^strB") ~ "strA/strB",
      str_detect(gene_clean, "^aac\\(6'\\)-Ib-cr") ~ "aac(6')-Ib-cr",
      str_detect(gene_clean, "^sul1") ~ "sul1",
      str_detect(gene_clean, "^sul2") ~ "sul2",
      str_detect(gene_clean, "^dfrA14") ~ "dfrA14",
      str_detect(gene_clean, "^dfrA12") ~ "dfrA12",
      TRUE ~ "Other_Dropped"
    )
  ) %>%
  # Filter out background noise and rare grouped elements
  filter(gene_family != "Other_Dropped")

# ---- 6. Calculate prevalence per ST ----
st_totals <- kleborate %>%
  count(ST, name = "n_total")

prevalence <- gene_classes %>%
  distinct(`Genome ID`, ST, drug_class, gene_family) %>%
  count(ST, drug_class, gene_family, name = "n_with_gene") %>%
  left_join(st_totals, by = "ST") %>%
  mutate(prevalence = 100 * n_with_gene / n_total)

# ---- 7. Complete grid for true implicit 0% values ----
heatmap_data <- prevalence %>%
  complete(ST, nesting(drug_class, gene_family), fill = list(prevalence = 0, n_with_gene = 0))

# ---- 8. Order Sequence Types and Drug Classes ----
st_order <- c("ST11", "ST15", "ST101", "ST147", "ST258", "ST307")
heatmap_data$ST <- factor(heatmap_data$ST, levels = st_order)

drug_class_order <- c(
  "Carbapenems", "Cephalosporins (ESBL)", "Other β-lactams",
  "Aminoglycosides", "Fluoroquinolones", "Sulphonamides", "Trimethoprim"
)
heatmap_data$drug_class <- factor(heatmap_data$drug_class, levels = drug_class_order)

# Order rows strictly by custom priority appearance
family_order <- c(
  "blaNDM-type", "blaKPC-type", "blaOXA-48-like", 
  "blaCTX-M-type", 
  "blaOXA-1/9-like", 
  "strA/strB", 
  "aac(6')-Ib-cr", 
  "sul1", "sul2", 
  "dfrA14", "dfrA12"
)
heatmap_data$gene_family <- factor(heatmap_data$gene_family, levels = family_order)

# ---- 9. Generate mathematical expressions for Y-axis italics ----

# 1. Iterates through all gene family factor levels using lapply().
# 2. Uses conditional string matching (startsWith) to identify beta-lactamase ('bla') genes.
# 3. For 'bla' genes: Isolates the 'bla' prefix, applies an italic font style to it, 
#    and attaches the specific allele/family suffix (e.g., KPC-type) in standard roman text
#    using bquote() and the plotmath non-spacing glue operator (*).
# 4. For non-bla genes (e.g., sul1, strA/strB): Applies a uniform italic style across the entire string.
# 5. Maps these rendered expressions back to their plain-text data keys using setNames().
gene_labels <- setNames(
  lapply(levels(heatmap_data$gene_family), function(x) {
    if (startsWith(x, "bla")) {
      base_name <- str_remove(x, "bla")
      bquote(italic("bla") * .(base_name))
    } else {
      bquote(italic(.(x)))
    }
  }),
  levels(heatmap_data$gene_family)
)

# ---- 10. Build the high-legibility heatmap ----
heatmap_plot <- ggplot(
  heatmap_data,
  aes(x = ST, y = gene_family, fill = prevalence)
) +
  geom_tile(colour = "white", linewidth = 0.8) + 
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.15))) + 
  scale_y_discrete(labels = gene_labels) + 
  scale_fill_viridis_c(
    name   = "Prevalence (%)",
    option = "viridis",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  facet_grid(
    drug_class ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +
  labs(
    x = "Sequence type (n = 50 per ST)", 
    y = NULL,
    title = bquote(atop(
  bold("Key Resistance Family Prevalence"),
  bold("Across High-Risk ") * bolditalic("K. pneumoniae") * bold(" Lineages")
)),
    subtitle = "Consolidated Analysis of 300 Genomes from PathogenWatch"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x          = element_text(face = "bold", size = 13, margin = margin(t = 6)),
    axis.text.y          = element_text(size = 13, colour = "black"), 
    panel.grid           = element_blank(),
    strip.placement      = "outside",
    strip.text.y.left    = element_text(angle = 0, face = "bold", hjust = 0, size = 16, margin = margin(r = 16)),
    strip.background     = element_blank(),
    
    # Optimized Legend Properties for High-Visibility Poster Contrast
    legend.position      = "right",
    legend.title         = element_text(face = "bold", size = 14, margin = margin(b = 8)), 
    legend.text          = element_text(face = "bold", size = 13),                       
    legend.key.height    = unit(1.6, "cm"),                                               
    legend.key.width     = unit(0.7, "cm"),                                               
    
    # Centers title perfectly across the whole image canvas (including left text margins)
    plot.title.position  = "plot", 
    
    plot.title           = element_text(face = "bold", size = 20, hjust = 0.5, lineheight = 1.1), 
    plot.subtitle        = element_text(size = 13, colour = "grey35", face = "bold", hjust = 0.5, margin = margin(t = 5, b = 10)),
    panel.spacing.y      = unit(0.7, "lines"),
    plot.margin          = margin(t = 15, r = 15, b = 15, l = 15) 
  )

# ---- 11. Save at smaller, optimized presentation dimensions ----
graphics.off() 

print(heatmap_plot)
ggsave("kp_resistance_heatmap_clean.png", heatmap_plot, width = 11.0, height = 6.8, dpi = 300)
