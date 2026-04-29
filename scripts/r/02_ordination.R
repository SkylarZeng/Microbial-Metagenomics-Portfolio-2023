# =============================================================================
# 02_ordination.R
# Beta-diversity analysis: Bray-Curtis dissimilarity, NMDS/PCoA ordination,
# and PERMANOVA testing of community composition by land-use group.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggplot2)
  library(ggrepel)
})

# ---------- paths ----------
BRACKEN_DIR  <- "data/bracken/species"
METADATA_CSV <- "data/metadata.csv"
OUTPUT_DIR   <- "figures"
# ---------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load data (re-uses rarefied matrix from 01_diversity_analysis.R) ─────
# Source the helper only if the object does not exist
if (!exists("otu_rarefied")) {
  source("scripts/r/01_diversity_analysis.R")
}

metadata <- read_csv(METADATA_CSV, show_col_types = FALSE) %>%
  mutate(land_use = factor(land_use, levels = c("forest", "agricultural")))

# Ensure metadata row order matches otu_rarefied
meta_ordered <- metadata %>%
  filter(sample_id %in% rownames(otu_rarefied)) %>%
  arrange(match(sample_id, rownames(otu_rarefied)))

# ── 2. Bray-Curtis distance ──────────────────────────────────────────────────

bc_dist <- vegdist(otu_rarefied, method = "bray")

# ── 3. PERMANOVA ─────────────────────────────────────────────────────────────

set.seed(42)
permanova_result <- adonis2(
  bc_dist ~ land_use,
  data  = meta_ordered,
  permutations = 9999,
  by = "margin"
)
print(permanova_result)

# PERMDISP (test for homogeneity of dispersions)
disp <- betadisper(bc_dist, meta_ordered$land_use)
disp_test <- permutest(disp, permutations = 9999)
print(disp_test)

# ── 4. PCoA ordination ───────────────────────────────────────────────────────

pcoa <- cmdscale(bc_dist, eig = TRUE, k = 2)
pct_var <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)

pcoa_df <- data.frame(
  sample_id = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  left_join(meta_ordered, by = "sample_id")

r2   <- round(permanova_result$R2[1], 3)
pval <- permanova_result$`Pr(>F)`[1]
anno <- paste0("PERMANOVA: R² = ", r2, ", p = ", pval)

p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2,
                               colour = land_use, shape = land_use)) +
  geom_point(size = 3.5, alpha = 0.9) +
  stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.8) +
  geom_text_repel(aes(label = sample_id), size = 2.5,
                   max.overlaps = 10, show.legend = FALSE) +
  scale_colour_manual(values = c("forest" = "#2d6a4f",
                                  "agricultural" = "#b5838d"),
                       name = "Land use") +
  scale_shape_manual(values = c("forest" = 16, "agricultural" = 17),
                      name = "Land use") +
  labs(
    x = paste0("PC1 (", pct_var[1], "% variance)"),
    y = paste0("PC2 (", pct_var[2], "% variance)"),
    title = "Bray-Curtis PCoA — Community Composition",
    caption = anno
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "beta_diversity_pcoa.png"),
       p_pcoa, width = 7, height = 6, dpi = 300)

message("Saved: ", file.path(OUTPUT_DIR, "beta_diversity_pcoa.png"))

# ── 5. NMDS (stress diagnostic) ──────────────────────────────────────────────

set.seed(42)
nmds <- metaMDS(otu_rarefied, distance = "bray", k = 2,
                trymax = 100, trace = 0)
message("NMDS stress: ", round(nmds$stress, 4))

nmds_df <- data.frame(
  sample_id = rownames(nmds$points),
  NMDS1 = nmds$points[, 1],
  NMDS2 = nmds$points[, 2]
) %>%
  left_join(meta_ordered, by = "sample_id")

p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2,
                               colour = land_use, shape = land_use)) +
  geom_point(size = 3.5, alpha = 0.9) +
  stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.8) +
  scale_colour_manual(values = c("forest" = "#2d6a4f",
                                  "agricultural" = "#b5838d"),
                       name = "Land use") +
  scale_shape_manual(values = c("forest" = 16, "agricultural" = 17),
                      name = "Land use") +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("Stress = ", round(nmds$stress, 3)),
           hjust = 1.1, vjust = -0.5, size = 3.5) +
  labs(title = "Bray-Curtis NMDS — Community Composition",
       caption = anno) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "beta_diversity_nmds.png"),
       p_nmds, width = 7, height = 6, dpi = 300)

# ── 6. Save distance matrix ──────────────────────────────────────────────────

saveRDS(bc_dist,  "data/bray_curtis_dist.rds")
write_csv(pcoa_df, "data/pcoa_coordinates.csv")

sessionInfo()
