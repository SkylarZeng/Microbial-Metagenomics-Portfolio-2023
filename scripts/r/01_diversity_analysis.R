# =============================================================================
# 01_diversity_analysis.R
# Alpha-diversity analysis of microbial communities from Bracken outputs.
# Computes Shannon H', Simpson 1-D, Chao1 richness, and Pielou evenness.
# Statistical comparison between land-use groups (forest vs. agricultural).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(phyloseq)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
})

# ---------- paths (edit as needed) ----------
BRACKEN_DIR  <- "data/bracken/species"   # directory of per-sample Bracken tables
METADATA_CSV <- "data/metadata.csv"
OUTPUT_DIR   <- "figures"
# --------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load metadata ────────────────────────────────────────────────────────

metadata <- read_csv(METADATA_CSV, show_col_types = FALSE) %>%
  mutate(land_use = factor(land_use, levels = c("forest", "agricultural")))

# ── 2. Aggregate Bracken species tables into OTU matrix ─────────────────────

read_bracken <- function(path, sample_id) {
  read_tsv(path, show_col_types = FALSE) %>%
    select(name, new_est_reads) %>%
    mutate(sample = sample_id)
}

bracken_files <- list.files(BRACKEN_DIR, pattern = "_bracken_species\\.txt$",
                             full.names = TRUE)

otu_long <- map2_dfr(
  bracken_files,
  tools::file_path_sans_ext(basename(bracken_files)) %>%
    str_remove("_bracken_species"),
  read_bracken
)

otu_wide <- otu_long %>%
  pivot_wider(names_from = sample, values_from = new_est_reads,
              values_fill = 0) %>%
  column_to_rownames("name")

# samples as rows for vegan
otu_mat <- t(as.matrix(otu_wide))

# ── 3. Rarefy to even depth ──────────────────────────────────────────────────

set.seed(42)
min_reads <- min(rowSums(otu_mat))
message("Rarefying to: ", min_reads, " reads per sample")
otu_rarefied <- rrarefy(otu_mat, sample = min_reads)

# ── 4. Compute alpha-diversity indices ───────────────────────────────────────

alpha_df <- data.frame(
  sample_id  = rownames(otu_rarefied),
  shannon    = diversity(otu_rarefied, index = "shannon"),
  simpson    = diversity(otu_rarefied, index = "simpson"),
  richness   = specnumber(otu_rarefied),
  stringsAsFactors = FALSE
) %>%
  left_join(metadata, by = "sample_id") %>%
  mutate(
    evenness = shannon / log(richness)   # Pielou's J
  )

# ── 5. Statistical tests (Wilcoxon rank-sum) ────────────────────────────────

stats_results <- alpha_df %>%
  pivot_longer(cols = c(shannon, simpson, richness, evenness),
               names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  wilcox_test(value ~ land_use) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

print(stats_results)

# ── 6. Plot ──────────────────────────────────────────────────────────────────

plot_alpha <- function(df, yvar, ylabel, stat_tbl) {
  sig_row <- stat_tbl %>% filter(metric == yvar)
  p_label <- paste0("p = ", signif(sig_row$p.adj, 2),
                    " (", sig_row$p.adj.signif, ")")

  ggplot(df, aes(x = land_use, y = .data[[yvar]], fill = land_use)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    scale_fill_manual(values = c("forest" = "#2d6a4f", "agricultural" = "#b5838d")) +
    labs(x = NULL, y = ylabel, caption = p_label) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none")
}

p1 <- plot_alpha(alpha_df, "shannon",  "Shannon H'",   stats_results)
p2 <- plot_alpha(alpha_df, "simpson",  "Simpson 1-D",  stats_results)
p3 <- plot_alpha(alpha_df, "richness", "Species richness", stats_results)
p4 <- plot_alpha(alpha_df, "evenness", "Pielou evenness",  stats_results)

final_plot <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
                         labels = "AUTO")

ggsave(file.path(OUTPUT_DIR, "alpha_diversity.png"),
       final_plot, width = 8, height = 7, dpi = 300)

message("Saved: ", file.path(OUTPUT_DIR, "alpha_diversity.png"))

# ── 7. Save tidy results ────────────────────────────────────────────────────

write_csv(alpha_df,      "data/alpha_diversity_indices.csv")
write_csv(stats_results, "data/alpha_diversity_stats.csv")

sessionInfo()
