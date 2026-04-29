# =============================================================================
# 03_differential_abundance.R
# DESeq2-based differential abundance testing of microbial taxa between
# forest and agricultural soils, using Bracken read counts as input.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ---------- paths ----------
BRACKEN_DIR  <- "data/bracken/genus"   # genus-level for DA testing
METADATA_CSV <- "data/metadata.csv"
OUTPUT_DIR   <- "figures"
# ---------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load and assemble count matrix ───────────────────────────────────────

metadata <- read_csv(METADATA_CSV, show_col_types = FALSE) %>%
  mutate(land_use = factor(land_use, levels = c("forest", "agricultural"))) %>%
  arrange(sample_id)

read_bracken_genus <- function(path, sample_id) {
  read_tsv(path, show_col_types = FALSE) %>%
    select(name, new_est_reads) %>%
    mutate(sample = sample_id)
}

bracken_files <- list.files(BRACKEN_DIR, pattern = "_bracken_genus\\.txt$",
                             full.names = TRUE)
sample_ids <- tools::file_path_sans_ext(basename(bracken_files)) %>%
  str_remove("_bracken_genus")

count_long <- map2_dfr(bracken_files, sample_ids, read_bracken_genus)

count_mat <- count_long %>%
  pivot_wider(names_from = sample, values_from = new_est_reads,
              values_fill = 0L) %>%
  column_to_rownames("name") %>%
  as.matrix() %>%
  round()   # DESeq2 requires integer counts

# Subset and order columns to match metadata
count_mat <- count_mat[, metadata$sample_id]

# Remove taxa with zero counts in all samples
count_mat <- count_mat[rowSums(count_mat) > 0, ]

# ── 2. DESeq2 ────────────────────────────────────────────────────────────────

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = metadata,
  design    = ~ land_use
)

# Filter low-count taxa: keep genera with ≥ 10 reads in ≥ 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
message("Taxa retained after low-count filter: ", nrow(dds))

dds <- DESeq(dds, quiet = TRUE)

# Agricultural vs. forest (reference)
res <- results(dds,
               contrast  = c("land_use", "agricultural", "forest"),
               alpha     = 0.05,
               pAdjustMethod = "BH")

res_df <- as.data.frame(res) %>%
  rownames_to_column("taxon") %>%
  arrange(padj) %>%
  mutate(
    sig     = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
    label   = ifelse(sig & rank(padj) <= 20, taxon, NA_character_),
    direction = case_when(
      sig & log2FoldChange >  1 ~ "Enriched in agricultural",
      sig & log2FoldChange < -1 ~ "Enriched in forest",
      TRUE ~ "Not significant"
    )
  )

# ── 3. Volcano plot ──────────────────────────────────────────────────────────

p_volcano <- ggplot(res_df,
                    aes(x = log2FoldChange, y = -log10(padj),
                        colour = direction)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_text_repel(aes(label = label), size = 2.8,
                   max.overlaps = 20, show.legend = FALSE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(
    values = c("Enriched in agricultural" = "#b5838d",
               "Enriched in forest"       = "#2d6a4f",
               "Not significant"          = "grey70"),
    name = NULL
  ) +
  labs(
    x     = "log₂ fold-change (agricultural / forest)",
    y     = "-log₁₀(adjusted p-value)",
    title = "Differential Abundance by Land Use (DESeq2, genus level)"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "differential_abundance_volcano.png"),
       p_volcano, width = 8, height = 6, dpi = 300)

message("Saved: ", file.path(OUTPUT_DIR, "differential_abundance_volcano.png"))

# ── 4. Heatmap of top differentially abundant genera ────────────────────────

top_genera <- res_df %>%
  filter(sig) %>%
  slice_min(padj, n = 40) %>%
  pull(taxon)

if (length(top_genera) > 0) {
  vst_mat <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  heat_mat <- vst_mat[top_genera, ]

  ann_col <- metadata %>%
    select(sample_id, land_use) %>%
    column_to_rownames("sample_id") %>%
    as.data.frame()

  ann_colors <- list(
    land_use = c(forest = "#2d6a4f", agricultural = "#b5838d")
  )

  pheatmap(
    heat_mat,
    annotation_col  = ann_col,
    annotation_colors = ann_colors,
    scale           = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean",
    show_rownames   = TRUE,
    show_colnames   = FALSE,
    fontsize_row    = 8,
    filename        = file.path(OUTPUT_DIR, "differential_abundance_heatmap.png"),
    width = 9, height = 10
  )
  message("Saved: ", file.path(OUTPUT_DIR, "differential_abundance_heatmap.png"))
}

# ── 5. Save results ──────────────────────────────────────────────────────────

write_csv(res_df, "data/deseq2_results.csv")

sessionInfo()
