# =============================================================================
# 04_mag_functional_heatmap.R
# Visualise the distribution of key functional gene categories (KEGG Orthology)
# across medium/high-quality MAGs recovered from the metagenomes.
# Highlights nitrogen cycling, carbon fixation, and stress-response genes.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(patchwork)
})

# ---------- paths ----------
EGGNOG_DIR   <- "data/eggnog"          # per-MAG eggNOG-mapper output (.annotations)
CHECKM_CSV   <- "data/checkm2_all_samples_summary.csv"
GTDBTK_CSV   <- "data/gtdbtk_summary.csv"
OUTPUT_DIR   <- "figures"
# ---------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Target KEGG modules (nitrogen & carbon cycling + stress) ─────────────────

TARGET_KO <- tribble(
  ~ko_id,     ~gene,    ~pathway,
  "K02588",   "nifH",   "N fixation",
  "K02586",   "nifD",   "N fixation",
  "K02591",   "nifK",   "N fixation",
  "K00370",   "narG",   "Denitrification",
  "K00371",   "narH",   "Denitrification",
  "K02567",   "napA",   "Denitrification",
  "K04561",   "norB",   "Denitrification",
  "K00376",   "nosZ",   "Denitrification",
  "K10944",   "amoA",   "Nitrification",
  "K10945",   "amoB",   "Nitrification",
  "K10946",   "amoC",   "Nitrification",
  "K00399",   "mcrA",   "Methanogenesis",
  "K00401",   "mcrB",   "Methanogenesis",
  "K14067",   "mmoX",   "Methanotrophy",
  "K00194",   "cooS",   "CO oxidation",
  "K01601",   "rbcL",   "C fixation (CBB)",
  "K01602",   "rbcS",   "C fixation (CBB)",
  "K14468",   "aclA",   "C fixation (rTCA)",
  "K03737",   "porA",   "C fixation (rTCA)"
)

# ── 1. Parse eggNOG-mapper annotations ───────────────────────────────────────

parse_eggnog <- function(path) {
  read_tsv(path,
           comment   = "#",
           col_names = TRUE,
           show_col_types = FALSE) %>%
    select(query = `#query`, KEGG_ko) %>%
    filter(!is.na(KEGG_ko), KEGG_ko != "-") %>%
    separate_rows(KEGG_ko, sep = ",") %>%
    mutate(KEGG_ko = str_trim(KEGG_ko),
           mag_id   = tools::file_path_sans_ext(basename(path)))
}

eggnog_files <- list.files(EGGNOG_DIR, pattern = "\\.annotations$",
                            full.names = TRUE)

annotation_df <- map_dfr(eggnog_files, parse_eggnog)

# ── 2. Presence/absence matrix for target KOs ────────────────────────────────

pa_mat <- annotation_df %>%
  filter(KEGG_ko %in% TARGET_KO$ko_id) %>%
  distinct(mag_id, KEGG_ko) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = KEGG_ko, values_from = present,
              values_fill = 0L) %>%
  column_to_rownames("mag_id")

# Add any missing target KOs as all-zero columns
missing_ko <- setdiff(TARGET_KO$ko_id, colnames(pa_mat))
pa_mat[missing_ko] <- 0L

pa_mat <- pa_mat[, TARGET_KO$ko_id]   # enforce column order

# ── 3. Load CheckM2 and GTDB-Tk metadata ────────────────────────────────────

checkm  <- read_csv(CHECKM_CSV,  show_col_types = FALSE)
gtdbtk  <- read_csv(GTDBTK_CSV, show_col_types = FALSE) %>%
  mutate(phylum = str_extract(classification, "(?<=p__)[^;]+"))

mag_meta <- checkm %>%
  left_join(gtdbtk %>% select(user_genome, phylum, classification),
            by = c("mag_id" = "user_genome")) %>%
  filter(mag_id %in% rownames(pa_mat))

# ── 4. ComplexHeatmap ────────────────────────────────────────────────────────

row_order <- mag_meta %>%
  arrange(phylum, desc(completeness)) %>%
  pull(mag_id)

row_anno <- mag_meta %>%
  select(mag_id, completeness, contamination, phylum) %>%
  column_to_rownames("mag_id")

row_anno_ht <- rowAnnotation(
  Completeness  = anno_barplot(row_anno[row_order, "completeness"],
                                gp = gpar(fill = "#52b788")),
  Contamination = anno_barplot(row_anno[row_order, "contamination"],
                                gp = gpar(fill = "#e76f51")),
  Phylum = anno_text(row_anno[row_order, "phylum"],
                      gp = gpar(fontsize = 7))
)

col_split <- TARGET_KO$pathway
col_labels <- TARGET_KO$gene

ht <- Heatmap(
  as.matrix(pa_mat[row_order, ]),
  name                  = "Presence",
  col                   = c("0" = "#f0f0f0", "1" = "#1d3557"),
  column_split          = col_split,
  column_labels         = col_labels,
  row_title             = NULL,
  show_row_names        = FALSE,
  cluster_rows          = FALSE,
  cluster_columns       = FALSE,
  column_names_gp       = gpar(fontsize = 9),
  left_annotation       = row_anno_ht,
  heatmap_legend_param  = list(
    title = "Gene present",
    labels = c("absent", "present"),
    at = c(0, 1)
  )
)

png(file.path(OUTPUT_DIR, "functional_heatmap.png"),
    width = 3000, height = 2400, res = 300)
draw(ht, merge_legend = TRUE)
dev.off()

message("Saved: ", file.path(OUTPUT_DIR, "functional_heatmap.png"))

# ── 5. MAG quality scatter ───────────────────────────────────────────────────

p_qual <- ggplot(checkm,
                 aes(x = completeness, y = contamination,
                     colour = completeness >= 90 & contamination <= 5)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_hline(yintercept = c(5, 10), linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = c(70, 90), linetype = "dashed", colour = "grey60") +
  scale_colour_manual(
    values = c("TRUE" = "#52b788", "FALSE" = "#adb5bd"),
    labels = c("TRUE" = "High quality", "FALSE" = "Medium quality"),
    name   = NULL
  ) +
  labs(
    x     = "Completeness (%)",
    y     = "Contamination (%)",
    title = "MAG Quality (CheckM2)"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "mag_completeness_contamination.png"),
       p_qual, width = 6, height = 5, dpi = 300)

message("Saved: ", file.path(OUTPUT_DIR, "mag_completeness_contamination.png"))

sessionInfo()
