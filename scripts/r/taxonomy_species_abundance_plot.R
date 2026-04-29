library(devtools)
library(wilkoxmisc)
library(reshape2)
library(dplyr)
library(readr)
library(ggplot2)
library(starsExtra)

OTU <- read.csv("kraken_report_all_species.csv")

OTU <- OTU %>%
  gather(Sample, Count, 2:ncol(OTU)) %>%
  filter(Count > 0)

names(OTU)[1] <- "Species"
OTU$Sample <- gsub("kraken.report.", "", OTU$Sample)
OTU$Sample <- gsub("_kraken_report", "", OTU$Sample)

Meta <- read_tsv("metadata_doorknob_skin.txt")
OTU <- OTU %>% filter(OTU$Sample %in% Meta$Sample)

OTU <- ddply(OTU, .(Sample, Species), summarise, count = sum(Count))
OTU <- ddply(OTU, .(Sample), mutate, RelativeAbundance = (count * 100) / sum(count))

write_tsv(OTU, "kraken_species_abundance.tidy.txt")

OTU <- read_tsv("kraken_species_abundance.tidy.txt")
Meta <- read_tsv("metadata_doorknob_skin.txt")
Meta2 <- Meta %>% filter(Location == "Residence 3")
OTUTable <- merge(OTU, Meta2, by = "Sample", all.x = TRUE)
OTU3 <- OTUTable %>% filter(Location == "Residence 3")

Top12 <- collapse_taxon_table(OTU3, n = 12, Rank = "Species")
Top12species <- merge(Top12, Meta2, by = "Sample", all.x = TRUE)
write_tsv(Top12species, "Top12species.txt")

SpeciesAbundance <- read_tsv("Top12Species.txt")

Plot <- ggplot(
  SpeciesAbundance,
  aes(
    x = factor(
      Sampling_day,
      levels = c(
        "Day 1 day", "Day 1 night", "Day 2 day", "Day 2 night",
        "Day 3 day", "Day 3 night", "Day 4 day", "Day 4 night",
        "Day 5 day", "Day 5 night", "Day 6 day", "Day 6 night",
        "Day 7 day", "Day 7 night", "Day 8 day", "Day 8 night",
        "Day 9 day", "Day 9 night", "Day 10 day", "Day 10 night"
      )
    ),
    y = RelativeAbundance,
    fill = factor(
      Species,
      levels = c(
        "Chryseobacterium taklimakanense", "Cutibacterium acnes",
        "Dietzia lutea", "Dietzia sp. oral taxon 368",
        "Gordonia bronchialis", "Gordonia sp. KTR9", "Gordonia terrae",
        "Janibacter indicus", "Kytococcus sedentarius",
        "Micrococcus luteus", "Moraxella osloensis",
        "Minor/Unclassified"
      )
    )
  )
)

Plot <- Plot + geom_bar(stat = "identity")
Plot <- Plot + scale_fill_brewer(palette = "Paired")
Plot <- Plot + theme_classic()
Plot <- Plot + facet_grid(Sample_type ~ Location, scales = "free", space = "free")
Plot <- Plot + scale_y_continuous(expand = c(0, 0))
Plot <- Plot + ylab("Relative Abundance (%)") + xlab("Sampling day") + labs(fill = "Species")
Plot <- Plot + theme(axis.text.x = element_text(hjust = 0, angle = 90))

ggsave("taxonomy_general.png", width = 12, height = 7)
ggsave("top12species_by_sample_type.pdf", dpi = 1086)

