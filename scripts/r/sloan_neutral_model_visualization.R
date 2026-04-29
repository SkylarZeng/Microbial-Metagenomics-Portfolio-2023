library(readr)
library(dplyr)
library(ggplot2)

skin <- read_tsv("skin_to_doorknob_residence3.txt")
names(skin)[1] <- "Species"
skin <- skin %>% mutate(Partition = ifelse(freq < pred.lwr, "Below", "Neutral"))
skin <- skin %>% mutate(Partition = ifelse(freq > pred.upr, "Above", Partition))
write_tsv(skin, "skin_to_doorknob_residence3.txt")

spp <- read_tsv("kraken_report_all_table_rarefied_skin_residence3.txt")
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m / N
p <- as.data.frame(p)
colnames(p) <- "mean_abundance"

write.table(
  p,
  "kraken_species_mean_abundance_skin_residence3.txt",
  sep = "\t",
  row.names = TRUE
)

abundance <- read_tsv("kraken_species_mean_abundance_skin_residence3.txt")
names(abundance)[1] <- "Species"
merged <- left_join(skin, abundance)
write_tsv(merged, "skin_to_doorknob_residence3_merge.txt")

cbPalette <- c("#537bff", "#ff3c6d", "#b248ff")
Table <- read_tsv("skin_to_doorknob_residence3.txt") %>%
  filter(!Species == "unclassified")

Plot <- ggplot(Table, aes(x = log(mean_abundance)))
Plot <- Plot + geom_point(aes(y = freq, colour = Partition), size = 1.5, alpha = 0.8)
Plot <- Plot + geom_line(aes(y = freq.pred))
Plot <- Plot + geom_ribbon(aes(ymin = pred.lwr, ymax = pred.upr), fill = "#969696", alpha = 0.5)
Plot <- Plot + scale_color_manual(values = cbPalette)
Plot <- Plot + theme_classic()
Plot <- Plot + theme(legend.title = element_blank())
Plot <- Plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
Plot <- Plot + ggtitle("Residence 3")
Plot <- Plot + xlab("log (mean relative abundance of species on skin)")
Plot <- Plot + ylab("Occurrence frequency of species on doorknob")

ggsave("sloan_neutral_model_skin_to_doorknob_residence3.png", width = 12, height = 7)

