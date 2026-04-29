library(reshape2)
library(wilkoxmisc)
library(ggplot2)
library(dplyr)
library(cultevo)
library(readr)
library(vegan)
library(MASS)
library(cowplot)

tax_cast <- read.table("kraken_species_abundance.cast.txt", header = TRUE, row.names = 1)
distance <- vegdist(tax_cast, method = "bray")
Matrix <- as.matrix(distance)
write.matrix(Matrix, sep = " ", "bray_curtis_distance_matrix.txt")

UniFrac <- read_tsv("bray_curtis_distance_matrix.txt")
UniFrac <- melt(UniFrac, value.name = "Distance")
names(UniFrac)[1:2] <- c("Sample1", "Sample2")

Samples <- read_tsv("metadata_doorknob_skin.txt") %>%
  select(Sample, Sample_type, Location, Sampling_day)

names(Samples)[1] <- "Sample2"
UniFrac <- merge(UniFrac, Samples, by = "Sample2", all.x = FALSE)
names(UniFrac)[4:6] <- c("Sample_type2", "Location2", "Sampling_day2")

names(Samples)[1] <- "Sample1"
UniFrac <- merge(UniFrac, Samples, by = "Sample1", all.x = FALSE)
names(UniFrac)[7:9] <- c("Sample_type1", "Location1", "Sampling_day1")

time <- read_tsv("sampling_time_hour_conversion.txt")
names(time)[1] <- "Sampling_day1"
merged <- left_join(UniFrac, time)
names(merged)[10] <- "Time_hour1"

names(time)[1] <- "Sampling_day2"
merged <- left_join(merged, time)
names(merged)[11] <- "Time_hour2"

UniFrac <- merged
UniFrac <- UniFrac[which(!UniFrac$Sample1 == UniFrac$Sample2), ]
UniFrac$HouseholdType <- ifelse(UniFrac$Location1 == UniFrac$Location2, "Same household", "Different household")

table <- UniFrac %>%
  filter(HouseholdType == "Same household") %>%
  filter(Sample_type2 == "Door knob") %>%
  filter(Sample_type1 %in% c("Right palm", "Left palm"))

table$Group <- paste0(table$Sample_type2, " vs. ", table$Sample_type1)
table <- table %>% mutate(time_decay = Time_hour1 - Time_hour2)

Plot1 <- ggplot(table, aes(x = time_decay, y = Distance))
Plot1 <- Plot1 + geom_point(color = "#455def", alpha = 0.4)
Plot1 <- Plot1 + geom_smooth(color = "#fd017d", fill = "#fd017d")
Plot1 <- Plot1 + theme_bw()
Plot1 <- Plot1 + geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5)
Plot1 <- Plot1 + ylab("Bray-Curtis dissimilarity between doorknob and palm")
Plot1 <- Plot1 + xlab("Time interval (hours)")
Plot1 <- Plot1 + theme(legend.position = "none")

cbPalette <- c("#590A30", "#90AA3C", "#EF6125")
P1 <- ggplot(table, aes(x = Distance, color = Group))
P1 <- P1 + geom_line(stat = "density", size = 2)
P1 <- P1 + scale_y_continuous(expand = c(0, 0), limits = c(0, 5))
P1 <- P1 + theme_classic()
P1 <- P1 + facet_wrap(~Location1)
P1 <- P1 + theme(panel.background = element_rect(colour = "Black"))
P1 <- P1 + scale_colour_manual(values = cbPalette)
P1 <- P1 + xlab("Normalized Bray Curtis Dissimilarity") + ylab("Density (%)")
P1 <- P1 + theme(legend.title = element_blank(), legend.position = "bottom")

Figure <- plot_grid(Plot1, P1, labels = LETTERS[1:2], ncol = 2)

ggsave(plot = Plot1, "Bray-Curtis_dissimilarity_between_doorknob_and_palm.png", width = 5, height = 5)
ggsave(plot = Figure, "Normalized_Bray_Curtis_Dissimilarity_combined_by_sample_type.png", width = 10, height = 5)
write_tsv(table, "doorkob_palm_comparison_Bray_Curtis.txt")

