library(wilkoxmisc)
library(reshape2)
library(dplyr)
library(readr)
library(ggplot2)
library(starsExtra)
library(tidyverse)
library(readxl)

Meta <- metadata_doorknob_skin
Meta2 <- Meta %>% filter(Location == "Residence 3")
Match <- match
MergeData <- merge(Meta2, Match, by = "Sample", all.x = TRUE)
write_tsv(MergeData, "Information.txt")

coverm_output <- coverm_results

plot <- ggplot(data = coverm_output)
plot <- plot + geom_point(
  mapping = aes(
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
    y = Sample_type,
    colour = Sample_type,
    size = RelativeAbundance
  )
)

plot <- plot + labs(
  x = "Sampling day",
  y = "Sample type",
  colour = "Sample type",
  size = "RelativeAbundance(%)"
)
plot <- plot + facet_wrap(~classification)
plot <- plot + theme(axis.text.x = element_text(hjust = 0, angle = 90))

ggsave("coverM_transmission.png", width = 12, height = 7)

