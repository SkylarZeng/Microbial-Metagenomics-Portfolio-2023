library(ggplot2)
library(dplyr)

dRep_rMAGs <- read.csv("dRep_rMAGs.csv")
dRep_rMAGs_filter <- dRep_rMAGs %>% filter(Genotype == "Cluster A")

plot <- ggplot(data = dRep_rMAGs_filter)
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
    colour = factor(Sample_type),
    size = 2
  )
)

plot <- plot + labs(x = "Sampling day", y = "Sample type", colour = "Sample type")
cbPalette <- c("#7998ff", "#ff6780", "#c371ff")
plot <- plot + scale_color_manual(values = cbPalette)
plot <- plot + facet_wrap(~classification)
plot <- plot + theme(axis.text.x = element_text(hjust = 0, angle = 90))
plot <- plot + guides(size = "none")

ggsave("taxonomy_day.png", width = 12, height = 7)

