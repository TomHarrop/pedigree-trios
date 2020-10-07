#!/usr/bin/Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(ggplot2)
library(data.table)

trios_file <- snakemake@input[["trios"]]

trios <- fread(trios_file)

trios[, c("queen", "drone", "worker") := 
        tstrsplit(`[4]Trio (mother,father,child)`, ",")]
trios[, freq_err := `[2]nBad`/(`# [1]nOK` + `[2]nBad`)]

# put in metadata
trios[, paste("queen", c("hive", "caste", "indiv"), sep = "_") :=
        tstrsplit(queen, "_")]
trios[, paste("drone", c("hive", "caste", "indiv"), sep = "_") :=
        tstrsplit(drone, "_")]
trios[, paste("worker", c("hive", "caste", "indiv"), sep = "_") :=
        tstrsplit(worker, "_")]

trios[, drone_hive_plot := paste("Drone from", drone_hive)]
trios[, queen_hive_plot := paste("Queen from", queen_hive)]

# plot
gp <- ggplot(trios,
       aes(x = worker, y = freq_err, colour = worker_hive)) +
  theme_grey(base_size = 6,
             base_family = "Lato") +
  theme(legend.key.size = unit(0.5, "lines")) +
  coord_flip() + 
  xlab(NULL) + ylab("Mendelian error frequency") +
  scale_colour_viridis_d(guide = guide_legend(title = "Worker hive")) +
  facet_grid(drone_hive_plot ~ queen_hive_plot,
             scales = "free_y") +
  geom_point(shape = 16,
             size = 1,
             position = position_jitter(width = 0.1),
             alpha = 0.75)

fh <- grid::convertUnit(
  grid::unit(227, "pt"),
  "in",
  valueOnly = TRUE)
fw <- grid::convertUnit(
  grid::unit(398, "pt"),
  "in",
  valueOnly = TRUE)

ggsave(snakemake@output[["plot"]], 
       gp,
       device = cairo_pdf,
       width = fw,
       height = fh,
       units = "in")


sessionInfo()