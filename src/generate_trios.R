#!/usr/bin/Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)

sample_info_file <- snakemake@input[["sample_info"]]
fai_file <- snakemake@input[["fai"]]

# dev
# sample_info_file <- "data/trios_gt.csv"
# fai_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai"

# make the regions file from the reference genome
fai_names <- c("name", "length", "offset", "line_bases", "line_bytes")

fai <- fread(fai_file, col.names = fai_names)

inheritance <- fai[startsWith(name, "NC_"), .(
  CHROM = paste0(name, ":1-", length),
  ploidy = "M/M + F > M/F")]

# make the trios file
sample_info <- fread(sample_info_file)
sample_info[, c("hive", "caste", "indiv") := tstrsplit(sample, "_")]

# make all possible trios
# mother == queen
# father == drone
# child  == worker
mothers <- sample_info[caste == "queen", unique(sample)]
fathers <- sample_info[caste == "drone", unique(sample)]
children <- sample_info[caste == "worker", unique(sample)]

all_trios <- unique(data.table(expand.grid(mothers,
                              fathers,
                              children,
                              stringsAsFactors = FALSE)))

fwrite(all_trios,
       snakemake@output[["trios"]],
       sep = ",",
       col.names = FALSE)

fwrite(inheritance,
       snakemake@output[["rules"]],
       sep = "\t",
       col.names = FALSE)

sessionInfo()