library(dplyr)
library(janitor)
library(readr)
set.seed(1)
setwd("~/github/2021-orpheum-refseq")

metadata <- read_tsv("inputs/metadata.tsv")

# remove 3 failed samples -- samples weren't detected as properly formated paired-end
# sequences by khmer
failed <- c("DRR164904", "ERR025451", "SRR11451571", "SRR2566971", "ERR1760467",
            "DRR183252", "DRR239351")

metadata <- metadata %>%
  filter(!Run %in% failed) 

write_tsv(metadata, "inputs/metadata.tsv")
