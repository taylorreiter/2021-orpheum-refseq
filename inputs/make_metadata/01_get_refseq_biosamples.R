# on 09/09/2021 downloaded information
library(dplyr)
library(janitor)
library(readr)
set.seed(1)
setwd("~/github/2021-orpheum-refseq")

refseq <- read_tsv("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt", 
                   skip = 2, col_names = "assembly_accession") %>%
  mutate(accession_number = gsub("GCF_", "", assembly_accession))

prok <- read_tsv("https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt") %>%
  clean_names() 

prok_refseq <- prok %>%
  mutate(accession_number = gsub("GCA_", "", assembly_accession)) %>%
  filter(accession_number %in% refseq$accession_number) 

write_tsv(prok_refseq[ , 'bio_sample_accession'], "inputs/make_metadata/refseq_prokaryotes_biosamples.txt", col_names = F)
