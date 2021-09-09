# on 09/09/2021 downloaded information
library(dplyr)
library(janitor)
library(readr)
set.seed(1)
setwd("~/github/2021-orpheum-refseq")

sra_refseq <- read_tsv("inputs/make_metadata/refseq_prokaryotes_metadata_sra.txt", 
                       col_names = c("Run", 'spots', 'bases', 'spots_with_mates',
                                     'avgLength', 'size_MB', 'download_path', 
                                     'Experiment', 'LibraryStrategy', 
                                     'LibrarySelection', 'LibrarySource', 'LibraryLayout',
                                     'InsertSize', 'InsertDev', 'Platform', 'Model', 
                                     'SRAStudy', 'BioProject', 'ProjectID', 'Sample', 
                                     'BioSample', 'SampleType', 'TaxID', 'ScientificName', 
                                     'SampleName', 'CenterName', 'Submission', 'Consent'))

sra_refseq_filt <- sra_refseq %>%
  filter(Platform == "ILLUMINA") %>%
  filter(LibraryLayout == "PAIRED") %>%
  filter(LibraryStrategy == "WGS") %>%
  filter(LibrarySelection == "RANDOM") %>%
  filter(LibrarySource == "GENOMIC") %>%
  group_by(BioSample) %>% 
  filter(n() == 1)

all <- left_join(sra_refseq_filt, prok_refseq, by = c("BioSample" = "bio_sample_accession"))

all_sub <- all %>% 
   group_by(group)  %>%
   sample_n(5, replace = T) %>%
   distinct()
all_sub <- all_sub %>% 
  filter(!is.na(ftp_path))

write_tsv(all_sub, "inputs/metadata.tsv")
