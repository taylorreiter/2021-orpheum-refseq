library(dplyr)
library(purrr)
library(readr)

total_bp <- unlist(snakemake@input[["bp"]]) %>%
    set_names() %>%
    map_dfr(read_tsv, col_names = "total_bp", .id = "assembly") %>%
    mutate(assembly = gsub("_total_bp.txt", "", assembly))

coding_aa <- unlist(snakemake@input[["aa"]]) %>%
    set_names() %>%
    map_dfr(read_tsv, col_names = "cds_aa", .id = "assembly") %>%
    mutate(assembly = gsub("_aa.txt", "", assembly)) %>%
    mutate(cds_bp = cds_aa * 3)

assembly_stats <- left_join(total_bp, coding_aa, by = "assembly")
write_tsv(assembly_stats, snakemake@output[['tsv']]) 
