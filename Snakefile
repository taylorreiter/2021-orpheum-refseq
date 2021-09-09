import pandas as pd
metadata = pd.read_csv("inputs/metadata.tsv", sep = "\t", header = 0)
SRR = metadata['Run']

ORPHEUM_DB = ["roary_with_megahit_and_isolates", "ruminococcusB", "f__Lachnospiraceae", "p__Firmicutes_A"]
# set constrained k sizes
dayhoff_ksizes = [11, 13, 15, 17]
protein_ksizes = [7, 10, 11]
# Snakemake will use the ALPHA_KSIZE wildcard from rule all to generate output file names
# Then, snakemake will back propagate the strings from the final file names to solve for
# the wildcards "alphabet" and "ksize" throughout the rest of the workflow. 
# The underscore for the chrs in the list ALPHA_KSIZE separates the alphabet string from 
# the ksize string, allowing snakemake to solve {alphabet}_{ksize} wildcard strings. 
# Therefore, the chrs in the ALPHA_KSIZE list also set the alphabet names as "dayhoff" and "protein".
ALPHA_KSIZE = expand('protein_ksize{k}', k=protein_ksizes)
ALPHA_KSIZE += expand('dayhoff_ksize{k}', k=dayhoff_ksizes)

rule all:
    input:
        expand("outputs/aa_paladin/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        expand("outputs/nuc_noncoding_bwa/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        expand("outputs/nuc_coding_bwa/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        expand("outputs/nuc_noncoding_bwa/{orpheum_db}/{alpha_ksize}/{}.nuc_noncoding.stat", alpha_ksize=ALPHA_KSIZE, orpheum_db = ORPHEUM_DB, ...),
        expand("outputs/nuc_coding_bwa/{orpheum_db}/{alpha_ksize}/{}.nuc_coding.stat", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE, ...),

rule download_sra:
    output: 
        r1="inputs/raw/{srr}_1.fq.gz",
        r2="inputs/{srr}_2.fq.gz"
    conda: "envs/sra-tools.yml"
    params: out_dir= "inputs/raw"
    resources: mem_mb = 2000
    threads: 1
    shell:"""
    fasterq-dump {wilcards.srr} -O {params.out_dir} -e {threads} -p
    gzip -9c {params.out_dir}/{wildcards.srr}_1.fastq > {output.r1}
    gzip -9c {params.out_dir}/{wildcards.srr}_2.fastq > {output.r2}
    rm {params.out_dir}/{wildcards.srr}_*.fastq
    """

rule fastp_sra:
    input: 
        r1 = "inputs/raw/{srr}_1.fq.gz",
        r2 = "inputs/raw/{srr}_2.fq.gz"
    output: 
        r1 = 'outputs/fastp/{srr}_R1.fastp.fq.gz',
        r2 = 'outputs/fastp/{srr}_R2.fastp.fq.gz',
        json = 'outputs/fastp/{srr}.json'
    conda: 'envs/fastp.yml'
    threads: 1
    resources: mem_mb = 2000
    shell:'''
    fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -q 4 -j {output.json} -l 31 -c
    '''

rule orpheum_translate_sgc_nbhds:        
    input: 
        ref="outputs/orpheum_index/{orpheum_db}_{alphabet}_ksize{ksize}.bloomfilter.nodegraph",
        fastq="outputs/{library}.fastq"
    output:
        pep="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}.coding.faa",
        nuc="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}.nuc_coding.fna",
        nuc_noncoding="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}.reads.nuc_noncoding.fna",
        csv="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}.coding_scores.csv",
        json="outputs/orpheum/{orpheum_db}/{alphabet}_ksize{ksize}/{library}.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum_translate_{library}_{orpheum_db}_{alphabet}_ksize{ksize}.txt"
    resources: mem_mb = 16000
    threads: 1
    shell:'''
    orpheum translate --alphabet {wildcards.alphabet} --peptide-ksize {wildcards.ksize}  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fastq} > {output.pep}
    '''

