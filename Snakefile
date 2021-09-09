import pandas as pd
metadata = pd.read_csv("inputs/metadata.tsv", sep = "\t", header = 0)
SRR = metadata['Run']
#ORPHEUM_DB = ["c__Clostridia", "f__Lachnospiraceae", "p__Firmicutes_A"]
ORPHEUM_DB = ["p__Firmicutes_A"]
# set constrained k sizes
dayhoff_ksizes = [14, 16, 18]
protein_ksizes = [7, 10, 11]
# Snakemake will use the ALPHA_KSIZE wildcard from rule all to generate output file names
# Then, snakemake will back propagate the strings from the final file names to solve for
# the wildcards "alphabet" and "ksize" throughout the rest of the workflow. 
# The underscore for the chrs in the list ALPHA_KSIZE separates the alphabet string from 
# the ksize string, allowing snakemake to solve {alphabet}_{ksize} wildcard strings. 
# Therefore, the chrs in the ALPHA_KSIZE list also set the alphabet names as "dayhoff" and "protein".
ALPHA_KSIZE = expand('protein-k{k}', k=protein_ksizes)
ALPHA_KSIZE += expand('dayhoff-k{k}', k=dayhoff_ksizes)

rule all:
    input:
        #expand("outputs/aa_paladin/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        #expand("outputs/nuc_noncoding_bwa/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        #expand("outputs/nuc_coding_bwa/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        #expand("outputs/nuc_noncoding_bwa/{orpheum_db}/{alpha_ksize}/{srr}.nuc_noncoding.stat", alpha_ksize=ALPHA_KSIZE, orpheum_db = ORPHEUM_DB, srr = SRR),
        #expand("outputs/nuc_coding_bwa/{orpheum_db}/{alpha_ksize}/{srr}.nuc_coding.stat", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE, srr = SRR),
        expand("outputs/orpheum/{orpheum_db}/{alpha_ksize}/{srr}.summary.json",orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE, srr = SRR),


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

rule kmertrim_sra:
    input: 
        r1 = 'outputs/fastp/{srr}_R1.fastp.fq.gz',
        r2 = 'outputs/fastp/{srr}_R2.fastp.fq.gz',
    output: "outputs/abundtrim/{srr}.abundtrim.fq.gz"
    conda: 'envs/orpheum.yml'
    threads: 1
    resources:
        mem_mb=64000
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

rule orpheum_translate_sgc_nbhds:        
    input: 
        ref="inputs/orpheum_index/{orpheum_db}.{alphabet}-k{ksize}.nodegraph",
        fastq="outputs/abundtrim/{srr}.abundtrim.fq.gz"
    output:
        pep="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.coding.faa",
        nuc="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.nuc_coding.fna",
        nuc_noncoding="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.reads.nuc_noncoding.fna",
        csv="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.coding_scores.csv",
        json="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum-translate-{srr}-{orpheum_db}-{alphabet}-k{ksize}.txt"
    resources: mem_mb = 16000
    threads: 1
    shell:'''
    orpheum translate --alphabet {wildcards.alphabet} --peptide-ksize {wildcards.ksize}  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fastq} > {output.pep}
    '''

