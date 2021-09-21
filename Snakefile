import pandas as pd
metadata = pd.read_csv("inputs/metadata.tsv", sep = "\t", header = 0)
SRR = metadata['Run'].to_list()
ACC = metadata['assembly_accession'].to_list()

#ORPHEUM_DB = ["c__Clostridia", "f__Lachnospiraceae", "p__Firmicutes_A"]
ORPHEUM_DB = ["p__Firmicutes_A"]
# set constrained k sizes
dayhoff_ksizes = [14, 16, 18]
protein_ksizes = [7, 10]
# Snakemake will use the ALPHA_KSIZE wildcard from rule all to generate output file names
# Then, snakemake will back propagate the strings from the final file names to solve for
# the wildcards "alphabet" and "ksize" throughout the rest of the workflow. 
# The underscore for the chrs in the list ALPHA_KSIZE separates the alphabet string from 
# the ksize string, allowing snakemake to solve {alphabet}_{ksize} wildcard strings. 
# Therefore, the chrs in the ALPHA_KSIZE list also set the alphabet names as "dayhoff" and "protein".
ALPHA_KSIZE = expand('protein-k{k}', k=protein_ksizes)
ALPHA_KSIZE += expand('dayhoff-k{k}', k=dayhoff_ksizes)

# Do something similar to constrain SRR:genome accessions
metadata['srr_acc'] = metadata['Run'] + "-" + metadata['assembly_accession']
SRR_ACC =  metadata['srr_acc'].to_list()

rule all:
    input:
        #expand("outputs/orpheum/{orpheum_db}/{alpha_ksize}/{srr}.summary.json",orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE, srr = SRR),
        expand("outputs/aa_paladin/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        expand("outputs/nuc_noncoding_bwa/{orpheum_db}/{alpha_ksize}/multiqc_report.html", orpheum_db = ORPHEUM_DB, alpha_ksize = ALPHA_KSIZE),
        "outputs/assembly_stats/all_assembly_stats.tsv",
        "outputs/assembly_diginorm_bwa/multiqc_report.html"

rule download_sra:
    output: 
        r1="inputs/raw/{srr}_1.fq.gz",
        r2="inputs/raw/{srr}_2.fq.gz"
    conda: "envs/sra-tools.yml"
    params: out_dir= "inputs/raw"
    resources:  mem_mb=lambda wildcards, attempt: attempt *2000
    threads: 1
    shell:"""
    fasterq-dump {wildcards.srr} -O {params.out_dir} -e {threads} -p
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
    resources:  mem_mb=lambda wildcards, attempt: attempt *8000
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
        mem_mb=128000
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 128e9 -V - -o {output}
    '''

rule diginorm_sra:
    input: "outputs/abundtrim/{srr}.abundtrim.fq.gz"
    output: "outputs/diginorm/{srr}.diginorm.fq.gz"
    conda: 'envs/orpheum.yml'
    threads: 1
    resources:
        mem_mb=64000
    shell:'''
    normalize-by-median.py -k 20 -C 20 -M 64e9 --gzip -o {output} {input}
    '''

rule count_reads_in_diginorm_sra:    
    input: "outputs/diginorm/{srr}.diginorm.fq.gz"
    output: "outputs/diginorm_fastp/{srr}.diginorm.fastp.json"
    conda: 'envs/fastp.yml'
    threads: 1
    resources:
        mem_mb=16000
    shell:'''
    fastp -i {input} --interleaved_in -j {output}
    '''

rule orpheum_translate_sra_reads:        
    input: 
        ref="inputs/orpheum_index/{orpheum_db}.{alphabet}-k{ksize}.nodegraph",
        fastq="outputs/diginorm/{srr}.diginorm.fq.gz"
    output:
        pep="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.coding.faa",
        nuc="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.nuc_coding.fna",
        nuc_noncoding="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.nuc_noncoding.fna",
        csv="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.coding_scores.csv",
        json="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum-translate-{srr}-{orpheum_db}-{alphabet}-k{ksize}.txt"
    resources:  mem_mb=lambda wildcards, attempt: attempt *16000
    threads: 1
    shell:'''
    orpheum translate --alphabet {wildcards.alphabet} --peptide-ksize {wildcards.ksize}  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fastq} > {output.pep}
    '''

##################################################
## PREPROCESS ASSEMBLIES
##################################################

rule download_assemblies:
    output: "inputs/assemblies/{acc}_genomic.fna.gz",
    threads: 1
    resources: mem_mb=1000
    run:
        row = metadata.loc[metadata['assembly_accession'] == wildcards.acc]
        assembly_ftp = row['ftp_path'].values
        assembly_ftp = assembly_ftp[0]
        assembly_ftp = assembly_ftp + "/*genomic.fna.gz"
        shell("wget -O {output} {assembly_ftp}")

rule gunzip_assemblies:
    input: "inputs/assemblies/{acc}_genomic.fna.gz",
    output: "inputs/assemblies/{acc}_genomic.fna",
    threads: 1
    resources: mem_mb=1000
    shell:'''
    gunzip -c {input} > {output}
    '''

rule prodigal_translate_assemblies:
    input: "inputs/assemblies/{acc}_genomic.fna",
    output: 
        logg= "outputs/assembly_prodigal/{acc}.out",
        proteins="outputs/assembly_prodigal/{acc}.proteins.faa",
        genes="outputs/assembly_prodigal/{acc}.genes.fa"
    conda: "envs/prodigal.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
    shell:"""
    prodigal -i {input} -o {output.logg} -a {output.proteins} -d {output.genes} 
    """

rule count_coding_bp_assemblies:
    input: "inputs/assemblies/{acc}_genomic.fna"
    output: "outputs/assembly_stats/{acc}_total_bp.txt"
    threads: 1
    resources: mem_mb=1000
    shell:'''
    grep -v ">" {input} | wc -m > {output}
    '''

rule count_total_aa_assemblies:
    input: "outputs/assembly_prodigal/{acc}.proteins.faa"
    output: "outputs/assembly_stats/{acc}_aa.txt"
    threads: 1
    resources: mem_mb=1000
    shell:'''
    grep -v ">" {input} | wc -m > {output}
    '''

rule summarize_coding_assemblies:
    input:
        aa=expand("outputs/assembly_stats/{acc}_aa.txt", acc = ACC),
        bp=expand("outputs/assembly_stats/{acc}_total_bp.txt", acc = ACC)
    output: tsv = "outputs/assembly_stats/all_assembly_stats.tsv"
    threads: 1
    resources: mem_mb = 4000
    conda: "envs/tidyverse.yml"
    script: "scripts/combine_assembly_stats.R"

rule index_nuc_cds_assemblies:
    input: "outputs/assembly_prodigal/{acc}.genes.fa"
    output: "outputs/assembly_prodigal/{acc}.genes.fa.bwt"
    conda: "envs/bwa.yml"
    resources: mem_mb = 2000
    threads: 1
    shell:'''
    bwa index {input}
    ''' 

rule map_nucleotide_reads_against_nucleotide_cds:
    input: 
        ref_nuc_cds= "outputs/assembly_prodigal/{acc}.genes.fa",
        ref_nuc_cds_bwt= "outputs/assembly_prodigal/{acc}.genes.fa.bwt",
        reads="outputs/diginorm/{srr}.diginorm.fq.gz"
    output: temp("outputs/assembly_diginorm_bwa/{srr}-{acc}.bam")
    conda: "envs/bwa.yml"
    threads: 4
    resources: mem_mb = 4000
    shell:'''
    bwa mem -p -t {threads} {input.ref_nuc_cds} {input.reads} | samtools sort -o {output} -
    '''

rule flagstat_map_nucleotide_reads_against_nucleotide_cds:
    input: "outputs/assembly_diginorm_bwa/{srr}-{acc}.bam"
    output: "outputs/assembly_diginorm_bwa/{srr}-{acc}.flagstat"
    conda: "envs/bwa.yml"
    resources: mem_mb = 2000
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule multiqc_flagstat_map_nucleotide_reads_against_nucleotide_cds:
    input: expand("outputs/assembly_diginorm_bwa/{srr_acc}.flagstat", srr_acc = SRR_ACC)
    output: "outputs/assembly_diginorm_bwa/multiqc_report.html"
    params: 
        iodir = "outputs/assembly_diginorm_bwa/"
    conda: "envs/multiqc.yml"
    resources: mem_mb = 8000
    threads: 1
    shell:'''
    multiqc {params.iodir} -o {params.iodir} 
    '''

####################################################
## Evaluate orpheum
####################################################

## Amino Acid ######################################

rule paladin_index_assembly_aa:
    input: "outputs/assembly_prodigal/{acc}.proteins.faa" 
    output:"outputs/assembly_prodigal/{acc}.proteins.faa.pro",
    conda: "envs/paladin.yml"
    resources: mem_mb = 8000
    threads: 1
    shell:'''
    paladin index -r3 {input}
    '''

rule paladin_align_aa:
    input: 
        ref="outputs/assembly_prodigal/{acc}.proteins.faa",
        idx="outputs/assembly_prodigal/{acc}.proteins.faa.pro",
        pep="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.coding.faa"
    output: temp("outputs/aa_paladin/{orpheum_db}/{alphabet}-k{ksize}/{srr}-{acc}.aa.sam")
    conda: "envs/paladin.yml"
    resources: mem_mb = 6000
    threads: 1
    shell:'''
    paladin align -t 1 -p {input.ref} {input.pep} > {output}
    '''

rule samtools_flagstat_paladin_aa:
    input: "outputs/aa_paladin/{orpheum_db}/{alphabet}-k{ksize}/{srr}-{acc}.aa.sam"
    output:"outputs/aa_paladin/{orpheum_db}/{alphabet}-k{ksize}/{srr}-{acc}.aa.flagstat"
    conda: "envs/paladin.yml"
    resources: mem_mb = 2000
    threads: 1
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule multiqc_samtools_flagstat_paladin:
    input: expand("outputs/aa_paladin/{{orpheum_db}}/{{alphabet}}-k{{ksize}}/{srr_acc}.aa.flagstat", srr_acc = SRR_ACC)
    output: "outputs/aa_paladin/{orpheum_db}/{alphabet}-k{ksize}/multiqc_report.html"
    resources: mem_mb = 8000
    threads: 1
    params: 
        iodir = lambda wildcards: "outputs/aa_paladin/" + wildcards.orpheum_db + "/" + wildcards.alphabet + "-k" + wildcards.ksize,
    conda: "envs/multiqc.yml"
    resources: mem_mb = 8000
    threads: 1
    shell:'''
    multiqc {params.iodir} -o {params.iodir} 
    '''

## Nucleotide noncoding ##########################

rule cut_dedup_nuc_noncoding_read_names:
    input: "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.nuc_noncoding.fna",
    output:  "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.nuc_noncoding.cut.dedup.fna.gz",
    resources: mem_mb = 8000
    threads: 1
    shell:'''
    sed '/^>/ s/__.*//g' {input} | awk '/^>/{{f=!d[$1];d[$1]=1}}f' | gzip > {output}
    '''

rule grab_cut_dedup_coding_read_names:
    input: "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.coding.faa", 
    output: "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.aa_names.cut.dedup.txt",  
    resources: mem_mb = 8000
    threads: 1
    shell:'''
    grep ">" {input} | sed '/^>/ s/__.*//g' | awk '/^>/{{f=!d[$1];d[$1]=1}}f' > {output}
    '''

rule isolate_noncoding_only_reads:
    input: 
        pep ="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.aa_names.cut.dedup.txt",  
        noncoding = "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.nuc_noncoding.cut.dedup.fna.gz",
    output: "outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.nuc_noncoding.cut.dedup.only.fna.gz",
    #conda: "envs/bioawk.yml"
    resources: mem_mb = 8000
    threads: 1
    shell:"""
    cat {input.noncoding} | paste - - | grep -v -F -f {input.pep} | tr "\t" "\n" | gzip > {output}
    """

rule map_nuc_noncoding_to_ref_nuc_cds: 
    input: 
        ref_nuc_cds= "outputs/assembly_prodigal/{acc}.genes.fa",
        ref_nuc_cds_bwt= "outputs/assembly_prodigal/{acc}.genes.fa.bwt",
        nuc_noncoding="outputs/orpheum/{orpheum_db}/{alphabet}-k{ksize}/{srr}.nuc_noncoding.cut.dedup.only.fna.gz",
    output: temp("outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}-k{ksize}/{srr}-{acc}.nuc_noncoding.bam")
    conda: "envs/bwa.yml"
    resources: mem_mb = 8000
    threads: 1
    shell:'''
    bwa mem -t {threads} {input.ref_nuc_cds} {input.nuc_noncoding} | samtools sort -o {output} -
    '''

rule flagstat_map_nuc_noncoding_to_ref_nuc_cds:
    input: "outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}-k{ksize}/{srr}-{acc}.nuc_noncoding.bam"
    output: "outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}-k{ksize}/{srr}-{acc}.nuc_noncoding.flagstat"
    conda: "envs/bwa.yml"
    resources: mem_mb = 8000
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule multiqc_flagstat_map_nuc_noncoding_to_ref_nuc_set:
    input: expand("outputs/nuc_noncoding_bwa/{{orpheum_db}}/{{alphabet}}-k{{ksize}}/{srr_acc}.nuc_noncoding.flagstat", srr_acc = SRR_ACC), 
    output: "outputs/nuc_noncoding_bwa/{orpheum_db}/{alphabet}-k{ksize}/multiqc_report.html"
    params: 
        iodir = lambda wildcards: "outputs/nuc_noncoding_bwa/" + wildcards.orpheum_db + "/" + wildcards.alphabet + "-k" + wildcards.ksize,
    conda: "envs/multiqc.yml"
    resources: mem_mb = 8000
    threads: 1
    shell:'''
    multiqc {params.iodir} -o {params.iodir} 
    '''
