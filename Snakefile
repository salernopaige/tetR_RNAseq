# RNA-Seq analysis for B. fragilis samples treated with and without tetracycline
# Paige Salerno July 13th, 2023

from pathlib import Path

# Load configuration from YAML file
configfile: "config.yaml"
fastq=Path(config['fastq_dir'])
trimmed=Path(config['trimmed_dir'])
alignment=Path(config['alignment_dir'])
htseq=Path(config['htseq_counts_dir'])
sample=config['samples']

# Target rule
rule all:
    input:
        "%s/all_counts.txt" %(htseq),
        expand("%s/{sample}_R1.fastp.fastq.gz" %(trimmed), sample=sample),
        expand("%s/{sample}_R2.fastp.fastq.gz" %(trimmed), sample=sample),
        expand("%s/{sample}.sorted.bam" %(alignment), sample=sample),
        expand("%s/{sample}.sorted.bam.bai" %(alignment), sample=sample),
        expand("%s/{sample}_counts.txt" %(htseq), sample=sample)

# Rule to run fastp
rule run_fastp:
    input:
        S1="%s/{sample}_R1.fastq.gz" %(fastq),
        S2="%s/{sample}_R2.fastq.gz" %(fastq)
    output:
        out_read1="%s/{sample}_R1.fastp.fastq.gz" %(trimmed),
        out_read2="%s/{sample}_R2.fastp.fastq.gz" %(trimmed)
    shell:
        "./code/run_fastp.sh {wildcards.sample}"

rule fastp_all:
     input:
        read1=expand("%s/{sample}_R1.fastp.fastq.gz" %(trimmed), sample=sample),
        read2=expand("%s/{sample}_R2.fastp.fastq.gz" %(trimmed), sample=sample)

# Rule to run BWA alignment
rule run_bwa:
    input:
        read1="%s/{sample}_R1.fastp.fastq.gz" %(trimmed),
        read2="%s/{sample}_R2.fastp.fastq.gz" %(trimmed)
    output:
        bam_output="%s/{sample}.sorted.bam" %(alignment),
        bam_index="%s/{sample}.sorted.bam.bai" %(alignment)
    shell:
        "./code/run_bwa.sh {wildcards.sample}"

rule bwa_all:
    input:
        bam_output=expand("%s/{sample}.sorted.bam" %(alignment), sample=sample),
        bam_index=expand("%s/{sample}.sorted.bam.bai" %(alignment), sample=sample)

# Rule to run HTSeq-count
rule run_htseq:
    input:
        bam_file="%s/{sample}.sorted.bam" %(alignment)
    output:
        counts_output="%s/{sample}_counts.txt" %(htseq)
    shell:
        "./code/run_htseq.sh {wildcards.sample}"

rule htseq_all:
    input:
        expand("%s/{sample}_counts.txt" %(htseq), sample=sample)

# Rule to run the expression matrix script
rule run_exp_matrix:
    input:
        counts_files=expand("%s/{sample}_counts.txt" %(htseq), sample=sample)
    output:
        expression_matrix="%s/all_counts.txt" %(htseq)
    shell:
        "./code/run_exp_matrix.sh {input.counts_files}"

