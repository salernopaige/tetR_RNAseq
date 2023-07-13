# RNA-Seq analysis for B. fragilis samples treated with and without tetracycline
# Paige Salerno July 13th, 2023

# Load configuration from YAML file
configfile: "config.yaml"

# Sample names
samples = config['samples']

# Target rule
rule all:
    input:
        "all_counts.txt",
	"{sample}_fastp_R1.fastq.gz",


# Makes links to the base library so that they are named consistently with the folder names
# Snakemake works really well if we have files that go {sample}/{sample}.fastq.gz
rule make_read_links:
    input:
        dir=ancient("%s/" %(fastq_dir))
    output:
        R1="%s/{sample}_R1.fastq.gz" %(fastq_dir),
        R2="%s/{sample}_R2.fastq.gz" %(fastq_dir)
    run:
        infiles_R1 = list(Path(input.dir).glob("*R1_001.fastq.gz"))
        infiles_R2 = list(Path(input.dir).glob("*R2_001.fastq.gz"))
        shell("ln -s %s {output.R1}" %(infiles_R1[0]))
        shell("ln -s %s {output.R2}" %(infiles_R2[0]))

# Rule to run fastp
rule run_fastp:
    input:
        S1=config['fastq_dir'] + "{sample}_R1_001.fastq.gz",
        S2=config['fastq_dir'] + "{sample}_R2_001.fastq.gz"
    output:
        out_read1=config['trimmed_dir'] + "{sample}_fastp_R1.fastq.gz",
        out_read2=config['trimmed_dir'] + "{sample}_fastp_R2.fastq.gz"
    shell:
        "./run_fastp.sh"

# Rule to run BWA alignment
rule run_bwa:
    input:
        read1=config['trimmed_dir'] + "{sample}_fastp_R1.fastq.gz",
        read2=config['trimmed_dir'] + "{sample}_fastp_R2.fastq.gz"
    output:
        alignment_output=config['alignment_dir'] + "{sample}_aligned.bam",
        bam_output=config['alignment_dir'] + "{sample}_aligned.sorted.bam",
        bam_index=config['alignment_dir'] + "{sample}_aligned.sorted.bam.bai"
    shell:
        "./run_bwa.sh"

# Rule to run HTSeq-count
rule run_htseq:
    input:
        bam_file=config['alignment_dir'] + "{sample}_aligned.sorted.bam"
    output:
        counts_output=config['htseq_counts_dir'] + "{sample}_counts.txt"
    shell:
        "./run_htseq.sh"

# Rule to run the expression matrix script
rule run_exp_matrix:
    input:
        counts_files=expand(config['htseq_counts_dir'] + "{sample}_counts.txt", sample=config['samples'])
    output:
        expression_matrix="all_counts.txt"
    shell:
        "./run_exp_matrix.sh"

