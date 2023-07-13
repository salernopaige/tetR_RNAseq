#!/bin/bash

# Set the path to the BWA index file
index="/dartfs/rc/lab/R/RossB/SalernoP/tetR_RNAseq/ref/bfrag"

# Set the path to the directory containing the FASTQ files
fastq_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/raw_reads/"

# Set the path to the directory for the output alignment files
alignment_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/alignment/"

# Loop through all the FASTQ files in the directory
for fastq_file in "$fastq_dir"*.fastp.fastq.gz; do
    # Get the sample name from the FASTQ file name
    sample=$(basename "$fastq_file" _R1_fastp_fastq.gz)
  
    # Set the paths to the input FASTQ files
    read1="${fastq_dir}${sample}_R1.fastp.fastq.gz"
    read2="${fastq_dir}${sample}_R2.fastp.fastq.gz"
  
    # Set the path to the output alignment file
    alignment_output="${alignment_dir}${sample}_aligned.sam"
  
    # Run BWA alignment using bwa mem
    bwa mem "$index" "$read1" "$read2" > "$alignment_output"

    # Completed alignment message
    echo "Alignment completed for ${sample}"

    # Sort SAM file using samtools 
    sort_file="${sample}_aligned.sorted.sam"
    samtools sort "$alignment_output" -o "$sort_file"

    # Convert SAM to BAM using samtools
    bam_output="${alignment_dir}${sample}_aligned.bam"
    samtools view -bS "$alignment_output" > "$bam_output"

    # Index the BAM file using samtools
    samtools index "$bam_output"

    # BAM processing message
    echo "Completed Samtools for ${sample}"

done
