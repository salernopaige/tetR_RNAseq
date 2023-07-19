#!/bin/bash

# Set the path to the BWA index file
index="/dartfs/rc/lab/R/RossB/SalernoP/tetR_RNAseq/ref/bfrag"

# Set the path to the directory containing the FASTQ files
trimmed_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/trimmed/"

# Set the path to the directory for the output alignment files
alignment_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/alignment/"

# Get the sample name from the FASTQ file name
sample="$1"

# Set the paths to the input FASTQ files
read1="${trimmed_dir}${sample}_R1.fastp.fastq.gz"
read2="${trimmed_dir}${sample}_R2.fastp.fastq.gz"

# Run BWA alignment using bwa mem
align_out="${alignment_dir}${sample}_aligned.sam"
bwa mem "$index" "$read1" "$read2" > "$align_out"

# Sort SAM file using samtools
sort_file="${alignment_dir}${sample}.sorted.sam"
samtools sort "$align_out" -o "$sort_file"

# Sort and convert SAM to BAM using samtools
alignment_output="${alignment_dir}${sample}.sorted.bam"
samtools view -Sb "$sort_file" > "$alignment_output"

# Index the BAM file using samtools
samtools index "$alignment_output"

# Remove the temporary SAM file
rm "$align_out"
rm "$sort_file"

# Completed alignment message
echo "Alignment completed for ${sample}"
