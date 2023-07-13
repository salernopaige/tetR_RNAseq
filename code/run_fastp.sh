#!/bin/bash

# Set the path to the fastq files
indir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/raw_reads/"

# Set the path to the output files
outdir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/trimmed/"

# Loop through all the input FASTQ files in the directory
for fastq in "$indir"*_001.fastq.gz; do
    # Get the sample name from the FASTQ file name
    sample=$(basename "$fastq" _R1_001.fastq.gz)

    # Set the paths to the input fastq files
    read1="${indir}${sample}_R1.fastq.gz"
    read2="${indir}${sample}_R2.fastq.gz"

    # Set the paths to the output FASTQ files after trimming
    out_read1="${outdir}${sample}_R1.fastp.fastq.gz"
    out_read2="${outdir}${sample}_R2.fastp.fastq.gz"

    # Run fastq to trim the reads
    fastp -i "$read1" -I "$read2" -o "$out_read1" -O "$out_read2"

    echo "Succesfully trimmed ${sample}"
    date

done