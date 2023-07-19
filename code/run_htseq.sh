#!/bin/bash

# Set the path to the directory containing the BAM files
bam_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/alignment/"

# Set the path to the GFF/GTF annotation file
annotation_file="/dartfs/rc/lab/R/RossB/SalernoP/tetR_RNAseq/ref/bfrag_NCTC_9343.gff"

# Set the path to the directory for the HTSeq-count output files
output_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/htseq_counts/"

# Get the sample name from the BAM file name
sample="$1"
bam_file="${bam_dir}${sample}.sorted.bam"

# Set the path to the HTSeq-count output file
counts_output="${output_dir}${sample}_counts.txt"

# Run HTSeq-counts using the sorted BAM file and the annotation file
htseq-count -f bam -r pos -s no -t gene -i ID "$bam_file" "$annotation_file" > "$counts_output"

# Count completion message
echo "HTSeq-count completed for ${sample}"
