#!/bin/bash

# Set the path to the directory containing the BAM files
bam_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/alignment/"

# Set the path to the GFF/GTF annotation file
annotation_file="/dartfs/rc/lab/R/RossB/SalernoP/tetR_RNAseq/ref/bfrag_NCTC_9343.gff"

# Set the path to the directory for the HTSeq-count output files
output_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/htseq_counts/"

# Loop through all the BAM files in the directory
for bam_file in "$bam_dir"*_aligned.bam; do
    # Get the sample name from the BAM file name
    sample=$(basename "$bam_file" _aligned.bam)

    # Set the path to the HTSeq-count output file
    counts_output="${output_dir}${sample}_counts.txt"

    # Run HTSeq-counts using the sorted BAM file and the annotation file
    htseq-count -f bam -r pos -s no -t gene -i ID "$sorted_bam" "$annotation_file" > "$counts_output"

    # Count completion message
    echo "HTSeq-count completed for ${sample}"

done
