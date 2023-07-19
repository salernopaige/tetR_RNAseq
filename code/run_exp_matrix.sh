#!/bin/bash

# Set the path to the directory for the HTSeq-count output files
output_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/htseq_counts/"

# Get access to snakemake input files
count_files=("$@")

# Set up an array that we will fill with shorthand sample names
myarray=()

for file in "${count_files[@]}"; do
    # Get the sample name from the counts file name
    sample=$(basename "${file}" _counts.txt)

    # Extract the second column of the file to get read counts only
    echo "Counts for ${sample} being extracted"
    cut -f 2 "${file}" > "${output_dir}${sample}.tmp.counts.txt"

    # Save shorthand sample names into an array
    myarray+=("${sample}")

    # Extraction message
    echo "Read counts extracted for ${sample}"

done

# Extract ENSG gene IDs from one of the files
cut -f 1 "${output_dir}tetR1_counts.txt" > "${output_dir}gene_IDs.txt"

# Use the paste command to combine gene IDs and raw counts for all files into one file
paste_cmd="paste ${output_dir}gene_IDs.txt "
for sample in "${myarray[@]}"; do
    paste_cmd+="${output_dir}${sample}.tmp.counts.txt "
done
paste_cmd+="> ${output_dir}tmp_all_counts.txt"

# Execute paste commant 
eval "$paste_cmd"

# Look at the contents of the array we made with shorthand sample names
echo ${myarray[@]}

# Print contents of array into text file with each element on a new line
printf "%s\n" "${myarray[@]}" > "${output_dir}names.txt"
cat "${output_dir}names.txt"

# Create a file to fill
touch "${output_dir}all_counts.txt"

# Use the 'cat' command to combine the tmp.counts.txt files into all_counts.txt
cat <(printf "%s\n" "${myarray[@]}" | sort | paste -s) "${output_dir}tmp_all_counts.txt" > "${output_dir}all_counts.txt"

# Remove the tmp files
rm -f "${output_dir}"*tmp*

# Count file message
echo "HTSeq-counts complete"
