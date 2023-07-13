#!/bin/bash

# Set the path to the directory for the HTSeq-count output files
output_dir="/dartfs/rc/lab/R/RossB/RobitailleS/rnaseq_230705/htseq_counts/"

# Set up an array that we will fill with shorthand sample names
myarray=()

# Loop over htseq-counts files and extract the 2nd column (the raw read counts) using the 'cut' command
for counts_file in "$output_dir"*.txt; do
    # Get the sample name from the counts file name
    sample=$(basename "$counts_file" _counts.txt)

    # Split up sample names to remove everything after "-"
    sname=$(echo "$sample" | cut -d "-" -f 1)

    # Extract the second column of the file to get read counts only
    echo "Counts for $sname being extracted"
    cut -f 2 "$counts_file" > "${sname}.tmp.counts"

    # Save shorthand sample names into an array
    myarray+=("$sname")

    # Extraction message
    echo "Read counts extracted for ${sample}"
done

# Extract ENSG gene IDs from one of the files
cut -f 1 "${output_dir}3492_1_S1_counts.txt" > gene_IDs.txt

# Use the paste command to combine gene IDs and raw counts for all files into one file
paste gene_IDs.txt "${output_dir}"*.tmp.counts > tmp_all_counts.txt

# Look at the contents of the array we made with shorthand sample names
echo ${myarray[@]}

# Print contents of array into text file with each element on a new line
printf "%s\n" "${myarray[@]}" > names.txt
cat names.txt

# Create a file to fill
touch all_counts.txt

# Use the 'cat' command to combine the tmp.counts.txt files into all_counts.txt
cat <(printf "%s\n" "${myarray[@]}" | sort | paste -s) tmp_all_counts.txt > all_counts.txt

# Remove the tmp files
rm -f *tmp*

# Count file message
echo "HTSeq-counts complete"
