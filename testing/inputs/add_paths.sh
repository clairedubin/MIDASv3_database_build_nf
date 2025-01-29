#!/bin/bash

# Define the input and output file paths
input_file="genomes.tsv"
output_file="genomes_with_paths.tsv"
fasta_dir="genomes"

# Check if the input file exists
if [[ ! -f "$input_file" ]]; then
  echo "Error: Input file '$input_file' not found."
  exit 1
fi

# Get the full path of the fasta directory
fasta_dir_full=$(realpath "$fasta_dir")

# Create the output file and add the header
header=$(head -n 1 "$input_file")
echo -e "$header\tfasta_path" > "$output_file"

# Read through the input file, skipping the header
while IFS=$'\t' read -r genome species representative genome_is_representative; do
  # Skip the header row
  if [[ "$genome" == "genome" ]]; then
    continue
  fi

  # Construct the expected fasta file path
  fasta_path="$fasta_dir_full/$genome.fasta"

  # Add the row to the output file, including the fasta path
  echo -e "$genome\t$species\t$representative\t$genome_is_representative\t$fasta_path" >> "$output_file"
done < "$input_file"

# Notify the user of completion
echo "The updated file with fasta paths has been saved to '$output_file'."

