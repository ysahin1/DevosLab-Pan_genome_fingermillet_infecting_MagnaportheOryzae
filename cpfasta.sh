# Define source and destination directories
SOURCE_DIR="Path_to_BAM"
DEST_DIR="Path_to_FASTAs"

# Loop through the files in the source directory
for file in "$SOURCE_DIR"/*E*_trimed.fastq; 
do
    # Extract the base name (e.g., T9) and the part number (e.g., 1 or 2)
    base_name=$(basename "$file" | cut -d'_' -f1)
    part_number=$(basename "$file" | cut -d'_' -f2)

    # Construct the new file name
    new_file="${base_name}.R${part_number}.fasta"

    # Copy and rename the file to the destination directory
    seqtk seq -a "$file" > "$DEST_DIR/$new_file"
done

