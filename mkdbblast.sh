SOURCE_DIR="Path_to_Fastas"

for file in $SOURCE_DIR/*.fasta;
do
    # Extract the base name (e.g., T9) and the part number (e.g., 1 or 2)
    base_name=$(basename "$file" | cut -d'.' -f1)
    part_number=$(basename "$file" | cut -d'.' -f2)

    # Construct the new file name
    new_file="${base_name}.${part_number}.fasta"

    # Copy and rename the file to the destination directory
    mkdir ../Straindb/"${base_name}_${part_number}"
    makeblastdb -in ../Fastas/"$new_file" -out ../Straindb/"${base_name}_${part_number}"/"${base_name}.${part_number}" -dbtype 'nucl'
    blastn -query ../pot2seq/transposon.fa -task blastn-short -db ../Straindb/"${base_name}_${part_number}"/"${base_name}.${part_number}" -reward 2 -penalty -3 -word_size 4 -gapopen 6 -gapextend 2 -num_alignments 20000 -out ../Pot2LocInStrains/r2p-3g6e2_O6_${part_number}_"${base_name}".txt -outfmt 6 -num_threads 10
done
