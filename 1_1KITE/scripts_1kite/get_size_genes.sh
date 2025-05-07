#!/bin/bash

#MUST ADJUST THESE PARAMETERS EVERYTIME
path_to_contig="/scratch/group/songlab/phylo/OR_TE/assembly/read_data_plate16_20/contig"
path_to_ortho_results="/scratch/group/songlab/phylo/OR_TE/environment/software/Orthograph/orthograph_results/orthograph_Taeniopoda"

#MUST ADJUST EVERYTIME = TYPE "y" FOR YES, TYPE "n" FOR NO
multiple_directories="n"
  

# STEP 1: Get the size of each contig

cd "$path_to_contig"

ls -lh . | awk '{print $9, $5}' > "$path_to_ortho_results/contig_size.tsv"

echo "Part one done"

# STEP 2a: Process multiple subdirectories

if [ "$multiple_directories" == "y" ]; then
    array=(1 2 3 4 5 6 7 8 9)

    for num in "${array[@]}" ; do
        cd "$path_to_ortho_results/orthograph_JBL_Ensifera_${num}R"
        for dir in */ ; do
            echo "$dir    $(ls "$dir/aa" | wc -l)" >> "$path_to_ortho_results/genes_recovered.tsv"
        done
        cd ..
    done

    echo "Part two done"

else

# STEP 2b: Process a single subdirectory

cd "$path_to_ortho_results/"

for dir in */ ; do
    echo "$dir    $(ls "$dir/aa" | wc -l)" >> "$path_to_ortho_results/genes_recovered.tsv"
done

echo "Part two done"

fi

