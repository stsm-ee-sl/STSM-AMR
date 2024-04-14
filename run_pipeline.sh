#!/bin/bash

CURRENT_DIR=$(pwd)
FASTA_DIR="${CURRENT_DIR}/analysis/spades/scaffolds_all"

######## AMRFinderPlus ##########################################################################################

source /home/kasutaja/miniconda3/etc/profile.d/conda.sh
conda_envs='/home/kasutaja/miniconda3/envs'
conda activate $conda_envs/amrfinder

#KMERFINDER_DIR="${CURRENT_DIR}/analysis/kmerfinder"
#KMERFINDER_TXT="${KMERFINDER_DIR}/kmerfinder_species_summary.txt"
KMERFINDER_CSV="${CURRENT_DIR}/analysis/kmerfinder/kmerfinder_species_summary.tsv"
#sed 's/,/\t/g' "$KMERFINDER_TXT" > "$KMERFINDER_TSV"


OUTPUT_DIR="${CURRENT_DIR}/analysis/amrfinder" && mkdir -p "$OUTPUT_DIR"


declare -A species_map=(
    ["Campylobacter jejuni"]="Campylobacter"
    ["Campylobacter coli"]="Campylobacter"
    ["Campylobacter hepaticus"]="Campylobacter"
    ["Escherichia coli"]="Escherichia"
    ["Salmonella enterica"]="Salmonella"
)

# Iterate over the kmerfinder CSV, processing each line
tail -n +2 "$KMERFINDER_CSV" | while IFS=$'\t' read -r filename species ncbi_assembly query_cov p_value; do
    #echo "Debug: query coverage for $filename is $query_cov"
    echo "Processing species '$species' for sample '$filename'"

    # Check if the species name exists in the map
    if [[ -z "${species_map["$species"]}" ]]; then
        echo "Species name '$species' not found in mapping. Skipping $filename."
        continue
    fi

    # Extract the AMR Finder species name from the map
    amr_species="${species_map["$species"]}"

    # Ensure query coverage is numeric and greater than 70
    if [[ "$query_cov" =~ ^[+-]?[0-9]+\.?[0-9]*$ ]] && (( $(echo "$query_cov > 70" | bc -l) )); then
        fasta_file="${FASTA_DIR}/${filename}_SPAdes_scaffolds.fasta"
        amrfinder -n "$fasta_file" -O "$amr_species" -o "$OUTPUT_DIR/${filename}_amrfinder_results.tsv" --plus
    #else
        #echo "Query coverage is not numeric or below threshold for $filename. Skipping."
    fi
done


#cd "$OUTPUT_DIR" && for file in *.txt; do sed 's/\t/,/g' "$file" > "${file%.txt}.csv"; done 
#&& find . -type f -name "*.txt" -delete

conda deactivate




######## CARD-RGI ################################################################################

conda activate $conda_envs/rgi

RGI_OUTPUT_DIR="${CURRENT_DIR}/analysis/rgi" && mkdir -p "$RGI_OUTPUT_DIR"

cd "$FASTA_DIR" && for i in *.fasta; do
    rgi main --input_sequence "$i" --output_file "${RGI_OUTPUT_DIR}/${i/.fasta/}_rgi_results.tsv" --clean
done



cd "$RGI_OUTPUT_DIR" && for file in *.txt; do sed 's/\t/,/g' "$file" > "${file%.txt}.tsv"; done && find . -type f -name "*.txt" -delete

conda deactivate
echo "RGI FINISHED"



######## ResFinder ###############################################################################

conda activate $conda_envs/resfinder

OUTPUT_DIR_resfinder="${CURRENT_DIR}/analysis/resfinder" && mkdir -p "$OUTPUT_DIR_resfinder"

#cd "$FASTA_DIR" && for i in *.fasta; do
 #   /home/tom/mambaforge/envs/resfinder/bin/run_resfinder.py -ifa "$i" -o "${OUTPUT_DIR_resfinder}/${i/.fasta/}" -s "campylobacter jejuni" -l 0.6 -t 0.8 --acquired --point
#done


# Your initial script parts remain unchanged

declare -A resfinder_map=(
    ["Campylobacter jejuni"]="campylobacter jejuni"
    ["Campylobacter coli"]="campylobacter coli"
    ["Campylobacter hepaticus"]="campylobacter"
    ["Escherichia coli"]="escherichia coli"
    ["Salmonella enterica"]="salmonella enterica"
)

# Iterate over the kmerfinder CSV, processing each line
tail -n +2 "$KMERFINDER_CSV" | while IFS=$'\t' read -r filename species ncbi_assembly query_cov p_value; do
    echo "Processing species '$species' for sample '$filename'"

    # Corrected the check for species name existence in the map
    if [[ -z "${resfinder_map["$species"]}" ]]; then
        echo "Species name '$species' not found in mapping. Skipping $filename."
        continue
    fi

    # The rest of your script appears correctly structured, assuming all variables are defined properly elsewhere
    resfinder_species="${resfinder_map["$species"]}"

    if [[ "$query_cov" =~ ^[+-]?[0-9]+\.?[0-9]*$ ]] && (( $(echo "$query_cov > 70" | bc -l) )); then
        fasta_file="${FASTA_DIR}/${filename}_SPAdes_scaffolds.fasta"
        python3 /home/kasutaja/miniconda3/envs/resfinder/bin/run_resfinder.py -ifa "$fasta_file" -o "${OUTPUT_DIR_resfinder}/${filename}_resfinder_results" -s "$resfinder_species" -l 0.6 -t 0.8 --acquired --point -db_res /home/kasutaja/resfinder_db/ -db_point /home/kasutaja/pointfinder_db/
    else
        echo "Query coverage is not numeric or below threshold for $filename. Skipping."
    fi
done

