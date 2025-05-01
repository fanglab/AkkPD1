#!/bin/bash
set -e  # Exit immediately on error

# Define directories
RAW_DIR=/pipeline/data/rawdata
KNEADDATA_DIR=/pipeline/data/kneaddata_results
MAPPING_DIR=/pipeline/data/mapping_results
DB_DIR=/pipeline/DB
PANGENOME_TSV=/pipeline/DB/Akkermansia_muciniphila_pangenome.tsv  
PANGENOME_IDX=/pipeline/DB/Akkermansia_muciniphila  

# Create output directories
mkdir -p ${RAW_DIR} ${KNEADDATA_DIR} ${MAPPING_DIR}

# Step 1: Download raw reads
if [ -z "$(ls ${RAW_DIR}/*.fastq* 2>/dev/null)" ]; then
  echo "‚ùå ERROR:raw data not found."
  exit 1
fi

# Step 2: Run Kneaddata (single-end)
echo "Ì†æÌ∑π Running Kneaddata (single-end mode)..."
for sample in $(ls ${RAW_DIR}/*.fastq | sed 's/\.fastq//' | xargs -n1 basename)
do
  echo "Ì†æÌ∑π Cleaning sample: ${sample}"
  kneaddata \
    -db ${DB_DIR} \
    -o ${KNEADDATA_DIR} \
    -i ${RAW_DIR}/${sample}.fastq \
    --trimmomatic /opt/conda/envs/akk_pipeline/share/trimmomatic-0.39-2
  
  # Move cleaned file to a standardized name
  mv ${KNEADDATA_DIR}/${sample}_kneaddata.repeats.removed.fastq ${KNEADDATA_DIR}/${sample}_kneaddata.fastq || true
  rm -f ${KNEADDATA_DIR}/${sample}_kneaddata.trimmed.fastq
  rm -f ${KNEADDATA_DIR}/reformatted_identifiers*
  
  echo "‚úÖ Cleaned sample saved as: ${KNEADDATA_DIR}/${sample}kneaddata.fastq"
done

# Step 3: Map reads to pangenome (single-end)
echo "Ì†æÌ∑¨ Mapping reads to pangenome..."
for sample in $(ls ${KNEADDATA_DIR}/*_kneaddata.fastq | sed 's/_kneaddata.fastq//' | xargs -n1 basename)
do
  echo "Ì†æÌ∑¨ Mapping sample: ${sample}"
  panphlan_map.py \
    -p ${PANGENOME_TSV} \
    --indexes ${PANGENOME_IDX} \
    -i ${KNEADDATA_DIR}/${sample}_kneaddata.fastq \
    -o ${MAPPING_DIR}/${sample}.csv
done

# Step 4: Profile gene presence/absence
echo "Ì†æÌ∑¨ Profiling genes across samples..."
panphlan_profiling.py \
  --min_coverage 1 \
  --left_max 1.7 \
  --right_min 0.3 \
  -i ${MAPPING_DIR} \
  -p ${PANGENOME_TSV} \
  --o_matrix /pipeline/data/result_gene_presence_absence.tsv \
  --add_ref

# Step 5: Phylogroup assignment
echo "Ì†ºÌºø Running phylogroup assignment..."
Rscript akk_phylo_assign/assign_phylogroups.R \
    /pipeline/data/result_gene_presence_absence.csv \
    /pipeline/data/phylogroup \
    4


echo "‚úÖ Phylogroup assignment completed! Results saved to /pipeline/data/phylogroup_assignments.csv"
echo "Ì†ºÌæâ Pipeline complete!"

