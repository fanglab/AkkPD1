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
  echo "ERROR:raw data not found."
  exit 1
fi

# ----------------------------
# Parallel workers (Option B)
# ----------------------------
# Defaults: use all cores if nproc exists, else 4. Override with: JOBS=8 ./run_akk_pipeline.sh
JOBS="${JOBS:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 4)}"

# Step 2: Run Kneaddata (single-end)
echo "Running Kneaddata (single-end mode)..."

process_one_kneaddata() {
  sample="$1"
  echo "Cleaning sample: ${sample}"
  kneaddata \
    -db ${DB_DIR} \
    -o ${KNEADDATA_DIR} \
    -i ${RAW_DIR}/${sample}.fastq \
    --threads 12 \
    --trimmomatic /opt/conda/envs/akk_pipeline/share/trimmomatic-0.39-2

  # Move cleaned file to a standardized name (unchanged behavior)
  mv ${KNEADDATA_DIR}/${sample}_kneaddata.repeats.removed.fastq ${KNEADDATA_DIR}/${sample}_kneaddata.fastq || true
  rm -f ${KNEADDATA_DIR}/${sample}_kneaddata.trimmed.fastq
  rm -f ${KNEADDATA_DIR}/reformatted_identifiers*_${sample}

  echo "Cleaned sample saved as: ${KNEADDATA_DIR}/${sample}kneaddata.fastq"
}
export -f process_one_kneaddata
export RAW_DIR KNEADDATA_DIR DB_DIR

# Build the exact same sample list as before (only .fastq, unchanged)
tmp_knead_list="$(mktemp)"
ls ${RAW_DIR}/*.fastq | sed 's/\.fastq//' | xargs -n1 basename > "${tmp_knead_list}"

# Prefer GNU parallel; fallback to xargs -P
if command -v parallel >/dev/null 2>&1; then
  parallel --will-cite -j "${JOBS}" process_one_kneaddata :::: "${tmp_knead_list}"
else
  xargs -n1 -P "${JOBS}" -I{} bash -lc 'process_one_kneaddata "$@"' _ {} < "${tmp_knead_list}"
fi
rm -f "${tmp_knead_list}"

# Step 3: Map reads to pangenome (single-end)
echo "Mapping reads to pangenome..."

process_one_map() {
  sample="$1"
  echo "Mapping sample: ${sample}"
  panphlan_map.py \
    -p ${PANGENOME_TSV} \
    --indexes ${PANGENOME_IDX} \
    -i ${KNEADDATA_DIR}/${sample}_kneaddata.fastq \
    -o ${MAPPING_DIR}/${sample}.csv
}
export -f process_one_map
export KNEADDATA_DIR MAPPING_DIR PANGENOME_TSV PANGENOME_IDX

# Recreate the same sample base names as before
tmp_map_list="$(mktemp)"
ls ${KNEADDATA_DIR}/*_kneaddata.fastq | sed 's/_kneaddata.fastq//' | xargs -n1 basename > "${tmp_map_list}"

if command -v parallel >/dev/null 2>&1; then
  parallel --will-cite -j "${JOBS}" process_one_map :::: "${tmp_map_list}"
else
  xargs -n1 -P "${JOBS}" -I{} bash -lc 'process_one_map "$@"' _ {} < "${tmp_map_list}"
fi
rm -f "${tmp_map_list}"

# Step 4: Profile gene presence/absence (unchanged)
echo "Profiling genes across samples..."
panphlan_profiling.py \
  --min_coverage 1 \
  --left_max 1.7 \
  --right_min 0.3 \
  -i ${MAPPING_DIR} \
  -p ${PANGENOME_TSV} \
  --o_matrix /pipeline/data/result_gene_presence_absence.tsv \
  --add_ref

# Step 5: Phylogroup assignment (unchanged)
echo "Running phylogroup assignment..."
Rscript akk_phylo_assign/assign_phylogroups.R \
    /pipeline/data/result_gene_presence_absence.tsv \
    /pipeline/data/phylogroup \
    4

echo "Phylogroup assignment completed! Results saved to /pipeline/data/phylogroup_assignments.csv"
echo "Pipeline complete!"

