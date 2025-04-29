#!/bin/bash
echo "Running minimal pipeline test..."


# Phylogroup assignment
echo "🌿 Running phylogroup assignment..."
Rscript akk_phylo_assign/assign_phylogroups.R \
    /pipeline/data/example_gene_presence_absence.csv \
    /pipeline/data/phylogroup \
    4

echo "✅ Phylogroup assignment completed! Results saved to /pipeline/data/phylogroup_assignments.csv"
echo "������ Pipeline complete!"
