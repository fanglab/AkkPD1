# Load required libraries
library(vegan)       # For Bray-Curtis distance and PCoA
library(ape)         # For dendrogram
library(dendextend)  # For colored dendrogram
library(ggplot2)     # For PCoA plots

# ---------------------------
# 0. Command-line arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_matrix.csv> <output_prefix> [k_clusters]\n")
}

input_file <- args[1]              # e.g., "/pipeline/data/result_gene_presence_absence.tsv"
output_prefix <- args[2]           # e.g., "/pipeline/data/phylogroup"
k <- ifelse(length(args) >= 3, as.numeric(args[3]), 4)  # Default 4 clusters

# ---------------------------
# 1a. Load and prepare data
# ---------------------------
gene_data <- read.table(input_file, row.names = 1, check.names = FALSE)
gene_matrix <- t(as.matrix(gene_data))

# ---------------------------
# 1b. Load reference metadata
# ---------------------------
ref_metadata_file <- paste0(dirname(input_file), "/reference_phylogroup_metadata.csv")
ref_metadata <- read.csv(ref_metadata_file, stringsAsFactors = FALSE)
ref_metadata$GenomeID <- paste0("REF_", ref_metadata$GenomeID)

# Ensure columns are named correctly
if (!all(c("GenomeID", "Phylogroup") %in% colnames(ref_metadata))) {
  stop("Reference metadata must have columns: GenomeID, Phylogroup")
}

reference_map <- setNames(ref_metadata$Phylogroup, ref_metadata$GenomeID)


# ---------------------------
# 2. Bray-Curtis distance and clustering
# ---------------------------
dist_matrix <- vegdist(gene_matrix, method = "bray")
hc <- hclust(dist_matrix, method = "complete") #change the method to "complete"

# ---------------------------
# 3. Cut into k clusters
# ---------------------------
raw_clusters <- cutree(hc, k = k)

# ---------------------------
# 4. Assign Phylogroup based on reference samples
# ---------------------------

# Initialize phylogroup assignment
phylogroups <- rep(NA, length(raw_clusters))
names(phylogroups) <- names(raw_clusters)

# For each cluster, check which reference genome is present
for (cluster_id in unique(raw_clusters)) {
  samples_in_cluster <- names(raw_clusters)[raw_clusters == cluster_id]
  
  # Keep reference phylogroup as-is
  known_refs <- intersect(samples_in_cluster, names(reference_map))
  unknown_samples <- setdiff(samples_in_cluster, names(reference_map))
  
  if (length(known_refs) > 0) {
    # Get majority phylogroup from known references
    assigned_group <- names(sort(table(reference_map[known_refs]), decreasing = TRUE))[1]
    # Assign to unknown samples only
    phylogroups[unknown_samples] <- assigned_group
  } else {
    # No reference in this cluster → mark all unknowns
    phylogroups[unknown_samples] <- paste0("UnknownCluster", cluster_id)
  }
}
# Now assign the known references directly (if not already done)
phylogroups[names(reference_map)] <- reference_map

# ---------------------------
# 5. Save phylogroup assignment
# ---------------------------
cluster_df <- data.frame(SampleID = names(raw_clusters),
                         ClusterNumber = raw_clusters,
                         Phylogroup = phylogroups)
cluster_df <- cluster_df[order(cluster_df$Phylogroup, cluster_df$SampleID), ]


# ---------------------------
# 6. Colored dendrogram (only!)
# ---------------------------
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = k)

png(paste0(output_prefix, "_colored_dendrogram.png"), width = 800, height = 600)
plot(dend, main = "Colored Dendrogram by Phylogroup")
dev.off()

# ---------------------------
# 7. Full PCoA plot (all samples)
# ---------------------------
pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)
pcoa_df <- data.frame(SampleID = rownames(gene_matrix),
                      PCoA1 = pcoa_result$points[,1],
                      PCoA2 = pcoa_result$points[,2],
                      Phylogroup = phylogroups)
# Add Label column: tested samples only (non-REF)
pcoa_df$Label <- ifelse(!grepl("^REF", pcoa_df$SampleID), pcoa_df$SampleID, NA)
pcoa_df$Label <- sub("(\\.[^.]+)+$", "", pcoa_df$Label)
pcoa_df$PointColor <- ifelse(grepl("^REF", pcoa_df$SampleID), pcoa_df$Phylogroup, "Tested")


pcoa_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = PointColor), size = 3) +
  geom_text(aes(label = Label), color = "black", size = 2.5, vjust = -1, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      "AmIa" = "red",     # or whatever color you want
      "AmIb" = "yellow",
      "AmII" = "blue",
      "AmIII" = "green",
	  "AmIV" = "purple",
      "Tested" = "black"  # <- tested samples always black point
    )
  ) +
  theme_classic() +
  labs(title = "PCoA of Samples (All Phylogroups)",
       x = paste0("PCoA1 (", round(pcoa_result$eig[1]/sum(pcoa_result$eig)*100, 1), "%)"),
       y = paste0("PCoA2 (", round(pcoa_result$eig[2]/sum(pcoa_result$eig)*100, 1), "%)")) +
  theme(legend.title = element_blank())

ggsave(paste0(output_prefix, "_PCoA_all.png"), plot = pcoa_plot, width = 8, height = 6)

# ---------------------------
# 8. PCoA plot only for AmIa and AmIb
# ---------------------------
# Extract AmIa and AmIb samples
amI_samples <- pcoa_df$SampleID[pcoa_df$Phylogroup %in% c("AmIa", "AmIb")]
gene_matrix_amI <- gene_matrix[amI_samples, ]
phylogroups_amI <- phylogroups[rownames(gene_matrix_amI)]

# Clean empty samples
gene_matrix_amI <- gene_matrix_amI[rowSums(gene_matrix_amI) > 0, ]
gene_matrix_amI <- gene_matrix_amI[, colSums(gene_matrix_amI) > 0]

cat("After cleaning: Samples =", nrow(gene_matrix_amI), " Features =", ncol(gene_matrix_amI), "\n")
dist_matrix_amI <- vegdist(gene_matrix_amI, method = "bray")

pcoa_result_amI <- cmdscale(dist_matrix_amI, eig = TRUE, k = 2)  # Classical MDS = PCoA
pcoa_df_amI <- data.frame(SampleID = rownames(gene_matrix_amI),
                      PCoA1 = pcoa_result_amI$points[,1],
                      PCoA2 = pcoa_result_amI$points[,2],
                      Phylogroup = phylogroups_amI)

pcoa_df_amI$Label <- ifelse(!grepl("^REF", pcoa_df_amI$SampleID), pcoa_df_amI$SampleID, NA)
pcoa_df_amI$Label <- sub("(\\.[^.]+)+$", "", pcoa_df_amI$Label)  
pcoa_df_amI$PointColor <- ifelse(grepl("^REF", pcoa_df_amI$SampleID), pcoa_df_amI$Phylogroup, "Tested")


# Step 1: Get coordinates of reference samples
ref_points <- pcoa_df_amI[pcoa_df_amI$SampleID %in% names(reference_map), ]
ref_points$Phylogroup <- reference_map[ref_points$SampleID]

# Step 2: Compute centroids of AmIa and AmIb
centroid_AmIa <- colMeans(ref_points[ref_points$Phylogroup == "AmIa", c("PCoA1", "PCoA2")], na.rm = TRUE)
centroid_AmIb <- colMeans(ref_points[ref_points$Phylogroup == "AmIb", c("PCoA1", "PCoA2")], na.rm = TRUE)

# Step 3: Assign unknowns based on 2D distance to centroids
unknown_samples <- pcoa_df_amI$SampleID[!(pcoa_df_amI$SampleID %in% names(reference_map))]

for (sid in unknown_samples) {
  coords <- pcoa_df_amI[pcoa_df_amI$SampleID == sid, c("PCoA1", "PCoA2")]
  
  if (any(is.na(coords))) {
    phylogroups[sid] <- "Unassigned"
    next
  }

  dist_to_AmIa <- sqrt(sum((coords - centroid_AmIa)^2))
  dist_to_AmIb <- sqrt(sum((coords - centroid_AmIb)^2))

  phylogroups[sid] <- ifelse(dist_to_AmIa < dist_to_AmIb, "AmIa", "AmIb")
}

# Overwrite and export the updated phylogroup assignments
cluster_df$Phylogroup <- phylogroups[cluster_df$SampleID]

write.csv(
  cluster_df,
  paste0(output_prefix, "_assignments.csv"),  # same as before
  row.names = FALSE,
  quote = FALSE
)


# Plot
pcoa_plot_amI <- ggplot(pcoa_df_amI, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = PointColor), size = 3) +
  geom_text(aes(label = Label), color = "black", size = 2.5, vjust = -1, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      "AmIa" = "red",
      "AmIb" = "blue",
      "Tested" = "black"
    )
  ) +
  theme_classic() +  # Better background
  labs(title = "PCoA (AmIa vs AmIb Only)",
       x = paste0("PCoA1 (", round(pcoa_result_amI$eig[1]/sum(pcoa_result_amI$eig)*100, 1), "%)"),
       y = paste0("PCoA2 (", round(pcoa_result_amI$eig[2]/sum(pcoa_result_amI$eig)*100, 1), "%)")) +
  theme(legend.title = element_blank())

ggsave(paste0(output_prefix, "_AmI_sub_PCoA_plot.png"), plot = pcoa_plot_amI, width = 8, height = 6)



# ---------------------------
# Done!
# ---------------------------
cat("Clustering and visualization complete!\n")
cat("Assignments saved to: ", paste0(output_prefix, "_assignments.csv"), "\n")
cat("Dendrograms and PCoA plot saved with prefix: ", output_prefix, "\n")

