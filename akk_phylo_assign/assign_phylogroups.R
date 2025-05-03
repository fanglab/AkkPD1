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

input_file <- args[1]              # e.g., "/pipeline/data/result_gene_presence_absence.csv"
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

# Ensure columns are named correctly
if (!all(c("GenomeID", "Phylogroup") %in% colnames(ref_metadata))) {
  stop("Reference metadata must have columns: GenomeID, Phylogroup")
}

reference_map <- setNames(ref_metadata$Phylogroup, ref_metadata$GenomeID)


# ---------------------------
# 2. Bray-Curtis distance and clustering
# ---------------------------
dist_matrix <- vegdist(gene_matrix, method = "bray")
hc <- hclust(dist_matrix, method = "ward.D2")

# ---------------------------
# 3. Cut into k clusters
# ---------------------------
raw_clusters <- cutree(hc, k = k)

# ---------------------------
# 4. Assign Phylogroup based on reference samples
# ---------------------------
# Define known references
reference_map <- list(
  "REF_GCA_000020225" = "AmIa",
  "REF_GCA_002884915" = "AmIII",
  "REF_GCA_002884995" = "AmII",
  "REF_GCA_002885155" = "AmIb"
)

# Initialize phylogroup assignment
phylogroups <- rep(NA, length(raw_clusters))
names(phylogroups) <- names(raw_clusters)

# For each cluster, check which reference genome is present
for (cluster_id in unique(raw_clusters)) {
  samples_in_cluster <- names(raw_clusters)[raw_clusters == cluster_id]
  
  matched_refs <- intersect(samples_in_cluster, names(reference_map))
  
  if (length(matched_refs) > 0) {
    assigned_group <- reference_map[[matched_refs[1]]]
    phylogroups[samples_in_cluster] <- assigned_group
  } else {
    phylogroups[samples_in_cluster] <- paste0("UnknownCluster", cluster_id)
  }
}

# ---------------------------
# 5. Save phylogroup assignment
# ---------------------------
cluster_df <- data.frame(SampleID = names(raw_clusters),
                         ClusterNumber = raw_clusters,
                         Phylogroup = phylogroups)
cluster_df <- cluster_df[order(cluster_df$Phylogroup, cluster_df$SampleID), ]

write.csv(cluster_df, paste0(output_prefix, "_assignments.csv"),
          row.names = FALSE, quote = FALSE)

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

