# Load required libraries
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
library("DESeq2")
library("pheatmap")
library("ggplot2")


# Set file paths and experimental groups
counts_file <- "/Users/annaboss/Desktop/RNA seq Breast cancer/gene_counts.txt"  
output_dir <- "/Users/annaboss/Desktop/RNA seq Breast cancer/"

# Read the raw count file
raw_counts <- read.table(counts_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(raw_counts)
# Remove unnecessary header lines (keeping rows where Geneid is non-NA)
clean_counts <- raw_counts[!is.na(raw_counts$Geneid), ]
# Remove unnecessary columns, keep only the "Geneid" and the columns containing sample counts
counts <- clean_counts[, c(1, 7:ncol(clean_counts))]
# Rename the columns for clarity (removing unnecessary path information)
colnames(counts)[2:ncol(counts)] <- gsub("X.data.users.aboss.rna_course.06_mapping_results.", "", colnames(counts)[2:ncol(counts)])
colnames(counts)[2:ncol(counts)] <- gsub("_aligned.sorted.bam", "", colnames(counts)[2:ncol(counts)])
# Save "Geneid" as row names
rownames(counts) <- counts$Geneid
counts <- counts[, -1]


# Create metadata within the script
metadata <- data.frame(
    Sample = c("HER21", "HER22", "HER23", "NonTNBC1", "NonTNBC2", "NonTNBC3", 
               "Normal1", "Normal2", "Normal3", "TNBC1", "TNBC2", "TNBC3"),
    Group = c("HER2_positive", "HER2_positive", "HER2_positive", 
              "Non_triple_negative", "Non_triple_negative", "Non_triple_negative", 
              "Normal", "Normal", "Normal", 
              "Triple_negative", "Triple_negative", "Triple_negative"), row.names = 1
)
metadata
rownames(metadata)
colnames(counts)

# Save metadata as a CSV file for reference (optional)
write.csv(metadata, "metadata.csv", row.names = TRUE)

# Check that rownames of metadata match colnames of counts
stopifnot(all(rownames(metadata) %in% colnames(counts)))

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Group)

# Run DESeq to normalize counts and estimate dispersion
dds <- DESeq(dds)

# Variance Stabilizing Transformation (vst)
vst_data <- vst(dds, blind = TRUE)

# PCA plot for sample clustering
pca <- plotPCA(vst_data, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

pca_plot <- ggplot(pca, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       title = "PCA of Samples") +
  theme_minimal()

ggsave(filename = file.path(output_dir, "PCA_plot_2.png"), plot = pca_plot, width = 8, height = 6)

# Hierarchical clustering with heatmap
sample_dist <- dist(t(assay(vst_data)))
sample_clust <- hclust(sample_dist)
heatmap <- pheatmap(as.matrix(sample_dist), 
                    clustering_distance_rows = sample_dist,
                    clustering_distance_cols = sample_dist,
                    main = "Sample Distance Heatmap")

# Save heatmap to file
png(file.path(output_dir, "Sample_Heatmap.png"), width = 800, height = 800)
print(heatmap)
dev.off()

