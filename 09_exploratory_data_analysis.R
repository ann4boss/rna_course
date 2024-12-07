# Load required libraries
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
BiocManager::install("reshape2")
library("DESeq2")
library("pheatmap")
library("ggplot2")
library("reshape2")
library("ggplot2")

#set working directory to source of script
setwd("~/Library/CloudStorage/OneDrive-UniversitédeFribourg/01_First Semester/467713 RNA sequencing/R_scripts_output")

# Set file paths and experimental groups
counts_file <- "./gene_counts.txt"  
output_dir <- "./output"

# Read the raw count file
raw_counts <- read.table(counts_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(raw_counts)
# Remove unnecessary header lines (keeping rows where Geneid is non-NA)
clean_counts <- raw_counts[!is.na(raw_counts$Geneid), ]
# Remove unnecessary columns, keep only the "Geneid" and the columns containing sample counts
counts <- clean_counts[, c(1, 7:ncol(clean_counts))]
# Rename the columns for clarity (removing unnecessary path information)
colnames(counts)[2:ncol(counts)] <- gsub("X.data.users.aboss.rna_course.06_mapping_results_uncleaned_stranded.", "", colnames(counts)[2:ncol(counts)])
colnames(counts)[2:ncol(counts)] <- gsub("_aligned.sorted.bam", "", colnames(counts)[2:ncol(counts)])
# Save "Geneid" as row names
rownames(counts) <- counts$Geneid
counts <- counts[, -1]
counts

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
write.csv(metadata, "./metadata.csv", row.names = TRUE)

# Check that rownames of metadata match colnames of counts
stopifnot(all(rownames(metadata) %in% colnames(counts)))

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Group)

# Run DESeq to normalize counts and estimate dispersion
dds <- DESeq(dds)

#to get all gene IDs
all_genes <- rownames(dds)
# Save to a file
write.table(all_genes, file.path(output_dir, "All_Genes_in_DESeq2.txt"), 
            quote = FALSE, row.names = FALSE, col.names = "GeneID")

# Variance Stabilizing Transformation (vst)
vst_data <- vst(dds, blind = TRUE)


# PCA plot for sample clustering
pca <- plotPCA(vst_data, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

# without labels
pca_plot <- ggplot(pca, aes(PC1, PC2, color=group, shape=type)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       title = "PCA of Samples") +
    coord_fixed()

# Add sample labels
pca$Sample <- rownames(pca)
pca_plot <- ggplot(pca, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 3) +       # Points for each sample
    geom_text(vjust = -0.8) +    # Sample labels slightly above the points
    labs(title = "PCA Plot", x = "PC1", y = "PC2") +
    theme_minimal() +
    theme_minimal()


ggsave(filename = file.path(output_dir, "PCA_plot.png"), plot = pca_plot, width = 8, height = 6)

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



###------------Differential Expression Analysis---------------------------------
resultsNames(dds)

# Define pairwise comparisons
comparisons <- list(
    c("Group", "Triple_negative", "HER2_positive"),
    c("Group", "Non_triple_negative", "HER2_positive"),
    c("Group", "Triple_negative", "Non_triple_negative")
)

# Initialize summary file
summary_file <- file.path(output_dir, "DE_summary_statistics.csv")
summary_data <- data.frame(
    Comparison = character(),
    Total_DE = integer(),
    Up_regulated = integer(),
    Down_regulated = integer(),
    stringsAsFactors = FALSE
)

# Loop through each comparison
for (comparison in comparisons) {
    # Extract the name of the comparison
    comp_name <- paste(comparison[3], "vs", comparison[2])
    
    # Get results for the comparison
    res <- results(dds, contrast = comparison)
    
    # Order results by adjusted p-value
    res <- res[order(res$padj), ]
    
    # Save all results to a file
    results_file <- file.path(output_dir, paste0("DE_results_", comp_name, ".csv"))
    write.csv(as.data.frame(res), results_file, row.names = TRUE)
    
    # Filter significant genes
    sig_res <- res[!is.na(res$padj) & res$padj < 0.05, ]
    
    # Count DE genes
    total_DE <- nrow(sig_res)
    up_regulated <- nrow(sig_res[sig_res$log2FoldChange > 0, ])
    down_regulated <- nrow(sig_res[sig_res$log2FoldChange < 0, ])
    
    # Append statistics to summary data
    summary_data <- rbind(summary_data, data.frame(
        Comparison = comp_name,
        Total_DE = total_DE,
        Up_regulated = up_regulated,
        Down_regulated = down_regulated
    ))
    
    # Print summary to console
    cat("Comparison: ", comp_name, "\n")
    cat("Total DE genes: ", total_DE, "\n")
    cat("Up-regulated genes: ", up_regulated, "\n")
    cat("Down-regulated genes: ", down_regulated, "\n")
}

# Write summary statistics to a file
write.csv(summary_data, summary_file, row.names = FALSE)
cat("Summary statistics saved to: ", summary_file, "\n")




# Function to create and save a volcano plot
create_volcano_plot <- function(res, comp_name) {
    res$Significant <- ifelse(!is.na(res$padj) & res$padj < 0.05, "Yes", "No")
    
    volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
        geom_point(alpha = 0.7) +
        scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
        labs(
            title = paste("Volcano Plot:", comp_name),
            x = "Log2 Fold Change",
            y = "-Log10(p-value)"
        ) +
        theme_minimal()
    
    # Save plot
    ggsave(
        filename = file.path(output_dir, paste0("Volcano_Plot_", comp_name, ".png")),
        plot = volcano,
        width = 8, height = 6
    )
    
    cat("Volcano plot saved for comparison:", comp_name, "\n")
}
# Generate volcano plots for each comparison
for (comparison in comparisons) {
    comp_name <- paste(comparison[3], "vs", comparison[2])
    res <- results(dds, contrast = comparison)
    create_volcano_plot(as.data.frame(res), comp_name)
}


# Function to create and save an MA plot
create_ma_plot <- function(res, comp_name) {
    ma <- ggplot(res, aes(x = baseMean, y = log2FoldChange, color = padj < 0.05)) +
        geom_point(alpha = 0.7) +
        scale_x_log10() +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
        labs(
            title = paste("MA Plot:", comp_name),
            x = "Mean Normalized Count",
            y = "Log2 Fold Change"
        ) +
        theme_minimal()
    
    # Save plot
    ggsave(
        filename = file.path(output_dir, paste0("MA_Plot_", comp_name, ".png")),
        plot = ma,
        width = 8, height = 6
    )
    
    cat("MA plot saved for comparison:", comp_name, "\n")
}
# Generate MA plots for each comparison
for (comparison in comparisons) {
    comp_name <- paste(comparison[3], "vs", comparison[2])
    res <- results(dds, contrast = comparison)
    create_ma_plot(as.data.frame(res), comp_name)
}


# Function to create and save a heatmap
create_heatmap <- function(dds, sig_genes, comp_name) {
    if (nrow(sig_genes) > 0) {
        # Get normalized counts for significant genes
        norm_counts <- counts(dds, normalized = TRUE)[rownames(sig_genes), ]
        heatmap_file <- file.path(output_dir, paste0("Heatmap_", comp_name, ".png"))
        
        # Plot and save heatmap
        png(heatmap_file, width = 800, height = 800)
        pheatmap(
            log2(norm_counts + 1),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            scale = "row",
            main = paste("Heatmap of Significant Genes:", comp_name)
        )
        dev.off()
        cat("Heatmap saved for comparison:", comp_name, "\n")
    } else {
        cat("No significant genes for heatmap in comparison:", comp_name, "\n")
    }
}
# Generate heatmaps for each comparison
for (comparison in comparisons) {
    comp_name <- paste(comparison[3], "vs", comparison[2])
    res <- results(dds, contrast = comparison)
    sig_genes <- res[!is.na(res$padj) & res$padj < 0.05, ]
    create_heatmap(dds, sig_genes, comp_name)
}

 