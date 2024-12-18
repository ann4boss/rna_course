#Library Installation and Import----------------------------------------------------------
# Install and load required libraries
required_packages <- c("BiocManager", 
                       "DESeq2", 
                       "pheatmap", 
                       "reshape2", 
                       "ggplot2", 
                       "rstudioapi", 
                       "biomaRt", 
                       "ggrepel", 
                       "dplyr", 
                       "clusterProfiler", 
                       "org.Hs.eg.db")

# Function to check and install missing packages
check_and_install <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (pkg == "BiocManager") {
            install.packages(pkg, repos = "https://cloud.r-project.org")
        } else {
            BiocManager::install(pkg, ask = FALSE)
        }
    }
}

# Install any missing packages
for (pkg in required_packages) {
    check_and_install(pkg)
}

# Load the packages
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggplot2)
library(rstudioapi)
library(biomaRt)
library(ggrepel)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

#Working Directory-------------------------------------------------------------------
# Dynamically set the working directory to the location of the script
set_working_directory <- function() {
    if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
        script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
        parent_dir <- dirname(script_dir)
        setwd(parent_dir)
        message("Working directory set to: ", getwd())
    } else {
        message("Unable to dynamically set the working directory. Please set it manually.")
    }
}

set_working_directory()

#Import of featureCount file and meta file-------------------------------------------------
# Define paths
counts_file <- "./data/adjusted_gene_counts.txt"  
metadata_file <- "./data/metadata.csv"
output_dir <- "./analysis"

# Read raw count file
raw_counts <- read.table(counts_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
nrow(raw_counts)
# Clean data: 
# Exclude rows which only have 0, 1 or 2 counts per sample for all samples
sample_columns <- c("HER21", "HER22", "HER23", "NonTNBC1", "NonTNBC2", 
                    "NonTNBC3", "Normal1", "Normal2", "Normal3", "TNBC1", "TNBC2")
clean_counts <- raw_counts[rowSums(raw_counts[sample_columns] > 2) > 0, ]
number_of_excluded_rows <- nrow(raw_counts) - nrow(clean_counts)
# Exclude not needed columns, keep gene ID and sample count for each sample
counts <- clean_counts[, c(1, 7:ncol(clean_counts))]
rownames(counts) <- counts$Geneid
counts <- counts[, -1]  # Remove the 'Geneid' column from the data

# Check if metadata file exists
if (!file.exists(metadata_file)) {
    message("Metadata file not found. Creating metadata...")
    metadata <- data.frame(
        Sample = c("HER21", "HER22", "HER23", "NonTNBC1", "NonTNBC2", "NonTNBC3", 
                   "Normal1", "Normal2", "Normal3", "TNBC1", "TNBC2", "TNBC3"),
        Group = c("HER2_positive", "HER2_positive", "HER2_positive", 
                  "Non_triple_negative", "Non_triple_negative", "Non_triple_negative", 
                  "Normal", "Normal", "Normal", 
                  "Triple_negative", "Triple_negative", "Triple_negative"),
        row.names = 1
    )
    if (!dir.exists("./data")) {
        dir.create("./data", recursive = TRUE)
    }
    write.csv(metadata, metadata_file, row.names = TRUE)
} else {
    message("Metadata file found. Loading metadata...")
    metadata <- read.csv(metadata_file, row.names = 1)
}

# Convert Group column to factor
metadata$Group <- as.factor(metadata$Group)
metadata$Group <- relevel(metadata$Group, ref = "Normal")

# Check metadata consistency
stopifnot(all(rownames(metadata) %in% colnames(counts)))

#DESeq Analysis-----------------------------------------------------------------
# DESeq Analysis
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~ Group)
# Fitting the model to the dds object and estimating size factors, dispersions, and gene-wise statistics
dds <- DESeq(dds)

# Normalize counts to account for differences in sequencing depth and differences in gene length
norm_counts <- counts(dds, normalized = TRUE)

#Comparison if VST or rlog transformation suits this dataset with 12 samples better
# Variance Stabilizing Transformation (assumes a negative binomial distribution)
vst_data <- vst(dds, blind = TRUE)

# PCA after rlog transformation
rlog_data <- rlog(dds, blind = TRUE)


# Compare the transformed data -> PC1 higher for vst, going with vst for further analysis
pca_rlog <- plotPCA(rlog_data, intgroup = "Group")
pca_rlog + ggtitle("PCA - rlog")
pca_vst <- plotPCA(vst_data, intgroup = "Group")
pca_vst + ggtitle("PCA - VST")


# PCA plot
pcaData <- plotPCA(vst_data, intgroup="Group", ntop=500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=Group)) +
    geom_point(size = 3) +
    geom_text(aes(label = rownames(pcaData)), vjust = -1, size = 3) +
    labs(
        x = paste0("PC1: ", percentVar[1], "% variance"),
        y = paste0("PC2: ", percentVar[2], "% variance"),
        title = "PCA Plot"
    ) + 
    coord_fixed() +
    theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path(output_dir, "PCA_plot.png"), plot = pca_plot, width = 8, height = 6)


# Heatmap
sample_dist <- dist(t(assay(vst_data)))
sample_clust <- hclust(sample_dist)
heatmap <- pheatmap(as.matrix(sample_dist), clustering_distance_rows = sample_dist, clustering_distance_cols = sample_dist, main = "Sample Distance Heatmap")
png(file.path(output_dir, "Sample_Heatmap.png"), width = 800, height = 800)
print(heatmap)
dev.off()


#Differential Expression Analysis-----------------------------------------------------------------
# Set pairwise comparison -> decision made with PCA analysis
comparisons <- list(
    c("Group", "Triple_negative", "HER2_positive"),
    c("Group", "Non_triple_negative", "HER2_positive"),
    c("Group", "Triple_negative", "Non_triple_negative")
)

# Initialize results data frame and summary table
results_df <- data.frame()
summary_table <- data.frame()


# Loop through comparisons
for (comp in comparisons) {
    # Get DESeq2 results for the current comparison
    res <- as.data.frame(results(dds, 
                                 contrast = comp,
                                 alpha = 0.05))
    res$comparison <- paste(comp[2], "vs", comp[3], sep = "_")
    res$gene_id <- rownames(res)  # Add gene IDs as a separate column
    
    # Add significance categorization based on adjusted pvalue - alpha level = 0.05
    res$Significance <- "Not Significant"  # Default category
    res$Significance[res$padj < 0.05 & res$log2FoldChange > 0] <- "Upregulated"
    res$Significance[res$padj < 0.05 & res$log2FoldChange < 0] <- "Downregulated"
    
    # Append to the results data frame
    results_df <- rbind(results_df, res)
    
    # Calculate summary statistics
    total_de <- sum(res$padj < 0.05, na.rm = TRUE)
    upregulated <- sum(res$Significance == "Upregulated", na.rm = TRUE)
    downregulated <- sum(res$Significance == "Downregulated", na.rm = TRUE)
    
    # Add to summary table
    summary_table <- rbind(summary_table, data.frame(
        Comparison = paste(comp[2], "vs", comp[3], sep = "_"),
        Total_DE = total_de,
        Up_regulated = upregulated,
        Down_regulated = downregulated
        
    ))
}


# Annotate with gene symbols using biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = unique(results_df$gene_id),
    mart = mart
)

# Merge annotations with results
results_df <- merge(results_df, annotations, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
# Replace empty strings with NA in the 'hgnc_symbol' column
results_df$hgnc_symbol[results_df$hgnc_symbol == ""] <- NA

# Save the results data frame and summary table
write.csv(results_df, file.path(output_dir, "DESeq2_results_all_comparisons.csv"), row.names = FALSE)
write.csv(summary_table, file.path(output_dir, "Summary_table.csv"), row.names = FALSE)


# Heatmap to investigate expression level per sample of top 20 genes
pdf(file.path(output_dir, "Heatmaps_all_comparisons.pdf"))
for (comp in comparisons) {
    # Create a label for the current comparison
    comp_label <- paste(comp[2], "vs", comp[3], sep = "_")
    data <- subset(results_df, comparison == comp_label)
    
    # Filter and sort data for the current comparison
    top_100_genes <- head(
        arrange(filter(data, Significance != "Not Significant" & !is.na(log2FoldChange)), 
                desc(abs(log2FoldChange))), 100
    )
    
    # Extract top 20 genes for heatmap using the normalized counts
    top_gene_ids <- head(top_100_genes$gene_id, 20)
    matched_rows <- match(top_gene_ids, rownames(norm_counts))
    top_gene_norm_counts <- norm_counts[na.omit(matched_rows), ]
    
    # Generate and save the heatmap
    pheatmap(
        log2(top_gene_norm_counts + 1),
        cluster_rows = TRUE,       # Enable clustering of rows (genes)
        cluster_cols = FALSE,      # Enable or not clustering of columns (samples)
        scale = "row",             # Scale genes to z-scores
        clustering_distance_rows = "euclidean",  # Set distance metric-> Euclidean distance measures the straight-line distance between two points
        #clustering_distance_cols = "correlation", # Use correlation for samples
        clustering_method = "average",  # Use average linkage for clustering
        show_rownames = TRUE,
        show_colnames = TRUE,
        #color = colorRampPalette(c("dodgerblue4", "cornsilk1", "firebrick4"))(50),
        main = paste("Heatmap -", comp_label)
    )
}
dev.off()

# Volcano plot for each comparison
# Loop through each comparison and create a volcano plot
pdf(file.path(output_dir, "Volcano_plots_all_comparisons.pdf"))
for (comp in comparisons) {
    # Filter data for the current comparison
    comp_label <- paste(comp[2], "vs", comp[3], sep = "_")
    data <- subset(results_df, comparison == comp_label)
    # Get top 100 genes based on absolute log2FoldChange and padj
    top_100_genes <- head(arrange(filter(data, !is.na(padj) & !is.na(log2FoldChange)), 
                                  desc(abs(log2FoldChange))), 100)
    
    # Highlight the top 20 genes from the top 100
    top_20_genes <- head(top_100_genes, n = 20)
    
    # Create the volcano plot
    volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
        # Points for all genes
        geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
        # Points for labeled genes with a lighter color
        geom_point(data = top_20_genes, aes(x = log2FoldChange, y = -log10(padj)),
                   color = "orange", alpha = 0.8, size = 2) +
        # Add labels with lines connecting to points
        geom_text_repel(data = top_20_genes, aes(label = hgnc_symbol), 
                        size = 3, max.overlaps = 10,
                        segment.color = "black", segment.size = 0.5) + 
        # Customize color palette
        scale_color_manual(values = c(
            "Upregulated" = "firebrick4", 
            "Downregulated" = "dodgerblue4", 
            "Not Significant" = "cornsilk1"
        )) +
        # Set the background to white
        theme_minimal() +
        theme(
            legend.position = "top",
            text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5),
        ) +
        # Add plot labels
        labs(
            title = paste("Volcano Plot -", comp_label),
            x = "Log2 Fold Change",
            y = "-Log10 Adjusted p-value",
            color = "Significance"
        )
    # Print the plot
    print(volcano_plot)
}
dev.off()

# Filter and select top 100 genes for each comparison
# Filter for significant genes and non-NA log2FoldChange
filtered_data <- filter(results_df, Significance != "Not Significant" & !is.na(log2FoldChange))
# Group by comparison
grouped_data <- group_by(filtered_data, comparison)
# Sort by absolute log2 fold change in descending order within each comparison
sorted_data <- arrange(grouped_data, comparison, desc(abs(log2FoldChange)))
# Select top 100 genes per comparison
top_100_genes_per_comparison <- slice_head(sorted_data, n = 100)


#Over-Representation analysis-------------------------------------------------------------------
# Load the background genes from the file
all_genes <- results_df$gene_id
length(all_genes)
# Define a function to perform GO enrichment analysis
perform_GO_enrichment <- function(gene_list, comparison, all_genes, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, pAdjustMethod = "BH") {
    enrich_result <- enrichGO(
        gene         = gene_list,
        universe     = all_genes,             
        OrgDb        = OrgDb,                 # Annotation database for human genes
        ont          = ont,                   # Use selected GO subontology ("BP", "MF", "CC", or "ALL")
        keyType      = "ENSEMBL",             # Gene ID format
        pvalueCutoff = pvalueCutoff,          # p-value cutoff
        pAdjustMethod = pAdjustMethod         # Method Benjamini-Hochberg for multiple testing correction
        #minGSSize    = 10,                   # Minimum number of genes in a GO term to include
        #maxGSSize    = 100,                  # Maximum number of genes in a GO term to include
        #qvalueCutoff = 0.2                   # FDR cutoff
    )
    return(enrich_result)
}

# Subset top 100 genes per comparison and store results
enrichment_results <- list()  # List to store enrichment results for each comparison

for (comparison_pair in comparisons) {
    # Construct the comparison name (e.g., "Triple_negative_vs_HER2_positive")
    comparison_name <- paste(comparison_pair[2], "vs", comparison_pair[3], sep = "_")
    
    # Subset the top 100 genes for each comparison
    top_100_genes <- subset(top_100_genes_per_comparison, comparison == comparison_name)
    gene_id_list <- top_100_genes$gene_id

    # Perform GO enrichment analysis
    enrich_result <- perform_GO_enrichment(gene_list = gene_id_list, comparison = comparison_name, all_genes = all_genes)
    
    # Convert the enrichment result to a data frame
    enrich_result_df <- as.data.frame(enrich_result)
    
    
    # Save the result as a CSV file
    write.csv(enrich_result_df, file.path(output_dir, paste0("GO_enrichment_", comparison_name, ".csv")), row.names = FALSE)
    
    # Store the result
    enrichment_results[[comparison_name]] <- enrich_result
}

#summary(enrich_result)

# Retrieve the result for "Triple_negative_vs_HER2_positive"
#enrich_result_Tripleneg_vs_HER2 <- enrichment_results[["Triple_negative_vs_HER2_positive"]]



pdf(file.path(output_dir, "GO_enrichment_dotplots_all_comparisons.pdf"))
# Loop over each comparison and create a dotplot
for (comparison_name in names(enrichment_results)) {
    # Get the GO enrichment result for the current comparison
    enrich_result <- enrichment_results[[comparison_name]]
    
    # Generate dotplot for the GO enrichment result
    # You can limit the number of terms displayed (e.g., top 20 GO terms)
    dot_plot <- dotplot(enrich_result, showCategory = 20) +
        ggtitle(paste("GO Enrichment for", comparison_name)) +
        theme_minimal()
    
    # Print the plot to the PDF
    print(dot_plot)
}
dev.off()


pdf(file.path(output_dir, "GO_enrichment_barcharts_all_comparisons.pdf"))
# Loop over each comparison and create a bar chart
for (comparison_name in names(enrichment_results)) {
    # Get the GO enrichment result for the current comparison
    enrich_result <- enrichment_results[[comparison_name]]
    
    # Generate a bar chart for the GO enrichment result
    # You can limit the number of terms displayed (e.g., top 20 GO terms)
    bar_plot <- barplot(enrich_result, showCategory = 20) +
        ggtitle(paste("GO Enrichment for", comparison_name)) +
        theme_minimal()
    
    # Print the plot to the PDF
    print(bar_plot)
}
dev.off()
    
