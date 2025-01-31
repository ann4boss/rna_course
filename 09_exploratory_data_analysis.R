# Library Installation and Import ----------------------------------------------------------

# Install and load required libraries
required_packages <- c(
    "BiocManager", "rstudioapi", "DESeq2", "clusterProfiler", "org.Hs.eg.db",
    "EnhancedVolcano", "enrichplot", "ggplot2", "pheatmap",
    "ggrepel", "DOSE", "gridExtra"
)

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
library(rstudioapi)         # To set working directory dynamically
library(DESeq2)             # For differential expression analysis
library(ggplot2)            # For building graphs
library(gridExtra)          # For arranging multiple plots
library(EnhancedVolcano)    # For creating volcano plots
library(ggrepel)            # For better text labeling in plots, used for volcano plots
library(clusterProfiler)    # For GO enrichment analysis
library(enrichplot)         # To plot GO enrichment results
library(org.Hs.eg.db)       # For human gene annotation, gene symbols
library(pheatmap)           # For heatmaps
library(DOSE)               # For disease ontology analysis

# Working Directory -----------------------------------------------------------------------

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

# Import Data -----------------------------------------------------------------------------

# Define paths
counts_file <- "./data/adjusted_gene_counts.txt"
metadata_file <- "./data/metadata.csv"
output_dir <- "./analysis"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Read raw count file
raw_counts <- read.table(counts_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message("Number of genes in raw counts: ", nrow(raw_counts))

# Clean data: Keep gene ID and sample counts
counts <- raw_counts[, c(1, 7:ncol(raw_counts))]
rownames(counts) <- counts$Geneid
counts <- counts[, -1]  # Remove the 'Geneid' column

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
    write.csv(metadata, metadata_file, row.names = TRUE)
} else {
    message("Metadata file found. Loading metadata...")
    metadata <- read.csv(metadata_file, row.names = 1)
}

# Convert Group column to factor and set "Normal" as the reference level
metadata$Group <- as.factor(metadata$Group)
metadata$Group <- relevel(metadata$Group, ref = "Normal")

# Check metadata consistency with counts table
stopifnot(all(rownames(metadata) %in% colnames(counts)))

# DESeq2 Analysis -------------------------------------------------------------------------

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Group)

# Filter low-count genes
initial_rows <- nrow(dds)
dds <- dds[rowSums(counts(dds)) > 10, ]
filtered_rows <- initial_rows - nrow(dds)
message("Number of genes filtered out: ", filtered_rows)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Normalize counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(output_dir, "Normalized_counts.csv"), row.names = TRUE)

# Compare VST and rlog transformations -> PC1 higher for vst, going with vst for further analysis
vst_data <- vst(dds, blind = TRUE)
rlog_data <- rlog(dds, blind = TRUE)
pca_rlog <- plotPCA(rlog_data, intgroup = "Group")
pca_rlog + ggtitle("PCA - rlog")
pca_vst <- plotPCA(vst_data, intgroup = "Group")
pca_vst + ggtitle("PCA - VST")


# PCA Analysis ----------------------------------------------------------------------------

# Define consistent colors for each group
group_colors <- c(
    "HER2_positive" = "dodgerblue4",
    "Non_triple_negative" = "orange",
    "Triple_negative" = "firebrick4",
    "Normal" = "darkolivegreen"
)

# Function to create PCA plots
create_pca_plot <- function(vst_data, title, tag) {
    pcaData <- plotPCA(vst_data, intgroup = "Group", ntop = 500, returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    ggplot(pcaData, aes(PC1, PC2, color = Group)) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(  # Use ggrepel for dynamic labeling
            aes(label = rownames(pcaData)), 
            size = 3, 
            box.padding = 0.3,
            arrow = NULL
        
        ) +
        labs(
            x = paste0("PC1: ", percentVar[1], "% variance"),
            y = paste0("PC2: ", percentVar[2], "% variance"),
            title = title,
            tag = tag
        ) +
        scale_color_manual(values = group_colors) +
        coord_fixed() +
        theme(
            plot.title = element_text(hjust = 0.5),
            plot.tag = element_text(size = 14, face = "bold")
        ) +
        xlim(c(-55, 80)) +
        ylim(c(-50,40))
    
}

# Create PCA plots
pca_plot_all <- create_pca_plot(vst_data, "PCA Plot (All Groups)", "A")
vst_data_subset <- vst_data[, vst_data$Group != "Normal"]
pca_plot_subset <- create_pca_plot(vst_data_subset, "PCA Plot (Excluding Normal Group)", "B")


# Combine and save PCA plots
combined_pca_plots <- grid.arrange(pca_plot_all, pca_plot_subset, ncol = 2)
ggsave(file.path(output_dir, "PCA_plot.png"), combined_pca_plots, width = 14, height = 8)
dev.off()


# Differential Expression Analysis --------------------------------------------------------

# Define pairwise comparisons
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
    res <- as.data.frame(results(dds, contrast = comp, alpha = 0.05, pAdjustMethod = "BH"))
    res$comparison <- paste(comp[2], "vs", comp[3], sep = "_")
    res$gene_id <- rownames(res)
    
    # Add significance categorization
    res$Significance <- "Not Significant"
    res$Significance[res$padj < 0.05 & res$log2FoldChange > 0] <- "Upregulated"
    res$Significance[res$padj < 0.05 & res$log2FoldChange < 0] <- "Downregulated"
    
    # Append to results data frame
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

results_df$gene_symbol <- mapIds(
    org.Hs.eg.db,
    keys = results_df$gene_id,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
)


# Replace NA gene symbols with the ensembl gene ID
results_df$gene_symbol <- ifelse(is.na(results_df$gene_symbol), results_df$gene_id, results_df$gene_symbol)

# or replace them with Not Available
#results_df$hgnc_symbol[results_df$hgnc_symbol == ""] <- "Not Available"

# Save results and adjust rownames to simple numbers
rownames(results_df) <- NULL
write.csv(results_df, file.path(output_dir, "DESeq2_results_all_comparisons.csv"), row.names = FALSE)
write.csv(summary_table, file.path(output_dir, "Summary_table.csv"), row.names = FALSE)

# Volcano Plots ---------------------------------------------------------------------------

png(file.path(output_dir, "Volcano_plots_enhanced_all_comparisons.png"), 
    width = 16, height = 16, units = "in", res = 300)

volcano_plots <- list()

# Create Volcano plot for each comparison
for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    comp_label <- paste(comp[2], "vs", comp[3], sep = "_")
    data <- subset(results_df, comparison == comp_label)
    
    volcano_plot <- EnhancedVolcano(
        data,
        lab = data$gene_symbol,
        x = 'log2FoldChange',
        y = 'padj',
        xlim = c(-10, 10),
        ylim = c(0, 40),
        title = comp_label,
        subtitle = "",
        subtitleLabSize = 1,
        legendLabels = c("Not Sig.", expression(Log[2] ~ FC), "adj. p-value", 
                         expression(adj.p-value ~ and ~ log[2] ~ FC)),
        ylab = bquote(~-Log[10] ~ italic(adj.p-value)),
        pCutoff = 0.05,
        FCcutoff = 2,
        pointSize = 2.0,
        labSize = 4.0,
        col = c('cornsilk1', 'dodgerblue4', 'orange', 'firebrick4'),
        colAlpha = 1,
        gridlines.major = FALSE,
        gridlines.minor = FALSE,
        drawConnectors = TRUE,
        widthConnectors = 0.75,
        colConnectors = 'black',
        max.overlaps = 20
    ) +
        labs(tag = LETTERS[i]) +
        theme(plot.tag = element_text(size = 14, face = "bold"))
    
    volcano_plots[[i]] <- volcano_plot
}

combined_volcano_plots <- grid.arrange(grobs = volcano_plots, ncol = 2)
print(combined_volcano_plots)
dev.off()

# Heatmap Analysis ------------------------------------------------------------------------
# Heatmap for all significant genes
# Set output file for heatmaps
png(file.path(output_dir, "Heatmaps.png"), 
    width = 20, height = 10, units = "in", res = 300) 

# Filter for all significant genes with absolute log2 fold change ≥ 2
significant_genes <- results_df[results_df$Significance != "Not Significant" & abs(results_df$log2FoldChange) >= 2, ]

# Extract normalized counts for all breast cancer samples (excluding normal samples)
norm_counts_without_normal <- norm_counts[, c(1:6, 10:12)]

# Subset normalized counts for significant genes
significant_norm_counts <- norm_counts_without_normal[significant_genes$gene_id, , drop = FALSE]

# Convert to matrix
significant_norm_counts_matrix <- as.matrix(significant_norm_counts)

# Create heatmap for all significant genes
all_genes_heatmap <- pheatmap(
    significant_norm_counts_matrix,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean", 
    clustering_distance_cols = "euclidean",  
    clustering_method = "complete",
    border_color = NA, 
    fontsize = 16, 
    scale = "row",
    treeheight_row = 0,
    treeheight_col = 20,
    main = "A"
)

# Heatmap for top 50 DEGs
top_significant_genes <- significant_genes[order(-abs(significant_genes$log2FoldChange)), ][1:50, ]
top_significant_genes_norm_counts <- norm_counts_without_normal[top_significant_genes$gene_id, , drop = FALSE]

# Convert to matrix
top_significant_genes_norm_counts_matrix <- as.matrix(top_significant_genes_norm_counts)

# Map gene IDs to gene symbols
gene_ids <- rownames(top_significant_genes_norm_counts_matrix)
gene_symbols <- mapIds(
    org.Hs.eg.db,       
    keys = gene_ids,     
    keytype = "ENSEMBL", 
    column = "SYMBOL",   
    multiVals = "first"
)

# Replace NA values with the original gene IDs
gene_labels <- ifelse(is.na(gene_symbols), gene_ids, gene_symbols)
rownames(top_significant_genes_norm_counts_matrix) <- gene_labels

# Create heatmap for top 50 DEGs
heatmap_top_50 <- pheatmap(
    top_significant_genes_norm_counts_matrix,  
    clustering_distance_rows = "euclidean", 
    clustering_distance_cols = "euclidean",  
    clustering_method = "complete",          
    show_rownames = TRUE,                    
    show_colnames = TRUE,
    scale = "row",
    treeheight_row = 20,
    treeheight_col = 10,
    fontsize_row = 11,                        
    fontsize_col = 16, 
    cellheight = 10,
    cellwidth = 20,
    main = "B"
)

# Arrange both heatmaps in one figure
combined_heatmap_plots <- grid.arrange(all_genes_heatmap[[4]], heatmap_top_50[[4]], ncol = 2)
print(combined_heatmap_plots)
dev.off()


# GO Enrichment Analysis ------------------------------------------------------------------

# Extract all unique gene IDs from the results dataframe to use as the background (universe) for GO analysis
all_genes <- unique(results_df$gene_id)
length(all_genes)

# Define a list of comparisons to analyze
comparisons_for_GO <- list(
    "Triple_negative_vs_HER2_positive",
    "Non_triple_negative_vs_HER2_positive",
    "Triple_negative_vs_Non_triple_negative"
)

# Loop through each comparison to perform GO enrichment analysis
for (comparison_name in comparisons_for_GO) {
    
    # Filter genes for the current comparison:
    # Keep only rows corresponding to the current comparison
    # Retain genes that are statistically significant (Significance != "Not Significant")
    # Retain genes with a log2 fold change of at least 2 (indicating strong differential expression)
    filtered_genes <- subset(results_df, 
                             comparison == comparison_name & 
                                 Significance != "Not Significant" & 
                                 abs(log2FoldChange) >= 2)
    
    # Extract the gene IDs of the filtered genes
    gene_ids <- filtered_genes$gene_id
    
    message("The count of genes inputed in enrich GO is ", length(gene_ids), " for ", comparison_name)
    
    # Perform GO enrichment analysis using the enrichGO function
    enrich_result <- enrichGO(
        gene = gene_ids,               # List of differentially expressed gene IDs
        universe = all_genes,          # Background set of genes
        OrgDb = "org.Hs.eg.db",        # Organism database (human)
        ont = "BP",                    # GO ontology: Biological Process
        keyType = "ENSEMBL",           # Gene ID type: Ensembl IDs
        pvalueCutoff  = 0.05,          # Significance threshold for GO terms
        pAdjustMethod = "BH"           # Method for p-value adjustment (Benjamini-Hochberg)
    )
    
    
    # Check if the enrichment analysis returned valid results
    if (!is.null(enrich_result)) {
        
        # Convert the enrichment results to a dataframe
        enrich_df <- as.data.frame(enrich_result)
        message("the number of GO terms ", nrow(enrich_df), " for the comparison: ", comparison_name)
        
        # Save the GO enrichment results to a CSV file
        write.csv(enrich_df, file.path(output_dir, paste0("GO_Enrichment_", comparison_name, ".csv")), row.names = FALSE)
        
        # Generate a bar plot of the top 20 enriched GO terms
        bar_plot <- barplot(enrich_result, showCategory = 20) +
            ggtitle(paste("GO Enrichment for", comparison_name))  # Add a title to the plot
        
        # Save the bar plot as a PNG file
        ggsave(file.path(output_dir, paste0("GO_Enrichment_", comparison_name, ".png")), bar_plot, width = 10, height = 6)
        
        # Convert Ensembl gene IDs to gene symbols for easier interpretation
        enrich_result_readable <- setReadable(enrich_result, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL")
        
        # Create a named vector of log2 fold changes for the filtered genes
        geneList <- filtered_genes$log2FoldChange
        names(geneList) <- filtered_genes$gene_id
        
        # Generate a gene-concept network plot (cnetplot)
        net_plot <- cnetplot(enrich_result_readable, foldChange = geneList, circular = TRUE, colorEdge = TRUE)
               
        # Save the gene-concept network plot as a PNG file
        ggsave(file.path(output_dir, paste0("Gene_Concept_Network_", comparison_name, ".png")), net_plot, width = 16, height = 8)
    }
}
 

