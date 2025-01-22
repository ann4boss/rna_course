# Library Installation and Import ----------------------------------------------------------

# Install and load required libraries
required_packages <- c(
    "BiocManager", "rstudioapi", "DESeq2", "clusterProfiler", "org.Hs.eg.db",
    "EnhancedVolcano", "enrichplot", "ggplot2", "biomaRt", "pheatmap", "reshape2",
    "ggrepel", "dplyr", "DOSE", "gridExtra"
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
library(rstudioapi)        # To set working directory dynamically
library(DESeq2)            # For differential expression analysis
library(ggplot2)           # For building graphs
library(gridExtra)         # For arranging multiple plots
library(EnhancedVolcano)   # For creating volcano plots
library(ggrepel)           # For better text labeling in plots
library(clusterProfiler)   # For GO enrichment analysis
library(enrichplot)        # To plot GO enrichment results
library(org.Hs.eg.db)      # For human gene annotation
library(biomaRt)           # For getting gene symbols
library(pheatmap)          # For heatmaps
library(reshape2)          # For data reshaping
library(dplyr)             # For data manipulation
library(DOSE)              # For disease ontology analysis

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

# Calculate gene variance and select top genes
gene_variance <- apply(norm_counts, 1, var)
variance_threshold <- quantile(gene_variance, 0.75)
top_genes <- names(gene_variance[gene_variance > variance_threshold])
counts_top <- norm_counts[top_genes, ]

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
    pcaData <- plotPCA(vst_data, intgroup = "Group", ntop = 5000, returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    ggplot(pcaData, aes(PC1, PC2, color = Group)) +
        geom_point(size = 3) +
        geom_text(aes(label = rownames(pcaData)), vjust = 1.9, size = 3) +
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
        xlim(c(-100, 150)) +
        ylim(c(-75,75))
    
}

# Create PCA plots
pca_plot_all <- create_pca_plot(vst_data, "PCA Plot (All Groups)", "A")
vst_data_subset <- vst_data[, vst_data$Group != "Normal"]
pca_plot_subset <- create_pca_plot(vst_data_subset, "PCA Plot (Excluding Normal Group)", "B")

# Combine and save PCA plots
combined_pca_plots <- grid.arrange(pca_plot_all, pca_plot_subset, ncol = 2)
ggsave(file.path(output_dir, "PCA_plot.png"), combined_pca_plots, width = 14, height = 10)

# Heatmap Analysis ------------------------------------------------------------------------

# Sample distance heatmap
sample_dist <- dist(t(assay(vst_data)))
sample_clust <- hclust(sample_dist)
heatmap <- pheatmap(as.matrix(sample_dist), clustering_distance_rows = sample_dist, 
                    clustering_distance_cols = sample_dist, main = "Sample Distance Heatmap")
png(file.path(output_dir, "Sample_Heatmap.png"), width = 800, height = 800)
print(heatmap)
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
    res <- as.data.frame(results(dds, contrast = comp, alpha = 0.05))
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
results_df$hgnc_symbol[results_df$hgnc_symbol == ""] <- "Not Available"

# Save results
write.csv(results_df, file.path(output_dir, "DESeq2_results_all_comparisons.csv"), row.names = FALSE)
write.csv(summary_table, file.path(output_dir, "Summary_table.csv"), row.names = FALSE)

# Heatmap of Top Genes --------------------------------------------------------------------

pdf(file.path(output_dir, "Heatmaps_all_comparisons.pdf"))
for (comp in comparisons) {
    comp_label <- paste(comp[2], "vs", comp[3], sep = "_")
    data <- subset(results_df, comparison == comp_label)
    
    # Filter and sort data for the current comparison
    top_100_genes <- head(
        arrange(filter(data, Significance != "Not Significant" & !is.na(log2FoldChange)), 
                desc(abs(log2FoldChange))), 100
    )
    
    # Extract top 20 genes for heatmap
    top_gene_ids <- head(top_100_genes$gene_id, 20)
    matched_rows <- match(top_gene_ids, rownames(norm_counts))
    top_gene_norm_counts <- norm_counts[na.omit(matched_rows), ]
    
    # Generate and save the heatmap
    pheatmap(
        log2(top_gene_norm_counts + 1),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        scale = "row",
        clustering_distance_rows = "euclidean",
        clustering_method = "average",
        show_rownames = TRUE,
        show_colnames = TRUE,
        main = paste("Heatmap -", comp_label)
    )
}
dev.off()

# Volcano Plots ---------------------------------------------------------------------------

png(file.path(output_dir, "Volcano_plots_enhanced_all_comparisons.png"), 
    width = 12, height = 8, units = "in", res = 300)

volcano_plots <- list()
for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    comp_label <- paste(comp[2], "vs", comp[3], sep = "_")
    data <- subset(results_df, comparison == comp_label)
    
    volcano_plot <- EnhancedVolcano(
        data,
        lab = data$hgnc_symbol,
        x = 'log2FoldChange',
        y = 'pvalue',
        xlim = c(-10, 10),
        ylim = c(0, 40),
        title = comp_label,
        pCutoff = 10e-5,
        FCcutoff = 2,
        pointSize = 2.0,
        labSize = 4.0,
        col = c('cornsilk1', 'dodgerblue4', 'orange', 'firebrick4'),
        colAlpha = 1,
        gridlines.major = FALSE,
        gridlines.minor = FALSE,
        drawConnectors = TRUE,
        widthConnectors = 0.75,
        colConnectors = 'black'
    ) +
        labs(tag = LETTERS[i]) +
        theme(plot.tag = element_text(size = 14, face = "bold"))
    
    volcano_plots[[i]] <- volcano_plot
}

combined_volcano_plots <- grid.arrange(grobs = volcano_plots, ncol = 2)
print(combined_volcano_plots)
dev.off()

# GO Enrichment Analysis ------------------------------------------------------------------

all_genes <- unique(results_df$gene_id)
comparisons_for_GO <- list(
    "Triple_negative_vs_HER2_positive",
    "Non_triple_negative_vs_HER2_positive",
    "Triple_negative_vs_Non_triple_negative"
)

for (comparison_name in comparisons_for_GO) {
    filtered_genes <- subset(results_df, comparison == comparison_name & Significance != "Not Significant")
    gene_ids <- filtered_genes$gene_id
    
    enrich_result <- enrichGO(
        gene = gene_ids,
        universe = all_genes,
        OrgDb = "org.Hs.eg.db",
        ont = "BP",
        keyType = "ENSEMBL",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
    )
    
    if (!is.null(enrich_result)) {
        enrich_df <- as.data.frame(enrich_result)
        write.csv(enrich_df, file.path(output_dir, paste0("GO_Enrichment_", comparison_name, ".csv")), row.names = FALSE)
        
        bar_plot <- barplot(enrich_result, showCategory = 20) +
            ggtitle(paste("GO Enrichment for", comparison_name))
        ggsave(file.path(output_dir, paste0("GO_Enrichment_", comparison_name, ".png")), bar_plot, width = 10, height = 6)
    }
}
