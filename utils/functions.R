################################################################################ 
######## Function to extract genes of interest from an annotation table ########
################################################################################ 

# This function extracts genes of interest based on a given annotation file.
# It matches the contig names in the DESeqDataSet with those in the annotation file
# and returns a list containing selected gene contigs and their corresponding orthology.

extract_interest_genes <- function(dds, annotation_file) {
  
  # Read the annotation table from the specified file
  annotation <- read.csv(annotation_file, header = TRUE, sep = "\t")
  
  # Retrieve contig names that are present both in DESeqDataSet and the annotation file
  select_genes <- rownames(dds)[!is.na(rownames(dds)) & rownames(dds) %in% annotation$contig_name]
  
  # Subset the annotation table to match the order of genes in the DESeqDataSet
  genes_in_order_dds <- annotation[match(select_genes, annotation$contig_name), ]
  
  # Extract the "orthology" column as a list for later use
  orthology_list_ordered <- as.list(genes_in_order_dds$orthology)
  
  # Return a list containing selected gene contigs and the ordered orthology list
  return(list(select_genes = select_genes, orthology_list_ordered = orthology_list_ordered))
}

######### How to use

# Example usage of the extract_interest_genes function
# annot <- extract_interest_genes(dds, "/path/to/tableau_annotation_OBP.txt")

# Accessing the selected genes
# genes_contigs <- annot$select_genes

# Accessing the ordered orthology list
# genes_in_order_dds <- unlist(annot$orthology_list_ordered)


################################################################################ 
####################### Function to perform DE analysis ########################
################################################################################ 

# This function performs differential expression analysis between two specified levels
# for a given factor (context) in the DESeqDataSet. It also applies p-value adjustment
# and writes the results to an Excel file.

perform_differential_analysis <- function(dds2, context, level1, level2, output_file) {
  
  # Perform differential expression analysis using the specified contrast
  contrast <- c(context, level1, level2)
  results_df <- results(dds2, contrast = contrast, test = "Wald")
  
  # Subset the results for genes of interest
  genes_contigs <- rownames(results_df) %in% genes_in_order_dds
  results_subset <- data.frame(genes_in_order_dds, rownames(results_df[genes_contigs, ]), results_df[genes_contigs, ])
  
  # Adjust p-values using Benjamini-Hochberg (BH) correction
  results_subset$pval_corrigees_BH <- p.adjust(results_subset[, "padj"], method = "BH")
  
  # Write the results to an Excel file
  write.xlsx(results_subset, output_file, colNames = TRUE, rownames = TRUE)
  
  # Return the results dataframe for further inspection
  return(results_subset)
}

######### How to use

# Example usage of the perform_differential_analysis function
# result <- perform_differential_analysis(dds2, "organe_stade", "pl_larvae", "pm_larvae", "res_pl_larvae_VS_pm_larvae.xlsx")
