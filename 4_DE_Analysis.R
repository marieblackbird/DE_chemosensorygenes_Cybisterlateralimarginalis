# Load necessary libraries
library(DESeq2)
library(tximport)
library(tximportData)
library(readr)
library(pheatmap)
library(openxlsx)
library(vegan)
library(ggplot2)
library(RColorBrewer)
source("Scripts/functions.R")

# Set the working directory
setwd("transcriptome_filtre_300_cleaned")

################################################################################ 
#################################### Import files ##############################
################################################################################ 

# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto

# List of h5 file names to import
file_names <- c(
  "larvae.a1_Clat_Ant.quantif/abundance.h5" ,
  "larvae.a6_Clat_Ant.quantif/abundance.h5" ,
  "larvae.a11_Clat_Ant.quantif/abundance.h5" ,
  "larvae.a5_Clat_legs.quantif/abundance.h5" ,
  "larvae.a10_Clat_legs.quantif/abundance.h5" ,
  "larvae.a15_Clat_legs.quantif/abundance.h5" ,
  "larvae.a4_Clat_Ma.quantif/abundance.h5" ,
  "larvae.a9_Clat_Ma.quantif/abundance.h5" ,
  "larvae.a14_Clat_Ma.quantif/abundance.h5" ,
  "larvae.a3_Clat_PL.quantif/abundance.h5" ,
  "larvae.a8_Clat_PL.quantif/abundance.h5" ,
  "larvae.a13_Clat_PL.quantif/abundance.h5" ,
  "larvae.a2_Clat_PM.quantif/abundance.h5" ,
  "larvae.a7_Clat_PM.quantif/abundance.h5" ,
  "larvae.a12_Clat_PM.quantif/abundance.h5" ,
  "adult.ClatAnt1_1.quantif/abundance.h5" ,
  "adult.ClatAnt2_2.quantif/abundance.h5" ,
  "adult.ClatAnt3_3.quantif/abundance.h5" ,
  "adult.ClatPL1_7.quantif/abundance.h5" ,
  "adult.ClatPL2_8.quantif/abundance.h5" ,
  "adult.ClatPL3_9.quantif/abundance.h5" ,
  "adult.ClatPm1_4.quantif/abundance.h5" ,
  "adult.ClatPm2_5.quantif/abundance.h5" ,
  "adult.ClatPm3_6.quantif/abundance.h5" ,
  "adult.ClatProF1_13.quantif/abundance.h5" ,
  "adult.ClatProF2_14.quantif/abundance.h5" ,
  "adult.ClatProF3_15.quantif/abundance.h5" ,
  "adult.ClatProM1_10.quantif/abundance.h5" ,
  "adult.ClatProM2_11.quantif/abundance.h5" ,
  "adult.ClatProM3_12.quantif/abundance.h5" ,
  "adult.ClatMus_16.quantif/abundance.h5" 
)

# Define sample names
names(file_names) <- c("larvae_ant_rep1", "larvae_ant_rep2", "larvae_ant_rep3", 
                       "larvae_legs_rep1", "larvae_legs_rep2", "larvae_legs_rep3",
                       "larvae_ma_rep1", "larvae_ma_rep2", "larvae_ma_rep3",
                       "larvae_pl_rep1", "larvae_pl_rep2", "larvae_pl_rep3",
                       "larvae_pm_rep1", "larvae_pm_rep2", "larvae_pm_rep3",
                       "adult_ant_rep1", "adult_ant_rep2", "adult_ant_rep3",
                       "adult_pl_rep1", "adult_pl_rep2", "adult_pl_rep3",
                       "adult_pm_rep1", "adult_pm_rep2", "adult_pm_rep3",
                       "adult_prof_rep1", "adult_prof_rep2", "adult_prof_rep3",
                       "adult_prom_rep1", "adult_prom_rep2", "adult_prom_rep3",
                       "adult_mus_rep1")

# Import Kallisto data using tximport
txi.kallisto <- tximport(file_names, type = "kallisto", txOut = TRUE)

# Load the sample information file
samples <- read.table("sample_info.txt",
                      header = TRUE, row.names = 1, sep = "\t")

# Convert 'organ' column to a factor
samples$organ <- factor(samples$organ)

# Convert 'stade' column to a factor
samples$stade <- factor(samples$stade)

# Create a new factor combining 'organ' and 'stade' columns, separated by an underscore
samples$organe_stade <- factor(paste(samples$organ, samples$stade, sep = "_"))

# Create DESeqDataSetFromTximport object 
dds <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = ~ stade + organ)

################################################################################ 
################################### DESeq ######################################    
################################################################################ 

# Perform differential expression analysis with DESeq2 (DESeqDataSet)
dds <- DESeq(dds)

# Save the 'dds' object since the analysis can take a long time to complete
saveRDS(dds, "dds_sup300_wo_filtre_cleaned.rds")

# Load the 'dds' object from the saved file
dds <- readRDS("dds_sup300_wo_filtre_cleaned.rds")

################################################################################ 
########################## Retrieve Annotation #################################    
################################################################################ 

# Use the function to extract genes of interest
annot <- extract_interest_genes(dds, "tableau_annotation_OR.txt")

# Access the selected genes
genes_contigs <- annot$select_genes

# Access the ordered list of orthologs
genes_in_order_dds <- unlist(annot$orthology_list_ordered)

################################################################################ 
########################## Display Count Numbers ##############################  
################################################################################

# For the entire transcriptome
counts_matrix <- counts(dds)
total_reads_per_sample <- colSums(counts_matrix)
print(total_reads_per_sample)

# Create a table specifically for ORs
counts_matrix <- counts(dds[genes_contigs,])
matches <- match(rownames(counts_matrix), odorants_receptors$contig_name)
rownames(counts_matrix) <- odorants_receptors$orthology[matches]

################################################################################ 
######################### Normalizations ######################################    
################################################################################

# Several normalizations
dds.vst <- varianceStabilizingTransformation(dds)
dds.vst.chemogenes <- dds.vst[genes_contigs,]

dds.normtransform <- assay(normTransform(dds))
dds.normtransform.chemogenes <- dds.normtransform[genes_contigs,]

################################################################################ 
######################## Heatmap ###############################################    
################################################################################

ordre_colonnes <- c(
  "adult_ant_rep1", "adult_ant_rep2", "adult_ant_rep3",
  "adult_pm_rep1", "adult_pm_rep2", "adult_pm_rep3",
  "adult_pl_rep1", "adult_pl_rep2", "adult_pl_rep3",
  "larvae_ant_rep1", "larvae_ant_rep2", "larvae_ant_rep3",
  "larvae_pm_rep1", "larvae_pm_rep2", "larvae_pm_rep3",
  "larvae_pl_rep1", "larvae_pl_rep2", "larvae_pl_rep3",
  "adult_mus_rep1"
)
# Reorganize the data table according to the column order
dds.normtransform.chemogenes <- dds.normtransform.chemogenes[, ordre_colonnes]

drows = dist(dds.normtransform.chemogenes[,c(1:9)] , method = "minkowski") 

# ma_palette_GR = colorRampPalette(c("#4575b4","#4575b4","#4575b4","#4575b4", rev(brewer.pal(n = 7, name = "RdYlBu"))))
# ma_palette_IR = colorRampPalette(c("#4575b4", rev(brewer.pal(n = 7, name ="RdYlBu"))))
# ma_palette_OBP = colorRampPalette(c("#4575b4", rev(brewer.pal(n = 7, name ="RdYlBu"))))

# Display the heatmap
pheatmap(
  dds.normtransform.chemogenes,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  labels_row = genes_in_order_dds,
  fontsize_row = 6,
  border_color = NA,
  clustering_distance_rows = drows,
  color = ma_palette_GR(100)
)

################################################################################ 
################ NMDS (Non-metric Multidimensional Scaling) ####################    
################################################################################

# Extract the relevant data (genes of interest) from the DESeq2 VST-transformed data
table <- assay(dds.vst)[genes_contigs,] 

# Transpose the table so that samples are rows and genes are columns
table_transposed <- t(table) 

# Perform NMDS (Non-metric Multidimensional Scaling) analysis using Euclidean distance
nmds_result <- metaMDS(table_transposed, distance = "euclidean", k=2)

# Create a data frame with NMDS coordinates
nmds_data <- data.frame(
  NMDS1 = nmds_result$points[, 1], # First NMDS axis
  NMDS2 = nmds_result$points[, 2]  # Second NMDS axis
)

# Create the color vector for different conditions
colors <- c(
  rep("red", 3), rep("pink", 3), rep("orange", 3),  # Colors for larvae_ant, larvae_legs, larvae_ma
  rep("yellow", 3), rep("green", 3), rep("lightblue", 3), # Colors for larvae_pl, larvae_pm, adult_ant
  rep("blue", 3), rep("darkblue", 3), rep("purple", 3),  # Colors for adult_pl, adult_pm, adult_prof
  rep("magenta", 3), "darkgrey" # Colors for adult_prom, adult_mus
)

# Define the order of conditions to match the color assignment
condition_order <- c(
  "larvae_ant_rep1", "larvae_ant_rep2", "larvae_ant_rep3",
  "larvae_legs_rep1", "larvae_legs_rep2", "larvae_legs_rep3",
  "larvae_ma_rep1", "larvae_ma_rep2", "larvae_ma_rep3",
  "larvae_pl_rep1", "larvae_pl_rep2", "larvae_pl_rep3",
  "larvae_pm_rep1", "larvae_pm_rep2", "larvae_pm_rep3",
  "adult_ant_rep1", "adult_ant_rep2", "adult_ant_rep3",
  "adult_pl_rep1", "adult_pl_rep2", "adult_pl_rep3",
  "adult_pm_rep1", "adult_pm_rep2", "adult_pm_rep3",
  "adult_prof_rep1", "adult_prof_rep2", "adult_prof_rep3",
  "adult_prom_rep1", "adult_prom_rep2", "adult_prom_rep3",
  "adult_mus_rep1"
)

# Assign the conditions to the NMDS data frame
nmds_data$Condition <- condition_order

# Reorganize the levels of the Condition factor to match the specified order
nmds_data$Condition <- factor(nmds_data$Condition, levels = condition_order)

# Create the NMDS plot using ggplot2
ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, color = Condition)) +
  geom_point() +  # Plot the points
  scale_color_manual(values = colors) +  # Apply the color scheme
  labs(title = "NMDS Plot") +  # Add a title to the plot
  theme_minimal()  # Use a minimal theme for the plot

# Stress plot (Shepard diagram) to assess the quality of NMDS model fit
stressplot(nmds_result, p.col = "black", cex = 0.5)

# Add stress value to the plot for reference
mtext(paste("Stress: ", nmds_result$stress), 
      side = 1, line = 4, cex = 0.8, 
      col = "darkgray") 


################################################################################ 
######################## DESeq2 by organe_stade ################################    
################################################################################

# Create a new variable in the design to represent the organ and stage
dds2 = dds

# Combine 'organ' and 'stade' into a single factor variable
dds2$organe_stade <- factor(paste(samples$organ, samples$stade, sep = "_"))

# Update the design formula to use the new 'organe_stade' variable
design(dds2) = ~ organe_stade

# Run the DESeq2 analysis (this step may take some time)
dds2 = DESeq(dds2)

# Check the names of the results available after the analysis
resultsNames(dds2)

# Save the DESeq2 analysis object to avoid re-running the previous step
saveRDS(dds2, "dds2.rds")

# Load the saved DESeq2 object if it has already been computed
dds2 <- readRDS("dds2.rds") 

##################### Collapse replicates #################################

# Collapse replicates based on the 'organe_stade' variable
ddsColl <- collapseReplicates(dds2, dds2$organe_stade)

# normTransform provides log2(n + 1) normalization, and assay extracts the matrix
tablecoll <- assay(normTransform(ddsColl))[genes_contigs,]

# Define the order of organs (according the order of segments/metameres)
ordre_colonnes <- c("ant_adult", "pm_adult", 
                    "pl_adult", "pro_adult",
                    "ant_larvae", "ma_larvae", 
                    "pm_larvae", "pl_larvae", 
                    "pro_larvae", "mus_adult"
)  

# Reorganize the data table according to the specified column order
tablecoll_reorganise <- tablecoll[, ordre_colonnes]

# Generate the heatmap (supplementary data)
pheatmap(
  tablecoll_reorganise,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  labels_row = genes_in_order_dds,
  fontsize_row = 6,
  color = ma_palette_GR(100)
  border_color = NA,
  clustering_distance_rows = drows,
)

################################################################################ 
######################## Clustering ############################################    
################################################################################

# Collapse the replicates based on the 'organe_stade' variable
ddsColl <- collapseReplicates(dds2, dds2$organe_stade)

# normTransform gives log2(n + 1) transformation, and assay extracts the matrix
tablecoll <- assay(normTransform(ddsColl))[genes_contigs,]

# Apply Variance Stabilizing Transformation (VST) to the collapsed data
dds.vst <- varianceStabilizingTransformation(ddsColl)

# Calculate sample-to-sample distances using the VST-transformed data for selected genes
sampleDists <- dist(t(assay(dds.vst)[genes_contigs,]))

# Convert the distance object to a matrix for visualization
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ddsColl$organe_stade)
colnames(sampleDistMatrix) <- paste(ddsColl$organe_stade)

# Perform hierarchical clustering on the sample distances
tree_clustering <- hclust(sampleDists)

# Plot the hierarchical clustering as a dendrogram
plot(tree_clustering)

################################################################################
############################# Statistics #######################################
################################################################################

######################## Larva Organ by Organ Comparison #######################

# Perform differential expression analysis: Compare ant larvae vs. pm larvae
res_larva_ant_vs_pm=results(dds2, contrast = c("organe_stade", "ant_larvae", "pm_larvae"))
res_larva_ant_vs_pm <- data.frame(genes_in_order_dds , rownames(res_larva_ant_vs_pm[genes_contigs, ]), res_larva_ant_vs_pm[genes_contigs, ])
res_larva_ant_vs_pm$pval_corrigees_BH <- p.adjust(res_larva_ant_vs_pm[, 7], method = "BH")
write.xlsx(res_larva_ant_vs_pm, "res_larva_ant_vs_pm.xlsx", colNames = TRUE, rownames = TRUE)

# Perform differential expression analysis: Compare ant larvae vs. pl larvae
res_larva_ant_vs_pl=results(dds2, contrast = c("organe_stade", "ant_larvae", "pl_larvae"))
res_larva_ant_vs_pl <- data.frame(genes_in_order_dds , rownames(res_larva_ant_vs_pl[genes_contigs, ]), res_larva_ant_vs_pl[genes_contigs, ])
res_larva_ant_vs_pl$pval_corrigees_BH <- p.adjust(res_larva_ant_vs_pl[, 7], method = "BH")
write.xlsx(res_larva_ant_vs_pl, "res_larva_ant_vs_pl.xlsx", colNames = TRUE, rownames = TRUE)

# Perform differential expression analysis: Compare pm larvae vs. pl larvae
res_larva_pm_vs_pl=results(dds2, contrast = c("organe_stade", "pm_larvae", "pl_larvae"))
res_larva_pm_vs_pl <- data.frame(genes_in_order_dds , rownames(res_larva_pm_vs_pl[genes_contigs, ]), res_larva_pm_vs_pl[genes_contigs, ])
res_larva_pm_vs_pl$pval_corrigees_BH <- p.adjust(res_larva_pm_vs_pl[, 7], method = "BH")
write.xlsx(res_larva_pm_vs_pl, "res_larva_pm_vs_pl.xlsx", colNames = TRUE, rownames = TRUE)

######################## Adult Organ by Organ Comparison #######################

# Perform differential expression analysis: Compare ant adult vs. pm adult
res_adult_ant_vs_pm=results(dds2, contrast = c("organe_stade", "ant_adult", "pm_adult"))
res_adult_ant_vs_pm <- data.frame(genes_in_order_dds , rownames(res_adult_ant_vs_pm[genes_contigs, ]), res_adult_ant_vs_pm[genes_contigs, ])
res_adult_ant_vs_pm$pval_corrigees_BH <- p.adjust(res_adult_ant_vs_pm[, 7], method = "BH")
write.xlsx(res_adult_ant_vs_pm, "res_adult_ant_vs_pm.xlsx", colNames = TRUE, rownames = TRUE)

# Perform differential expression analysis: Compare ant adult vs. pl adult
res_adult_ant_vs_pl=results(dds2, contrast = c("organe_stade", "ant_adult", "pl_adult"))
res_adult_ant_vs_pl <- data.frame(genes_in_order_dds , rownames(res_adult_ant_vs_pl[genes_contigs, ]), res_adult_ant_vs_pl[genes_contigs, ])
res_adult_ant_vs_pl$pval_corrigees_BH <- p.adjust(res_adult_ant_vs_pl[, 7], method = "BH")
write.xlsx(res_adult_ant_vs_pl, "res_adult_ant_vs_pl.xlsx", colNames = TRUE, rownames = TRUE)

# Perform differential expression analysis: Compare pm adult vs. pl adult
res_adult_pm_vs_pl=results(dds2, contrast = c("organe_stade", "pm_adult", "pl_adult"))
res_adult_pm_vs_pl <- data.frame(genes_in_order_dds , rownames(res_adult_pm_vs_pl[genes_contigs, ]), res_adult_pm_vs_pl[genes_contigs, ])
res_adult_pm_vs_pl$pval_corrigees_BH <- p.adjust(res_adult_pm_vs_pl[, 7], method = "BH")
write.xlsx(res_adult_pm_vs_pl, "res_adult_pm_vs_pl.xlsx", colNames = TRUE, rownames = TRUE)
