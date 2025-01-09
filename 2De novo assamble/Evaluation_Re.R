############################################################
## Topic: Evaluate biological replicates                  ##
## Author: Olga Andrea Hernandez Miranda, Miranda H       ##
## Date: 01/05/2024                                       ##
## Note: Various analyses to evaluate data with DeSeq2    ##
############################################################

# Set working directory
directory <- "C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutoral 2/TranscriptomaE/Prueba tabla ED etapas/EvaluarReplicas"
setwd(directory)

library(DESeq2)

# Read the counts file
counts_matrix <- read.table("trinity.isoform.counts.matrix", header = TRUE, row.names = 1)
counts_matrix

# Read the metadata file
sample_metadata <- read.csv("sample_metadata.csv", row.names = 1)
sample_metadata

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix), 
                              colData = sample_metadata, 
                              design = ~ condition)

# Convert the 'condition' column to a factor
sample_metadata$condition <- as.factor(sample_metadata$condition)

# Create the DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(counts_matrix), 
                              colData = sample_metadata, 
                              design = ~ condition)

# Normalization and analysis
dds <- DESeq(dds)

# Visualize dispersion
plotDispEsts(dds)

# 1. Evaluate replicate similarity with PCA
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

# 2. Distance matrix
sampleDists <- dist(t(assay(vsd)))
heatmap(as.matrix(sampleDists), symm = TRUE)

# 3. Correlation matrix
correlation_matrix <- cor(assay(vsd))
heatmap(correlation_matrix, symm = TRUE)

# 4. Pearson correlation between the first two replicates
correlation_pearson <- cor(assay(vsd)[, 1], assay(vsd)[, 2], method = "pearson")
cat("Pearson correlation between replicates 1 and 2: ", correlation_pearson, "\n")

# 5. Spearman correlation between the first two replicates
correlation_spearman <- cor(assay(vsd)[, 1], assay(vsd)[, 2], method = "spearman")
cat("Spearman correlation between replicates 1 and 2: ", correlation_spearman, "\n")

# 6. Pearson correlation across all samples
correlation_matrix_pearson <- cor(assay(vsd), method = "pearson")
cat("Pearson correlation matrix: \n")
print(correlation_matrix_pearson)

#install.packages("corrplot")

# Plot the correlation matrix with viridis
library(corrplot)
library(viridis)

corrplot(correlation_matrix_pearson, method = "color", type = "upper", 
         col = viridis(200),           # Apply viridis palette
         tl.col = "black", tl.srt = 45, # Label color and orientation
         addCoef.col = "black",        # Display numerical values in black
         number.cex = 0.7,             # Number size
         mar = c(0, 0, 1, 0))

# 7. Spearman correlation across all samples
correlation_matrix_spearman <- cor(assay(vsd), method = "spearman")
cat("Spearman correlation matrix: \n")
print(correlation_matrix_spearman)
