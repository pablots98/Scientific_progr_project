# Install and load required libraries
#install.packages("pheatmap")
#install.packages("openxlsx")
#install.packages("dplyr")
library(pheatmap)
library(openxlsx)
library(dplyr)

# Ensure your data is in a data.frame format with rows as genes and columns as samples.
# If your data is in a different format, you'll need to adjust this.

# Read a list of your 125 genes of interest
genes_interes_df = read.table("genes.txt", header = FALSE, sep = "\t")
genes_interes = as.vector(genes_interes_df$V1)  # Convert genes_interest into a vector
head(genes_interes)

# Read your expression data
data = read.xlsx("Merged_data.xlsx")
head(data)

# Filter the expression matrix to only include your genes of interest
data_filtered = data[data$HGNC_GeneSymbol_x %in% genes_interes,]
dim(data_filtered)

# Handle duplicate gene names
# Count the frequency of each name in the column
name_freq = table(data_filtered$HGNC_GeneSymbol_x)
# Find names that appear more than once
rep_names = name_freq[name_freq > 1]
# Print repeated names
print(rep_names)

# Remove 'TNF' to avoid errors
exp_data_filtered = data_filtered %>% filter(HGNC_GeneSymbol_x != 'TNF')
# Recalculate name frequencies
name_freq = table(exp_data_filtered$HGNC_GeneSymbol_x)
# Find names that appear more than once
rep_names = name_freq[name_freq > 1]
# Print repeated names
print(rep_names)

# Set row names and remove unnecessary columns
rownames(exp_data_filtered) = exp_data_filtered$HGNC_GeneSymbol_x
exp_data_filtered = exp_data_filtered[, -c(1:6)]

# Log-transform normalized counts from RNA-seq data
exp_data_filtered.log = log2(exp_data_filtered + 1)
exp_data_filtered.log$sum = rowSums(exp_data_filtered.log)
exp_data_filtered.log = exp_data_filtered.log[exp_data_filtered.log$sum > 0.5,]
exp_data_filtered.log = as.data.frame(t(exp_data_filtered.log))

print(exp_data_filtered.log)

# Create the heatmap
pheatmap(
  exp_data_filtered.log[c(1:48),],
  main = "Heatmap of related genes",
  cluster_cols = TRUE,
  cluster_rows = TRUE
)
