install.packages("data.table")
library(data.table)
library(umap)
library(ggplot2)
library(tidyverse)
library(reshape2)

expression_data <- read.table('C:/Users/manum/Desktop/lung.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata <- phenotype_lung

# Ensure that rownames of expression data are sample IDs
rownames(expression_data) <- expression_data[, 1]
expression_data <- expression_data[, -1]

# Ensure there are no duplicate row names and make them unique if necessary
if(any(duplicated(rownames(expression_data)))) {
  rownames(expression_data) <- make.unique(rownames(expression_data))
}

# Align metadata with expression data based on SUBJID
metadata <- metadata[metadata$SUBJID %in% rownames(expression_data), ]
metadata <- metadata[match(rownames(expression_data), metadata$SUBJID), ]

# Run UMAP on the expression data
umap_result <- umap(expression_data)

# Convert the UMAP result to a data frame
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$SUBJID <- rownames(expression_data)

# Add metadata to the UMAP data frame
umap_df$SEX <- metadata$SEX
umap_df$AGE <- metadata$AGE
umap_df$DTHHRDY <- metadata$DTHHRDY

# Plot UMAP result using ggplot2, color by SEX
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = as.factor(SEX))) +
  geom_point(size = 3) +
  labs(title = "UMAP of Gene Expression Data in Lung Tissue", x = "UMAP1", y = "UMAP2") +
  theme_minimal()

# Alternatively, color by AGE
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = AGE)) +
  geom_point(size = 3) +
  labs(title = "UMAP of Gene Expression Data in Lung Tissue", x = "UMAP1", y = "UMAP2") +
  theme_minimal()

# Alternatively, color by DTHHRDY
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = as.factor(DTHHRDY))) +
  geom_point(size = 3) +
  labs(title = "UMAP of Gene Expression Data in Lung Tissue", x = "UMAP1", y = "UMAP2") +
  theme_minimal()


# Calculate the median expression for each gene across all samples
median_expression <- apply(expression_data, 2, median)

# Convert the median expression to a data frame for plotting
median_expression_df <- data.frame(gene = names(median_expression), median_expression = median_expression)

# Plot the histogram of median expression values
ggplot(median_expression_df, aes(x = median_expression)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Median Gene Expression", x = "Median Expression", y = "Frequency") +
  theme_minimal()


# Calculate mean and median of median expressions
mean_median_expression <- mean(median_expression)
median_of_median_expression <- median(median_expression)

# Plot the histogram of median expression values
ggplot(median_expression_df, aes(x = median_expression)) +
  geom_histogram(binwidth = 0.02, fill = "blue", color = "black", alpha = 0.7) +
  geom_density(color = "red", adjust = 1.5, linewidth = 1) +
  geom_vline(aes(xintercept = mean_median_expression), color = "green", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = median_of_median_expression), color = "orange", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Histogram of Median Gene Expression in Lung tissue",
    x = "Median Expression",
    y = "Frequency",
    subtitle = paste("Mean:", round(mean_median_expression, 3), "Median:", round(median_of_median_expression, 3))
  ) +
  theme_minimal()

# Count the occurrences of each gene type
gene_type_counts <- table(finaldatalung$type)
gene_type_df <- as.data.frame(gene_type_counts)
colnames(gene_type_df) <- c("GeneType", "Count")

# Plot the bar plot of gene types
ggplot(gene_type_df, aes(x = GeneType, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Bar Plot of Gene Types", x = "Gene Type", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the bar plot of gene types
ggplot(gene_type_df, aes(x = GeneType, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Bar Plot of Gene Types in lung tissue", x = "Gene Type", y = "Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold")
  ) +
  scale_y_log10() # Use log scale for y-axis

#for the significant genes
# Calculate counts of each gene type
gene_type_counts <- bon_significant_lung %>%
  group_by(type) %>%
  summarize(Count = n())

# Print counts to the console
print(gene_type_counts)

# Plot the histogram with counts
ggplot(data = gene_type_counts, aes(x = type, y = Count)) +
  geom_bar(stat = "identity", fill = "dodgerblue", color = "black") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
  theme_minimal() +
  labs(title = "Distribution of Gene Types in Regression Results",
       x = "Gene Type",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#for the significant genes
# Calculate counts of each gene type
gene_type_count <- bon_significant_muscle %>%
  group_by(type) %>%
  summarize(Count = n())

# Print counts to the console
print(gene_type_count)

# Plot the histogram with counts
ggplot(data = gene_type_count, aes(x = type, y = Count)) +
  geom_bar(stat = "identity", fill = "dodgerblue", color = "black") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
  theme_minimal() +
  labs(title = "Distribution of Gene Types in Regression Results",
       x = "Gene Type",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#for the significant genes
# Calculate counts of each gene type
gene_type_c <- Bon_significant_skin %>%
  group_by(type) %>%
  summarize(Count = n())

# Print counts to the console
print(gene_type_c)

# Plot the histogram with counts
ggplot(data = gene_type_c, aes(x = type, y = Count)) +
  geom_bar(stat = "identity", fill = "dodgerblue", color = "black") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
  theme_minimal() +
  labs(title = "Distribution of Gene Types in Regression Results",
       x = "Gene Type",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
