library(readr)
phenotype<-read.delim("C:/Users/manum/Desktop/muscle_skeletal.txt", header=TRUE, sep = "\t")

#install libraries
library(dplyr)
library(tidyr)
library(stringr)

#transpose the dataframe
data<- t(phenotype)


class(data)
if(! is.data.frame(data)){
  data<-as.data.frame(data)
}


#match the gene data with the expression data for the gene information
#give column name for the expression id
library(tibble)
data1<- data %>% 
  rownames_to_column(var="id")


# Pivot the data to long format
data_long <- data1 %>%
  pivot_longer(
    cols = -id,        # Exclude the 'Sample' column from being pivoted
    names_to = "SUBJID",     # Name for the new 'key' column
    values_to = "Expression"  # Name for the new 'value' column
  )


#merge data long with phenotype skin table
merge1 <- left_join(data_long, phenotype_skin, by = c("SUBJID"= "SUBJID") )

#merge the table with gencode26 table
merge1<- left_join(merge1, gencode26, by = c("id" = "id"))


#for the regression 
# Make sure the 'expression' column is numeric
merge1$Expression <- as.numeric(merge1$Expression)

# Create an empty list to store regression results
regression_results <- list()

# Perform linear regression for each gene
for(gene in unique(merge1$name)) {
  # Filter the data for the current gene
  gene_data <- merge1 %>% filter(name == gene)
  
  # Perform linear regression: Value ~ AGE
  model <- lm(Expression ~ AGE, data = gene_data)
  
  # Store the model summary
  regression_results[[gene]] <- summary(model)
}

# Create an empty dataframe to store p-values and estimates
results <- data.frame(
  Gene = character(),
  Estimate = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through the stored regression results
for(gene in names(regression_results)) {
  # Get the summary of the model
  summary_model <- regression_results[[gene]]
  
  # Extract the p-value and estimate for AGE
  if ("AGE" %in% rownames(summary_model$coefficients)) {
    estimate <- summary_model$coefficients["AGE", "Estimate"]
    p_value <- summary_model$coefficients["AGE", "Pr(>|t|)"]
  } else {
    estimate <- NA
    p_value <- NA
  }
  
  # Store the results in the dataframe
  results <- rbind(results, data.frame(Gene = gene, Estimate = estimate, P_Value = p_value))
}

# Print the results dataframe
print(results)

# Save the results as a CSV file
write.csv(results, "age_regression_results.csv", row.names = FALSE)







#find the significant p values
# Define the significance threshold
significance_threshold <- 0.05

# Filter for significantly differentially expressed genes
significant_genes <- results[results$P_Value < significance_threshold, ]

# Display the significant genes
print(significant_genes)

#highest p value with gen
minpvalue<- min(significant_genes$P_Value)
genewithpvalue<- significant_genes[significant_genes$P_Value == minpvalue,]


#save the results in the excel file

#significant genes table
write_xlsx(significant_genes, "age_significantgene_skin.xlsx")

#single gene linear

# Specify the gene of interest
gene_of_interest <- "LHFPL4"

# Filter the data for the specific gene
specific_gene <- merge1[merge1$name == "LHFPL4", ]

write_xlsx(specific_gene, "age_single_skin.xlsx")

# Perform linear regression: Value ~ AGE
model1 <- lm(Expression ~ AGE, data = specific_gene)

# Get the summary of the model
summary_model1 <- summary(model1)

# Extract the p-value and estimate for SEX
if ("AGE" %in% rownames(summary_model1$coefficients)) {
  estimate <- summary_model1$coefficients["AGE", "Estimate"]
  p_value <- summary_model1$coefficients["AGE", "Pr(>|t|)"]
} else if ("AGE" %in% rownames(summary_model1$coefficients)) { # In case it's named "AGE" directly
  estimate <- summary_model1$coefficients["AGE", "Estimate"]
  p_value <- summary_model1$coefficients["AGE", "Pr(>|t|)"]
} else {
  estimate <- NA
  p_value <- NA
}

# Print the results for the specific gene
cat("Gene:", gene_of_interest, "\n")
cat("Estimate for SEX:", estimate, "\n")
cat("P-Value for SEX:", p_value, "\n")


install.packages("tidyverse")
library(tidyverse)


#creating plots
#histogram of p value

install.packages("ggplot2")
library(ggplot2)

# Ensure p-values are numeric
age_regression_results_skeletal$P_Value <- as.numeric(age_regression_results_skeletal$P_Value)

# Create the histogram
age_skeletal <- ggplot(age_regression_results_skeletal, aes(x = P_Value)) +
  geom_histogram(binwidth = 0.05, fill = "red", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histogram of P-Values - muscle tissue",
       x = "P-Value",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(age_skeletal)


#volcano plot


# Create a new column for -log10(p-value)
age_regression_results_skeletal$logP <- -log10(age_regression_results_skeletal$P_Value)

# Define significance thresholds
significance_threshold <- 0.05
logP_threshold <- -log10(significance_threshold)
library(ggplot2)
# Create the volcano plot
volcano_plot <- ggplot(age_regression_results_skeletal, aes(x = Estimate, y = logP)) +
  geom_point(aes(color = P_Value < significance_threshold), alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "orange", "FALSE" = "black")) +
  theme_minimal() +
  labs(title = "Volcano Plot of p-values in muscles",
       x = "Estimate",
       y = "-log10(p-value)",
       color = "Significant") +
  geom_hline(yintercept = logP_threshold, linetype = "dashed", color = "blue") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(volcano_plot)


#q values
install.packages("BiocManager")
BiocManager::install("qvalue")
library(qvalue)
pvalue<-age_regression_results_skeletal$P_Value

# Calculate q-values from p-values
qobj <- qvalue(p = age_regression_results_skeletal$P_Value)

# Add the q-values to your data frame
age_regression_results_skeletal$q_value <- qobj$qvalues

# View the updated data frame
head(age_regression_results_skeletal)



# Define the q-value threshold for significance
q_value_threshold <- 0.05

# Add a new column to indicate significance
age_regression_results_skeletal$Significant <- age_regression_results_skeletal$q_value < q_value_threshold

# View the updated data frame
print(age_regression_results_skeletal)

# Filter the data frame to include only significant genes
significant_genes <- age_regression_results_skeletal[age_regression_results_skeletal$Significant == TRUE, ]

# Optionally, save the significant genes to a new CSV file
write.csv(significant_genes, "significant_genes_qvalues_muscle.csv", row.names = FALSE)



#create manhatten plot

#combine the regression with gencode file for chr and position of snp

merge2<- merge(age_regression_results_skeletal, gencode26, by.x="Gene", by.y= "name")

write.csv(merge2, "age_regression_muscle.csv", row.names = FALSE)

install.packages("qqman")
library(qqman)

# Remove 'chr' prefix
merge2$chromosome <- gsub("chr", "", merge2$chromosome)

# Replace 'X', 'Y', 'M' with numeric values
merge2$chromosome[merge2$chromosome == "X"] <- "23"
merge2$chromosome[merge2$chromosome == "Y"] <- "24"
merge2$chromosome[merge2$chromosome == "M"] <- "25"

# Convert the chromosome column to numeric
merge2$chromosome <- as.numeric(merge2$chromosome)



manhattan(merge2, chr = "chromosome", bp = "position", p = "P_Value", snp = "Gene",
          genomewideline = -log10(5e-8), suggestiveline = -log10(1e-5), highlight = NULL,
          col = c("blue4", "red3"), chrlabs = NULL, main = "Manhattan Plot-muscle tissue",
          ylim = c(0, max(-log10(merge2$P_Value) + 2)))


#adjsuting the pvalue to reduce the number of hits
p_value_threshold<- 0.01
significant_hits<- subset(age_regression_results_skeletal, P_Value < p_value_threshold)

bonferroni_threshold <- 0.05 / length(merge2$P_Value)
bon_significant <- subset(merge2, P_Value < bonferroni_threshold)


# Apply Benjamini-Hochberg correction
age_regression_results_skeletal$BH_P_Value <- p.adjust(age_regression_results_skeletal$P_Value, method = "BH")

# Define a BH-corrected p-value threshold
bh_threshold <- 0.05

# Filter the results based on the BH-corrected threshold
significant_hits1 <- subset(age_regression_results_skeletal, BH_P_Value < bh_threshold)

write.csv(bon_significant, "bon_significant_muscle.csv", row.names = FALSE)

library(biomaRt)
library(tidyverse)


#2nd approach

# Extract the gene lists
total_gene_list <- age_regression_muscle$Gene
significant_genes_list <- bon_significant_muscle$Gene
gwas_association_data <- EFO_0005856_associations_export
gwas_genes_list <- unique(unlist(strsplit(gwas_association_data$mappedGenes, ",")))

# Ensure the gene names are trimmed of any extra spaces
total_gene_list <- trimws(total_gene_list)
gwas_genes_list <- trimws(gwas_genes_list)

# Find common genes between total genes and GWAS genes
common_genes <- intersect(total_gene_list, gwas_genes_list)

# Define the sets
significant_genes_list <- unique(significant_genes_list)
common_genes <- unique(common_genes)

# Calculate A, B, C, D
A <- length(intersect(significant_genes_list, common_genes))
B <- length(setdiff(significant_genes_list, common_genes))
C <- length(setdiff(common_genes, significant_genes_list))
total_genes <- unique(total_gene_list)  # Use the original total genes list
D <- length(total_genes) - A - B - C

# Create the contingency table
contingency_table <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE,
                            dimnames = list("Differentially Expressed" = c("Yes", "No"),
                                            "In GWAS List" = c("Yes", "No")))

print(contingency_table)

# Perform Fisher's Exact Test
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)

# Create a data frame for plotting
plot_data <- data.frame(
  Category = c("In DE and In GWAS", "In DE and Not in GWAS", "Not in DE and In GWAS", "Not in DE and Not in GWAS"),
  Count = c(A, B, C, D)
)

# Create the bar plot
library(ggplot2)
p <- ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Overlap Between Differentially Expressed Genes and GWAS Genes",
       x = "Category",
       y = "Gene Count") +
  scale_fill_manual(values = c("dodgerblue", "goldenrod1", "forestgreen", "gray"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))

# Display the bar plot
print(p)


# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Import gene list from file (assuming CSV format)
# Replace 'converted_gene_ids.csv' with the path to your file
gene_list <- muscle_converted

# Check the imported data
head(gene_list)
# Extract Entrez IDs from the imported data
entrez_ids <- gene_list$To

# Perform KEGG enrichment analysis
kegg_enrichment <- enrichKEGG(gene = entrez_ids, 
                              organism = 'hsa', 
                              keyType = 'ncbi-geneid', 
                              pAdjustMethod = "BH", 
                              pvalueCutoff = 0.05)

# Check the dimensions of the enrichment result
dim(kegg_enrichment)

# Plot top KEGG pathways
# Plot the barplot with smaller fonts
barplot(kegg_enrichment, showCategory = 20, title = "Top 20 Enriched KEGG Pathways in Muscle Tissue") +
  theme(
    plot.title = element_text(size = 12),  # Adjust title font size
    axis.title.x = element_text(size = 10),  # Adjust x-axis title font size
    axis.title.y = element_text(size = 10),  # Adjust y-axis title font size
    axis.text.x = element_text(size = 8),  # Adjust x-axis text font size
    axis.text.y = element_text(size = 8)   # Adjust y-axis text font size
  )

# Create a dotplot for visualization
dotplot(kegg_enrichment, showCategory = 20, title = "Top 20 Enriched KEGG Pathways")

# Convert results to a data frame for further inspection
kegg_enrichment_df <- as.data.frame(kegg_enrichment)

# Optionally, print the data frame to inspect the results
print(kegg_enrichment_df)

#log odds ratio vs -log10(pvalues)

# Function to extract numerator from ratio
extract_numerator <- function(ratio) {
  as.numeric(strsplit(ratio, "/")[[1]][1])
}

# Function to extract denominator from ratio
extract_denominator <- function(ratio) {
  as.numeric(strsplit(ratio, "/")[[1]][2])
}

# Apply the functions to extract numerators and denominators
data <- kegg_enrichment_df %>%
  mutate(
    GeneRatioNum = sapply(GeneRatio, extract_numerator),
    BgRatioNum = sapply(BgRatio, extract_numerator),
    GeneRatioDenom = sapply(GeneRatio, extract_denominator),
    BgRatioDenom = sapply(BgRatio, extract_denominator)
  )

# Total DEGs and background genes
total_DEGs <- 1491
total_background_genes <- 8842

# Calculate components for odds ratio
data <- data %>%
  mutate(
    a = Count,
    b = total_DEGs - Count,
    c = BgRatioNum - GeneRatioNum,
    d = total_background_genes - BgRatioNum,
    odds_ratio = (a * d) / (b * c),
    log_odds_ratio = log(odds_ratio),
    neg_log_pvalue = -log10(pvalue)
  )

# Select the top 20 pathways
top20_data <- data %>% arrange(pvalue) %>% head(20)


# Plotting
ggplot(top20_data, aes(x = log_odds_ratio, y = reorder(Description, log_odds_ratio), fill = neg_log_pvalue)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Top 20 KEGG Pathway Enrichment in Muscle tissue", x = "Log Odds Ratio", y = "Pathway", fill = "-Log10(p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))



