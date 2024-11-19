library(readr)
phenotype<-read.delim("C:/Users/manum/Desktop/lung.txt", header=TRUE, sep = "\t")

#install libraries
library(dplyr)
library(tidyr)
library(stringr)

median_expression <- apply(phenotype, 2, median)
# Convert the result to a data frame for plotting
median_df <- data.frame(median_expression = median_expression)

library(ggplot2)
ggplot(median_df, aes(x = median_expression))+
  geom_histogram(binwidth = 0.5, fill= "red", color = "black")+
  labs(title = "Median expression in Lung Tissue", x = "Median Expression", y= "frequency")+
  theme_minimal()

#transpose the dataframe#transpose the dataframemedian
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
merge1 <- left_join(data_long, phenotype_lung, by = c("SUBJID"= "SUBJID") )

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
library(writexl)


#single gene linear

# Specify the gene of interest
gene_of_interest <- "CDKN2A"

# Filter the data for the specific gene
specific_gene <- merge1[merge1$name == "CDKN2A", ]

write_xlsx(specific_gene, "age_single_lung.xlsx")

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



#creating plots
#histogram of p value

install.packages("ggplot2")
library(ggplot2)

# Ensure p-values are numeric
age_regression_lung$P_Value <- as.numeric(age_regression_lung$P_Value)

# Create the histogram
age_lung <- ggplot(age_regression_lung, aes(x = P_Value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histogram of P-Values",
       x = "P-Value",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(age_lung)


#volcano plot


# Create a new column for -log10(p-value)
age_regression_results_lung$logP <- -log10(age_regression_results_lung$P_Value)

# Define significance thresholds
significance_threshold <- 0.05
logP_threshold <- -log10(significance_threshold)

# Create the volcano plot
volcano_plot <- ggplot(age_regression_results_lung, aes(x = Estimate, y = logP)) +
  geom_point(aes(color = P_Value < significance_threshold), alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
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
pvalue<-age_regression_results_lung$P_Value

# Calculate q-values from p-values
qobj <- qvalue(p = age_regression_results_lung$P_Value)

# Add the q-values to your data frame
age_regression_results_lung$q_value <- qobj$qvalues

# View the updated data frame
head(age_regression_results_lung)



# Define the q-value threshold for significance
q_value_threshold <- 0.05

# Add a new column to indicate significance
age_regression_results_lung$Significant <- age_regression_results_lung$q_value < q_value_threshold

# View the updated data frame
print(age_regression_results_lung)

# Filter the data frame to include only significant genes
significant_genes <- age_regression_results_lung[age_regression_results_lung$Significant == TRUE, ]

# Optionally, save the significant genes to a new CSV file
write.csv(significant_genes, "significant_genes_qvalues.csv", row.names = FALSE)



#create manhatten plot

#combine the regression with gencode file for chr and position of snp

merge2<- merge(age_regression_results_lung, gencode26, by.x="Gene", by.y= "name")

write.csv(merge2, "age_regression_lung.csv", row.names = FALSE)

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
          col = c("blue4", "orange3"), chrlabs = NULL, main = "Manhattan Plot-lung tissue",
          ylim = c(0, max(-log10(merge2$P_Value) + 2)))



#subset of females and males

females <- subset(finaldatalung, SEX == "2")
males <- subset(finaldatalung, SEX == "1")


#separate regression for males and females
library(qvalue)

# Function to run regression and calculate q-values
run_regression <- function(males) {
  # Create an empty dataframe to store results
  results <- data.frame(
    Gene = character(),
    Estimate = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Perform linear regression for each gene
  for (gene in unique(males$name)) {
    # Filter data for the current gene
    gene_data <- subset(males, name == gene)
    
    # Perform linear regression: Expression ~ AGE
    model <- lm(Expression ~ AGE, data = gene_data)
    
    # Extract p-value and estimate for AGE
    estimate <- coef(summary(model))["AGE", "Estimate"]
    p_value <- coef(summary(model))["AGE", "Pr(>|t|)"]
    
    # Store the results in the dataframe
    results <- rbind(results, data.frame(Gene = gene, Estimate = estimate, P_Value = p_value))
  }
  
  # Calculate q-values from p-values
  qobj <- qvalue(p = results$P_Value)
  results$q_value <- qobj$qvalues
  
  return(results)
}

# Run regression for male and female data
male_results <- run_regression(males)
female_results <- run_regression(females)


# Define the q-value threshold for significance
q_value_threshold <- 0.05

# Add a significance column
male_results$Significant <- male_results$q_value < q_value_threshold
female_results$Significant <- female_results$q_value < q_value_threshold

# View the results
print(male_results)
print(female_results)

# Optionally, save the results to CSV files
write.csv(male_results, "male_regression_results_lung.csv", row.names = FALSE)
write.csv(female_results, "female_regression_results_lung.csv", row.names = FALSE)


# p value adjustment using the bonferroni correction
Bon_threshold<- 0.05/ length (merge2$P_Value)

bon_significant <- subset(merge2, P_Value<Bon_threshold)
write.csv(bon_significant, "bon_significant_lung.csv", row.names = FALSE)


write.csv(gencode26, "genecode.csv", row.names = FALSE)



#for the pathway analysis,
# Install and load necessary packages
if (!requireNamespace("biomaRt", quietly = TRUE))
  install.packages("biomaRt")
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

library(biomaRt)
library(tidyverse)

#2nd approach

# Extract the gene lists
significant_genes_list <- bon_significant_lung$Gene
gwas_genes_list <- unique(unlist(strsplit(EFO_0000341_associations_export$mappedGenes, ",")))

# Ensure the gene names are trimmed of any extra spaces
gwas_genes_list <- trimws(gwas_genes_list)

# Define the sets
significant_genes_list <- unique(significant_genes_list)
gwas_genes_list <- unique(gwas_genes_list)

# Calculate A, B, C, D
A <- length(intersect(significant_genes_list, gwas_genes_list))
B <- length(setdiff(gwas_genes_list, significant_genes_list))
C <- length(setdiff(significant_genes_list, gwas_genes_list))
total_genes <- unique(c(significant_genes_list, gwas_genes_list))
D <- length(total_genes) - A - B - C

# Create the contingency table
contingency_table <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE,
                            dimnames = list("Differentially Expressed" = c("Yes", "No"),
                                            "In GWAS List" = c("Yes", "No")))

print(contingency_table)

# Perform Fisher's Exact Test
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)

# Optional: Save the table to a CSV file
write.csv(contingency_table, "contingency_table_lung.csv")


#2nd approach

# Extract the gene lists
total_gene_list <- age_regression_lung$Gene
significant_genes_list <- bon_significant_lung$Gene
gwas_association_data <- EFO_0000341_associations_export
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
# Calculate odds ratio and log odds ratio
odds_A <- A / B
odds_B <- C / D
odds_ratio <- odds_A / odds_B
log_odds_ratio <- log(odds_ratio)

print(contingency_table)

# Perform Fisher's Exact Test
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)
p_value <- fisher_test_result$p.value

# Create a data frame for plotting
plot_data <- data.frame(
  Category = c("In DE and In GWAS", "In DE and Not in GWAS", "Not in DE and In GWAS", "Not in DE and Not in GWAS"),
  Count = c(A, B, C, D)
)

# Calculate expected values
total_genes_count <- length(total_genes)
expected_A <- (length(significant_genes_list) * length(common_genes)) / total_genes_count
expected_B <- length(significant_genes_list) - expected_A
expected_C <- length(common_genes) - expected_A
expected_D <- total_genes_count - expected_A - expected_B - expected_C

# Create a data frame for plotting
plot_data <- data.frame(
  Category = c("In DE and In GWAS", "In DE and Not in GWAS", "Not in DE and In GWAS", "Not in DE and Not in GWAS"),
  Observed = c(A, B, C, D),
  Expected = c(expected_A, expected_B, expected_C, expected_D)
)

# Create the bar plot
library(ggplot2)
p <- ggplot(plot_data, aes(x = Category)) +
  geom_bar(aes(y = Observed, fill = "Observed"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = Expected, fill = "Expected"), stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Overlap Between Differentially Expressed Genes and GWAS Genes",
       x = "Category",
       y = "Gene Count") +
  scale_fill_manual(values = c("dodgerblue", "goldenrod1", "forestgreen", "gray", "red", "blue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

# Create a data frame for plotting
plotdata <- data.frame(
  Category = rep(c("In DE and In GWAS", "In DE and Not in GWAS", "Not in DE and In GWAS", "Not in DE and Not in GWAS"), each = 2),
  Count = c(A, B, C, D, expected_A, expected_B, expected_C, expected_D),
  Type = rep(c("Observed", "Expected"), times = 4)
)

# Create the bar plot
library(ggplot2)
q <- ggplot(plotdata, aes(x = Category, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Overlap Between Differentially Expressed Genes and GWAS Genes",
       x = "Category",
       y = "Gene Count") +
  scale_fill_manual(values = c("dodgerblue", "goldenrod1")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add statistical annotation
q + geom_text(aes(label = round(Count, 1)), position = position_dodge(width = 0.9), vjust = -0.25)




# Prepare data for plotting
data <- data.frame(
  log_odds_ratio = log_odds_ratio,
  neg_log10_pvalue = -log10(p_value)
)

# Create the plot
ggplot(data, aes(x = log_odds_ratio, y = neg_log10_pvalue)) +
  geom_point() +
  labs(title = "Log Odds Ratio vs -log10(p-value)",
       x = "Log Odds Ratio",
       y = "-log10(p-value)") +
  theme_minimal()


# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Import gene list from file (assuming CSV format)
# Replace 'converted_gene_ids.csv' with the path to your file

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
barplot(kegg_enrichment, showCategory = 20, title = "Top 20 Enriched KEGG Pathways")

# Create a dotplot for visualization
dotplot(kegg_enrichment, showCategory = 20, title = "Top 20 Enriched KEGG Pathways")

# Convert results to a data frame for further inspection
kegg_enrichment_df <- as.data.frame(kegg_enrichment)

write.csv(kegg_enrichment_df, "kegg-lung.csv", row.names = FALSE)

# Optionally, print the data frame to inspect the results
print(kegg_enrichment_df)

significance_threshold <- 0.002
# Filter for significantly differentially expressed genes
pathways <- kegg_enrichment_df[kegg_enrichment_df$p.adjust < significance_threshold, ]



#######################################################
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)

pheno<- phenotype_skin
pheno$SEX <- factor(pheno$SEX, levels = c(1,2), labels = c("Male", "Female"))

#age and sex
ggplot(data = pheno, aes(x = as.factor(SEX), y = AGE))+
  geom_violin( fill = "blue", color = "black", alpha = 0.7)+
  labs (title = "Age Distribution by sex for Muscle tissue", x = "Sex", y = "Age")+
  theme_minimal()

#sex count
ggplot(pheno, aes(x= SEX))+
  geom_bar(fill = c("lightblue", "pink"))+
  labs(title = "Counts of Individuals by Sex for Muscle tissue", x = "Sex", y = "count")+
  theme_minimal()

#Convert DTHHRDY to a factor with labels
pheno$DTHHRDY <- factor(pheno$DTHHRDY, levels = c(0, 1, 2, 3, 4), 
                        labels = c("Ventilator Case", "Violent and fast death", "Fast death of natural causes", 
                                   "Intermediate death", "Slow death"))

# Create a bar plot for DTHHRDY count
ggplot(pheno, aes(x = DTHHRDY)) +
  geom_bar(fill = "lightblue", color = "darkblue") +
  labs(title = "Count of DTHHRDY Categories for Muscle", x = "Death Category", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text for better readability
##############################################################################################################

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
total_DEGs <- 1223
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
  labs(title = "Top 20 KEGG Pathway Enrichment in lung tissue", x = "Log Odds Ratio", y = "Pathway", fill = "-Log10(p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


