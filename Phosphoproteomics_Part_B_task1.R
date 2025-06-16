# Declaration:
# I, Samuel Rice, declare that I have not used generative AI
# in the creation of this R script.

# Read in the metdata
metadata_df = read.csv("TMTchannelAssignment.csv")

# Check IDs column
metadata_df$nID

# Read in peptide quantification data
df = read.csv("filtered_mod_peptides.csv")

# Check column names
colnames(df)


# Rename columns 15-30 of pepetide dataset with the metadata ids
colnames(df)[15:30] = metadata_df$nID

library(dplyr) # load dplyr for data manipulation 
# Select the relevant columns for analysis
# and the intensity value columns for 36 hour mark
df = select(df, Gene.Names, PEP, Proteins, Modifications, 
            S1.36.A, S1.36.B, S1.36.C, S2.36.A, S2.36.B,S2.36.C,
            mock.36.A, mock.36.B, mock.36.C)
colnames(df)     

# Replace 0 intensity values with 0.1
# this prevents log transformation issues
df[5:13][df[5:13] == 0] <- 0.1 
head(df)


# Load libraries for visualisation and data manipulation
library(ggplot2)
library(tidyr)

# log transform each intensity column to normalise distribution
df_logged <- df %>%
  mutate(
    S1.36.A_log = log(S1.36.A),
    S1.36.B_log = log(S1.36.B),
    S1.36.C_log = log(S1.36.C),
    S2.36.A_log = log(S2.36.A),
    S2.36.B_log = log(S2.36.B),
    S2.36.C_log = log(S2.36.C),
    mock.36.A_log = log(mock.36.A),
    mock.36.B_log = log(mock.36.B),
    mock.36.C_log = log(mock.36.C)
  )

# Reshape the data to long format
# so that it is suitable for boxplots
# gathers multiple values into 2 columns
# column for variables, column for variable values
df_long <- df_logged %>%
  gather(key = "column", value = "value", 
         S1.36.A_log, S1.36.B_log, S1.36.C_log, 
         S2.36.A_log, S2.36.B_log, S2.36.C_log, 
         mock.36.A_log, mock.36.B_log, mock.36.C_log)

# Create a boxplot
# visualise differencesd between log transformed intensity distributions
ggplot(df_long, aes(x = column, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Log-Transformed Columns", 
       x = "Columns", 
       y = "Log-transformed values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

colnames(df)
colnames(df_logged)
# Calculate a reference channel median
# as a baseline for normalisation
median_ref_col <- median(df_logged[,"mock.36.A_log"])

# Loop through columns 13 to 22
# normalise data by subtracting the median difference
# for each log transformed column relative to the refference chanel median 
for (i in 13:22) {
  # Get column name
  new_col_name <- paste("norm_", colnames(df_logged)[i], sep = "")
  
  # Calculate global difference
  global_diff <- median(df_logged[, i]) - median_ref_col
  
  # Create new normalized column
  df_logged[[new_col_name]] <- df_logged[, i] - global_diff
}

# Check first few rows and column names 
head(df_logged)
colnames(df_logged)

# Reshape normalized data for boxplot after being normalised
df_long_norm <- df_logged %>%
  gather(key = "column", value = "value", norm_S1.36.A_log, norm_S1.36.B_log, norm_S1.36.C_log,
         norm_S2.36.A_log, norm_S2.36.B_log, norm_S2.36.C_log,
         norm_mock.36.A_log, norm_mock.36.B_log, norm_mock.36.C_log)
colnames(df_long_norm)

# Create a new boxplot to check normalization
# verify the normalisation effect
ggplot(df_long_norm, aes(x = column, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Normalized Log-Transformed Columns",
       x = "Columns",
       y = "Normalized Log-transformed values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


colnames(df_logged)
# Perform unpaired t-test for each row
# S1 compared to mock treatment
for (i in 1:nrow(df_logged)) {
  ttest_result <- t.test(df_logged[i, 24:26], df_logged[i, 30:32], var.equal = TRUE)
  df_logged[i, "s1_pvalue"] <- ttest_result$p.value
}
# Perform unpaired t-test for each row
# S2 compared to mock treatment
for (i in 1:nrow(df_logged)) {
  ttest_result_s2 <- t.test(df_logged[i, 27:29], df_logged[i, 30:32], var.equal = TRUE)
  df_logged[i, "s2_pvalue"] <- ttest_result_s2$p.value
}
# print the summarys
summary(df_logged$s1_pvalue)
summary(df_logged$s2_pvalue)

# Filter to just keep phosphopeptides and save to csv
df_logged <- df_logged[grep("Phospho", df_logged$Modifications), ]
write.csv(df_logged, "filtered_phosphopeptides.csv", row.names = FALSE)

# Count significant phosphopeptides and unique proteins
# for S1
s1_significant_proteins <- df_logged$Proteins[df_logged$s1_pvalue < 0.05]
s1_split_proteins <- unlist(strsplit(s1_significant_proteins, ";"))
s1_protein_count <- length(unique(s1_split_proteins))

# For s2
s2_significant_proteins <- df_logged$Proteins[df_logged$s2_pvalue < 0.05]
s2_split_proteins <- unlist(strsplit(s2_significant_proteins, ";"))
s2_protein_count <- length(unique(s2_split_proteins))

# Print the results to check everything runs smoothly
print(paste("Significant S1 Phosphopeptides:", length(unique(s1_significant_proteins)), "from", s1_protein_count, "unique proteins"))
print(paste("Significant S2 Phosphopeptides:", length(unique(s2_significant_proteins)), "from", s2_protein_count, "unique proteins"))

# Store significant phosphopeptides
# p vaue < 0.05
s1_significant_phospho <- df_logged[df_logged$s1_pvalue < 0.05, ]
s2_significant_phospho <- df_logged[df_logged$s2_pvalue < 0.05, ]

# Save the filtered data to CSV
write.csv(s1_significant_phospho, "s1_significant_phosphopeptides.csv", row.names = FALSE)
write.csv(s2_significant_phospho, "s2_significant_phosphopeptides.csv", row.names = FALSE)


# Calculate log fold change and -log10 p-value
# across all entrys for the volcano plot
df_logged <- df_logged %>%
  mutate(
    mean_S1 = rowMeans(select(df_logged, norm_S1.36.A_log, norm_S1.36.B_log, norm_S1.36.C_log), na.rm = TRUE),
    mean_Mock = rowMeans(select(df_logged, norm_mock.36.A_log, norm_mock.36.B_log, norm_mock.36.C_log), na.rm = TRUE),
    s1_log_fold_change = mean_S1 - mean_Mock,  # Include all genes
    neg_log10_pvalue = -log10(s1_pvalue),
    significance = case_when(
      s1_pvalue < 0.05 & s1_log_fold_change > 0.25 ~ "Upregulated",
      s1_pvalue < 0.05 & s1_log_fold_change < -0.25 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Generate a volcano Plot using all phosphopeptides
# label points by upregulation/ downregulation for clarity 
# uses y threshold -log10 (0.05)
ggplot(df_logged, aes(x = s1_log_fold_change, y = neg_log10_pvalue, color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "cyan", "Not significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "-log10 p-values vs log fold change for phosphopeptides in human cell line infected with SARS-CoV1",
    x = "Log Fold Change",
    y = "-Log10 (p-value)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Significance threshold
  xlim(-2, 2)  # Limit x-axis for better visualization




















# Task 2 part one
# extracting protein list for DAVID

# 📌 Extract Significant Proteins for S1 and S2 (p < 0.05)
# Take only the first protein from multiple mappings (split by ";")
extract_first_protein <- function(protein_column) {
  sapply(strsplit(as.character(protein_column), ";"), `[`, 1)
}

# Extract for S1
s1_significant_proteins <- df_logged$Proteins[df_logged$s1_pvalue < 0.05]
s1_first_proteins <- unique(extract_first_protein(s1_significant_proteins))

# Extract for S2
s2_significant_proteins <- df_logged$Proteins[df_logged$s2_pvalue < 0.05]
s2_first_proteins <- unique(extract_first_protein(s2_significant_proteins))

# 📌 Extract Background List (All phosphopeptides)
background_proteins <- unique(extract_first_protein(df_logged$Proteins))

# 📥 Save the protein lists for DAVID
write.table(s1_first_proteins, "s1_significant_proteins.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(s2_first_proteins, "s2_significant_proteins.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(background_proteins, "background_proteins.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
