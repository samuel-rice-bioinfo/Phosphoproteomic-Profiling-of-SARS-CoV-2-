# Declaration:
# I, Samuel Rice, declare that I have used Chat GPT-3.5
# to enhance error handling within this R script.
# Install packages
install.packages(c("readxl", "ggplot2", "writexl", "htmlTable"))

# Load libraries
library(readxl)
library(ggplot2)
library(writexl)
library(htmlTable)
library(dplyr)
library(tidyverse)
# Load data and check for missing values
df <- read.csv("Animal_20200325_TM_HStmtpro_CoV12_ph2_fr3.csv", skip = 1, header=TRUE)
print(sum(is.na(df$protein))) # missing proteins
print(sum(is.na(df$e.value))) # missing evals

# check column names
colnames(df)

# Order by e value in ascending order
df <- df[order(df$e.value), ]


# Initialize empty vectors to store tp/fp/fdr/phosphopeptide boolean value (TRUE/FALSE)
fp_v = vector(mode = "integer", length = nrow(df))  # False positives
tp_v = vector(mode = "integer", length = nrow(df))  # True positives
fdr_v = vector(mode = "numeric", length = nrow(df))  # False Discovery Rate
isphospho_v = vector(mode = "logical", length = nrow(df))  # Phosphopeptide identifier

# Initialize FDR threshold vectors
# for hits
FDR_01_matches = numeric(nrow(df))
FDR_05_matches = numeric(nrow(df))
FDR_1_matches = numeric(nrow(df))
# for phsophopeptides
FDR_01_phos = numeric(nrow(df))
FDR_05_phos = numeric(nrow(df))
FDR_1_phos = numeric(nrow(df))

# Initialize fp and total to zero
fp = 0
total = 0

# For each row, check if its a decoy
# loop to calculate FDR and identify phosphopeptides
for (i in 1:nrow(df)) {
  prot_acc = df[i, "protein"]
  
  # Identify false positives 
  # if the line starts with decoy then add one to FP
  # else add one to the total
  if (substr(prot_acc, 1, 5) == "DECOY") {
    fp = fp + 1
  } else {
    total = total + 1
  }
  
  # Identify phosphopeptides with 79.96 modification
  pep_mod = df[i, "modifications"]
  isphospho_v[i] = grepl("79.96", pep_mod, fixed = TRUE)
  
  # Count tp and fp counts
  fp_v[i] = fp
  tp_v[i] = total - fp
  
  # Compute FDR using given equation
  # uses if statement to handle errors (0 handling)
  if ((tp_v[i] + fp_v[i]) == 0) {
    fdr_v[i] = 0
  } else {
    fdr_v[i] = fp_v[i] / (tp_v[i] + fp_v[i])
  }
  
  # Debugging step
  # print every 1000th row to track progress
  if (i %% 1000 == 0) {
    print(paste("Row:", i, "| TP:", tp_v[i], "| FP:", fp_v[i], "| FDR:", round(fdr_v[i], 4)))
  }
  
  # Identify hits within FDR thresholds
  # for each threshold, 
  # if FDR<threshold, add one to fdr matches
  # if it is a phosphopeptide, add one to FDR phos
  if (fdr_v[i] <= 0.05) {
    FDR_05_matches[i] = 1
    if (isphospho_v[i]) {
      FDR_05_phos[i] = 1
    }
  }
  if (fdr_v[i] <= 0.1) {
    FDR_1_matches[i] = 1
    if (isphospho_v[i]) {
      FDR_1_phos[i] = 1
    }
  }
  if (fdr_v[i] <= 0.01) {
    FDR_01_matches[i] = 1
    if (isphospho_v[i]) {
      FDR_01_phos[i] = 1
    }
  }
}

# Assign counted values to dataframe
df$count_fp = fp_v
df$count_tp = tp_v
df$is_phospho = isphospho_v
df$fdr = fdr_v
# Add matches and phosphopeptides to
# the dataframe in new columns
df$FDR_01_matches = FDR_01_matches
df$FDR_05_matches = FDR_05_matches
df$FDR_1_matches = FDR_1_matches
df$FDR_01_phos = FDR_01_phos
df$FDR_05_phos = FDR_05_phos
df$FDR_1_phos = FDR_1_phos

# Calculate total counts for FDR thresholds
# for summary table
df$FDR_01_matches[1] = sum(fdr_v <= 0.01, na.rm = TRUE)
df$FDR_05_matches[1] = sum(fdr_v <= 0.05, na.rm = TRUE)
df$FDR_1_matches[1] = sum(fdr_v <= 0.1, na.rm = TRUE)
df$FDR_01_phos[1] = sum(fdr_v <= 0.01 & isphospho_v, na.rm = TRUE)
df$FDR_05_phos[1] = sum(fdr_v <= 0.05 & isphospho_v, na.rm = TRUE)
df$FDR_1_phos[1] = sum(fdr_v <= 0.1 & isphospho_v, na.rm = TRUE)

# Create summary table that can show the counts
summary_table <- data.frame(
  FDR_Threshold = c("1%", "5%", "10%"),
  PSM_Count = c(df$FDR_01_matches[1], df$FDR_05_matches[1], df$FDR_1_matches[1]),
  Phosphopeptide_Count = c(df$FDR_01_phos[1], df$FDR_05_phos[1], df$FDR_1_phos[1])
)
# print the table
print("✅ Summary Table:")
print(summary_table)

# Save summary table as xlsx
write_xlsx(summary_table, "final_FDR_summary.xlsx")

# Generate HTML table for saving to the report
html_table <- htmlTable(summary_table, caption = "Final FDR Summary Table")
write(html_table, file = "final_FDR_summary.html")

# Save df to new file
write.csv(df, "final_Phosphopeptide_results.csv", row.names = FALSE)

# Debugging step
# ensuring FDR values are within expected range
# range of 0-1
if (max(df$fdr, na.rm = TRUE) > 1 || min(df$fdr, na.rm = TRUE) < 0) {
  stop("FDR values are outside expected range (0 to 1)")
} else {
  print("FDR values within expected range")
}

# Create scatter plot for fdr vs true positives
ggplot(df, aes(x = count_tp, y = fdr)) +
  geom_point(color = "blue") +
  labs(title = "FDR vs True Positives",
       x = "True Positive Count",
       y = "FDR") +
  theme_minimal()

# Save plot as PNG
ggsave("final_FDR_plot.png")


#####################
# statistics for results section
# Calculate additional statistics
# count number of decoys, tps, fp, psms, and phosphopeptides
total_decoys <- sum(grepl("^DECOY", df$protein))  # Count of decoy hits
total_true_positives <- max(df$count_tp)  # Total true positives
total_false_positives <- max(df$count_fp)  # Total false positives
total_psms <- nrow(df)  # Total PSMs
total_phosphopeptides <- sum(df$is_phospho)  # Count of phosphopeptides



# Percent stat (%/100) to 2 decimal places
percent_phosphopeptides <- round((total_phosphopeptides / total_psms) * 100, 2)
percent_decoys <- round((total_decoys / total_psms) * 100, 2)
percent_true_positives <- round((total_true_positives / total_psms) * 100, 2)
percent_false_positives <- round((total_false_positives / total_psms) * 100, 2)

# Print stats
print(paste("Total PSMs:", total_psms))
print(paste("Total Decoy Hits:", total_decoys, "| Percentage of Decoys:", percent_decoys, "%"))
print(paste("Total True Positives:", total_true_positives, "| Percentage of True Positives:", percent_true_positives, "%"))
print(paste("Total False Positives:", total_false_positives, "| Percentage of False Positives:", percent_false_positives, "%"))
print(paste("Total Phosphopeptides Identified:", total_phosphopeptides, "| Percentage of Phosphopeptides:", percent_phosphopeptides, "%"))

# Save stats to dataframe
final_stats <- data.frame(
  Metric = c("Total PSMs", "Total Decoy Hits", "Total True Positives", "Total False Positives", "Total Phosphopeptides"),
  Count = c(total_psms, total_decoys, total_true_positives, total_false_positives, total_phosphopeptides),
  Percentage = c(100, percent_decoys, percent_true_positives, percent_false_positives, percent_phosphopeptides)
)

# Save stats as excel and HTML for varied output options
write_xlsx(final_stats, "final_Phosphopeptide_Statistics.xlsx")
write(htmlTable(final_stats, caption = "Final Phosphopeptide Statistics"), file = "final_Phosphopeptide_Statistics.html")
 