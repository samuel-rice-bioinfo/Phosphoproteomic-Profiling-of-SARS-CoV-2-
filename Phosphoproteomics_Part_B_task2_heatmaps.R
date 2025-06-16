# Declaration:
# I, Samuel Rice, declare that I have used Chat GPT-3.5
# to enhance error handling within this R script.
# Load libraries
library(dplyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(knitr)

# Define a function to generate visualisations from cleaned david output
# produces heatmaps and barplots
generate_heatmap <- function(file_path, title_prefix = "SARS-CoV") {
  # Load cleaned DAVID file
  david_cleaned <- read.csv(file_path)
  
  # Filter significant terms by pvalue
  significant_terms <- david_cleaned %>%
    filter(PValue < 0.05) %>%
    arrange(PValue)
  
  # Select top 10 significant terms
  top_enriched_terms <- significant_terms %>%
    head(10)
  

  # Create a simplified and clearly labeled summary table
  top_terms_summary <- top_enriched_terms %>%
    transmute(
      `Cluster ID` = Cluster_ID,
      `Biological Term` = Term,
      `Fold Enrichment (Term-specific)` = round(Fold_Enrichment, 2),
      `Cluster Enrichment Score (Overall Cluster)` = round(EnrichmentScore, 2),
      `E-value (Significance)` = round(PValue, 2)
    )
  
  # Save the table as an HTML file
  kable(top_terms_summary, format = "html", table.attr = "style='width:100%;'") %>%
    writeLines(paste0(title_prefix, "_top_terms_summary.html"))
  
  # Prepare data for heatmap by converting to wide format using spread
  heatmap_data <- significant_terms %>%
    select(Term, Cluster_ID, Fold_Enrichment) %>%
    spread(key = Cluster_ID, value = Fold_Enrichment, fill = 0)
  
  # Convert to matrix and set the term as row names
  row.names(heatmap_data) <- heatmap_data$Term
  heatmap_matrix <- as.matrix(heatmap_data[, -1])
  
  # Plot the heatmap to visualise fold enrichment scores of terms by cluster
  pheatmap(heatmap_matrix,
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           main = paste("Fold Enrichment Scores of Specific Biological Terms by Functional Cluster Associated with", title_prefix, "Infection"),
           color = colorRampPalette(c("blue", "white", "red"))(50),
           angle_col = 45)
  
  # Summarize clusters by mean pval for general stats
  top_clusters <- significant_terms %>%
    group_by(Cluster_ID) %>%
    summarise(Average_PValue = mean(PValue), Top_Term = Term[which.min(PValue)]) %>%
    arrange(Average_PValue)
  
  # Return the summarized top cluster data
  return(top_clusters)
}

# calling the function for S1 and S2 cleaned david data
sars_cov1_clusters <- generate_heatmap("s1_david_final_cleaned.csv", "SARS-CoV-1")
sars_cov2_clusters <- generate_heatmap("s2_david_final_cleaned.csv", "SARS-CoV-2")
