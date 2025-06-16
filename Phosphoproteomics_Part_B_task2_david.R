# Declaration:
# I, Samuel Rice, declare that I have used Chat GPT-3.5
# to enhance error handling within this R script.
# i used it to help deal with error messages when processing the david output
# Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)

# Define function to process DAVID output
process_david <- function(file_path) {
  # Read raw DAVID output
  david_raw <- read_excel(file_path, col_names = FALSE) %>% 
    mutate_all(as.character)
  
  # Identify rows with Enrichment Scores and Empty Rows
  # helps split up the seperate tables and identify where they start and end
  # create a new column in the tables for enrichment scores
  david_raw <- david_raw %>%
    # Mutate adds new columns to the dataframe
    mutate(
      # the first column is checked for the enrichment score
      EnrichmentScore = ifelse(grepl("Enrichment Score", `...1`), 
                               # string extract takes the numeric score assocated with enrichment score
                               # \\d+ contains the digits
                               # \\. is the decimal place
                               # as numeric is used extract numeric type for stats
                               as.numeric(str_extract(`...1`, "\\d+\\.\\d+")), NA),
      # Error handling step to identify empty rows, preventing error
      # this is as there are empty rows between the clusters
      IsEmpty = ifelse(is.na(`...1`) | `...1` == "", TRUE, FALSE)
    )
  
  # Cluster_ID assigned based on which cluster it is part of
  # this helps retain the origincal format of seperate clusters
  david_raw <- david_raw %>%
    # cumsum computes the cumulative sum everytime a new cluster is detected
    # this creates the cluster id so that each cluster is labelled with a unique identifier
    mutate(Cluster_ID = cumsum(!is.na(EnrichmentScore) | IsEmpty)) %>%
    # The enrichment score is filled into
    # the most recent NA (cluster tables)
    fill(EnrichmentScore)  
  
  # Remove unnecessary enrichment score lines, empty rows or headers
  david_cleaned <- david_raw %>%
    # grepl checks for the presence of headers or empty lines that need to be removed 
    filter(!grepl("Annotation Enrichment Score|Category|Term|Count|%", `...1`) & !IsEmpty)
  
  # change column names
  colnames(david_cleaned) <- c("Category", "Term", "Count", "Percent", "PValue", "Genes", 
                               "List_Total", "Pop_Hits", "Pop_Total", "Fold_Enrichment", 
                               "Bonferroni", "Benjamini", "FDR", "EnrichmentScore", "IsEmpty", "Cluster_ID")
  
  # Convert numeric columns for stats processing
  numeric_cols <- c("Count", "Percent", "PValue", "Fold_Enrichment", "Bonferroni", "Benjamini", "FDR", "EnrichmentScore", "Cluster_ID")
  # mutate all the specified columns to convert them
  david_cleaned <- david_cleaned %>% mutate(across(all_of(numeric_cols), as.numeric))
  
  # Extract and assign enrichment scores from 'Term' column
  # this is as they saved incorrectly 
  # it uses the same logic as earlier
  david_cleaned <- david_cleaned %>%
    mutate(Extracted_Score = ifelse(str_starts(Term, "Enrichment Score"), 
                                    as.numeric(str_extract(Term, "\\d+\\.\\d+")), NA)) %>%
    # the data is grouped by cluster id to ensure seperate processing
    group_by(Cluster_ID) %>%
    # fill is use to fill the rest of the cluster with the cluster enrichment score
    fill(Extracted_Score, .direction = "downup") %>%
    # this ungroups the data 
    ungroup() %>%
    # if the extracted score is abvailable, replace the enrichment score column
    # if not, keep the original enrichment core
    mutate(EnrichmentScore = ifelse(!is.na(Extracted_Score), Extracted_Score, EnrichmentScore)) %>%
    select(-Extracted_Score, -IsEmpty)  # Remove unnecessary columns
  
  # Remove rows where Term starts with Enrichment Score to further remove unnecessary rows
  david_final <- david_cleaned %>%
    filter(!str_starts(Term, "Enrichment Score"))
  
  return(david_final) # return the cleaned version of the david output
}

# run the function on both david outputs
david_final1 <- process_david("s1_david_clustering_output.xlsx")
david_final2 <- process_david("s2_david_clustering_output.xlsx")

# Save cleaned datasets to csvs
write.csv(david_final1, "s1_david_final_cleaned.csv", row.names = FALSE)
write.csv(david_final2, "s2_david_final_cleaned.csv", row.names = FALSE)

# Function to understand david results
# function will summarise and visualise the results
analyse_david_results <- function(david_data, title_prefix) {
  # Count unique clusters
  n_clusters <- length(unique(david_data$Cluster_ID))
  print(paste(title_prefix, "- Total Clusters Identified:", n_clusters))
  
  # sort the top 10 most significant terms by pvalue
  top_significant_terms <- david_data %>%
    arrange(PValue) %>%
    select(Cluster_ID, Term, PValue, EnrichmentScore, Fold_Enrichment) %>%
    head(10)
  print(top_significant_terms)
  
  # Summarise fold enrichment scores
  print(summary(david_data$Fold_Enrichment))
  
  # Bar plot of enrichment scores by cluster ID
  # visualises the enrichment of overall clusters
  ggplot(david_data %>% distinct(Cluster_ID, EnrichmentScore), aes(x = factor(Cluster_ID), y = EnrichmentScore)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    labs(title = paste(title_prefix, "- Enrichment Score by Cluster"), x = "Cluster ID", y = "Enrichment Score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Run the analyse david results function for both datasets
analyse_david_results(david_final1, "SARS-CoV-1")
analyse_david_results(david_final2, "SARS-CoV-2")

