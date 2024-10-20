# This script creates individual sets files for each tissue
# Load necessary libraries
library(data.table)
library(readxl)
library(dplyr)  # Load dplyr for data manipulation and pipe operator

# Define file paths and threshold value
input_dir <- "output"
ad_sets_dir <- "ad_sets"
all_sets_dir <- "all_sets"
supplementary_file_path <- file.path("data", "supplementary_file.xlsx")
t <- 0.5
n <- 500 # number of genes from supplementary file to take

# Create the output directories if they don't exist
if (!dir.exists(ad_sets_dir)) {
  dir.create(ad_sets_dir)
}
if (!dir.exists(all_sets_dir)) {
  dir.create(all_sets_dir)
}

# Read the supplementary file and extract top n gene symbols by absolute sign value
supplementary_data <- read_excel(supplementary_file_path, sheet = "genes AD portraits")

# Assuming columns B and C correspond to `Gene.symbol...2` and `sign1...3`
top_genes <- supplementary_data %>%
  dplyr::select(GeneSymbol = `Gene.symbol...2`, Sign = `sign1...3`) %>%
  mutate(AbsSign = abs(Sign)) %>%
  arrange(desc(AbsSign)) %>%
  head(n) %>%
  pull(GeneSymbol)

# Function to determine the switching age group for a gene column
get_switching_age_group <- function(column) {
  changes <- diff(column)
  change_index <- which(changes != 0)[1]
  if (!is.na(change_index)) {
    return(change_index + 1)
  } else {
    return(NA)
  }
}

# Initialize the 2D array to store values for each tissue and age group
age_groups <- c("30-39", "40-49", "50-59", "60-69", "70-79")

# Function to process each file
process_file <- function(file_path, output_dir, gene_filter = NULL) {
  # Read the input file
  data <- fread(file_path)
  
  # Apply the threshold to create a binary matrix
  binary_data <- data
  binary_data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data), with = FALSE], function(x) ifelse(as.numeric(x) >= t, 1, 0))
  
  # Filter the binary data if a gene filter is provided
  if (!is.null(gene_filter)) {
    gene_columns <- which(names(binary_data) %in% gene_filter)
    binary_data <- binary_data[, c(1, gene_columns), with = FALSE]
  }
  
  # Determine the switching age group for each gene
  switching_age_groups <- sapply(binary_data[, 2:ncol(binary_data), with = FALSE], get_switching_age_group)
  
  # Convert the switching age group indices to age group labels
  age_group_labels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99")[switching_age_groups]
  
  # Create a data frame with gene symbols and their corresponding age group
  gene_symbols <- colnames(binary_data)[-1]
  switching_genes_data <- data.frame(GeneSymbol = gene_symbols, AgeGroup = age_group_labels, stringsAsFactors = FALSE)
  
  # Filter to include only the relevant age groups
  switching_genes_data <- switching_genes_data[switching_genes_data$AgeGroup %in% age_groups, ]
  
  # Prepare the output for each age group
  output_lines <- sapply(age_groups, function(age_group) {
    genes <- switching_genes_data$GeneSymbol[switching_genes_data$AgeGroup == age_group]
    paste(genes, collapse = " ")
  })
  
  # Prepare the output file path
  tissue_name <- sub("normalized_", "", sub(".txt$", "", basename(file_path)))
  output_file_path <- file.path(output_dir, paste0(tissue_name, "_sets.txt"))
  
  # Write the output lines to the file
  writeLines(output_lines, output_file_path)
  
  # Print the counts for this tissue (for immediate feedback)
  cat(tissue_name, ": ", paste(output_lines, collapse = " | "), "\n")
}

# Process all normalized files in the output directory for top n genes
normalized_files <- list.files(input_dir, pattern = "normalized_.*\\.txt", full.names = TRUE)
for (file_path in normalized_files) {
  process_file(file_path, ad_sets_dir, top_genes)
}

# Process all normalized files in the output directory for all genes
for (file_path in normalized_files) {
  process_file(file_path, all_sets_dir)
}
