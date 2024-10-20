# This script takes normalized files and creates two bargraphs for each, along with updating pearson.txt
# Load necessary libraries
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)  # Load dplyr for data manipulation and pipe operator

# Define file paths and threshold value
input_dir <- "output"
output_dir <- "graphs"
supplementary_file_path <- file.path("data", "supplementary_file.xlsx")
t <- 0.5
n <- 500 # number of genes from supplementary file to take

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

# Function to process each file
process_file <- function(file_path) {
  # Read the input file
  data <- fread(file_path)
  
  # Apply the threshold to create a binary matrix
  binary_data <- data
  binary_data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data), with = FALSE], function(x) ifelse(as.numeric(x) >= t, 1, 0))
  
  # Determine the switching age group for each gene
  switching_age_groups <- sapply(binary_data[, 2:ncol(binary_data), with = FALSE], get_switching_age_group)
  
  # Convert the switching age group indices to age group labels
  age_groups <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99")
  switching_age_labels <- age_groups[switching_age_groups]
  
  # Create a data frame for plotting the complete set
  plot_data <- data.frame(AgeGroup = factor(switching_age_labels, levels = age_groups))
  
  # Remove NA values
  plot_data <- na.omit(plot_data)
  
  # Create the bar graph for the complete set
  plot <- ggplot(plot_data, aes(x = AgeGroup)) +
    geom_bar() +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      axis.text.x = element_text(size = 16),  # Increase font size for x-axis
      axis.text.y = element_text(size = 16)   # Increase font size for y-axis
    )
  
  
  # Print the plot
  print(plot)
  
  # Save the plot to a file with the updated naming convention
  tissue_name <- sub("normalized_", "", sub(".txt$", "", basename(file_path)))
  plot_file_name <- paste0("all_", tissue_name, "_bargraph_t", t, ".png")
  plot_file_path <- file.path(output_dir, plot_file_name)
  ggsave(plot_file_path, plot)
  
  cat("Bar graph for all genes generated and saved as:", plot_file_path, "\n")
  
  # Filter the binary data to include only the top n genes
  top_gene_columns <- which(names(binary_data) %in% top_genes)
  binary_data_top <- binary_data[, c(1, top_gene_columns), with = FALSE]
  
  # Determine the switching age group for each of the top n genes
  switching_age_groups_top <- sapply(binary_data_top[, 2:ncol(binary_data_top), with = FALSE], get_switching_age_group)
  
  # Convert the switching age group indices to age group labels
  switching_age_labels_top <- age_groups[switching_age_groups_top]
  
  # Create a data frame for plotting the top n genes
  plot_data_top <- data.frame(AgeGroup = factor(switching_age_labels_top, levels = age_groups))
  
  # Remove NA values
  plot_data_top <- na.omit(plot_data_top)
  
  # Create the bar graph for the top n genes
  plot_top <- ggplot(plot_data_top, aes(x = AgeGroup)) +
    geom_bar() +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      axis.text.x = element_text(size = 16),  # Increase font size for x-axis
      axis.text.y = element_text(size = 16)   # Increase font size for y-axis
    )
  
  # Print the plot
  print(plot_top)
  
  # Save the plot for top n genes to a file with the updated naming convention
  plot_top_file_name <- paste0("ad_", tissue_name, "_bargraph_t", t, ".png")
  plot_top_file_path <- file.path(output_dir, plot_top_file_name)
  ggsave(plot_top_file_path, plot_top)
  
  cat("Bar graph for top", n, "genes generated and saved as:", plot_top_file_path, "\n")
  
  # Calculate Pearson correlation coefficient
  counts_all <- table(plot_data$AgeGroup)
  counts_top <- table(plot_data_top$AgeGroup)
  
  # Ensure both tables have the same levels
  all_levels <- union(names(counts_all), names(counts_top))
  counts_all <- as.numeric(counts_all[all_levels])
  counts_top <- as.numeric(counts_top[all_levels])
  
  # Check if standard deviations are zero
  if (sd(counts_all) != 0 && sd(counts_top) != 0) {
    # Compute Pearson correlation coefficient
    pearson_r <- cor(counts_all, counts_top, method = "pearson")
  } else {
    pearson_r <- NA
  }
  
  # Append the Pearson correlation value to pearson.txt
  pearson_file_path <- file.path(output_dir, "pearson.txt")
  entry <- sprintf("%s, %.3f, %.3f\n", tissue_name, t, pearson_r)
  
  # Read the existing entries, if any
  if (file.exists(pearson_file_path)) {
    existing_entries <- readLines(pearson_file_path)
  } else {
    existing_entries <- character(0)
  }
  
  # Add the new entry and sort alphabetically
  all_entries <- c(existing_entries, entry)
  sorted_entries <- sort(all_entries)
  
  # Write the sorted entries back to the file
  writeLines(sorted_entries, pearson_file_path)
  
  cat("Pearson correlation value appended to:", pearson_file_path, "\n")
}

# Process all normalized files in the output directory
normalized_files <- list.files(input_dir, pattern = "normalized_.*\\.txt", full.names = TRUE)
for (file_path in normalized_files) {
  process_file(file_path)
}
