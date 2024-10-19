# This script takes normalized files and creates plots showing switching gene count for the top 500 AD genes
# Load necessary libraries
library(data.table)
library(ggplot2) # for plotting
library(readxl)
library(dplyr)  # for data manipulation

# Define input and output directories
input_dir <- "output"
plot_output_dir <- file.path("plots", "AD")
supplementary_file_path <- file.path("data", "supplementary_file.xlsx")
n <- 500  # Number of top genes to take

# Create plots output directory if it doesn't exist
if (!dir.exists(plot_output_dir)) {
  dir.create(plot_output_dir)
}

# Read the supplementary file and extract top 500 gene symbols by absolute sign value
supplementary_data <- read_excel(supplementary_file_path, sheet = "genes AD portraits")

# Assuming columns B and C correspond to `Gene.symbol...2` and `sign1...3`
top_genes <- supplementary_data %>%
  dplyr::select(GeneSymbol = `Gene.symbol...2`, Sign = `sign1...3`) %>%
  mutate(AbsSign = abs(Sign)) %>%
  arrange(desc(AbsSign)) %>%
  head(n) %>%
  pull(GeneSymbol)

# Get a list of all normalized files in the output directory
normalized_files <- list.files(input_dir, pattern = "normalized_.*\\.txt", full.names = TRUE)

# Function to generate and save the plot for a given file
generate_plot <- function(file_path) {
  # Extract tissue name from the file name
  tissue_name <- sub("normalized_(.*)\\.txt", "\\1", basename(file_path))
  
  # Set the threshold increment value
  t_increment <- 0.005
  
  # Read the input file
  data <- fread(file_path)
  
  # Filter the data to include only columns for the top 500 AD genes that exist in the current file
  existing_genes <- intersect(names(data), top_genes)
  if (length(existing_genes) == 0) {
    cat("No top genes found in:", tissue_name, "\n")
    return(NULL)
  }
  
  data <- data[, c("Age", existing_genes), with = FALSE]
  
  # Function to check if a gene column has "switched"
  has_switched <- function(column) {
    changes <- diff(column)
    switch_count <- sum(changes != 0)
    return(switch_count == 1)
  }
  
  # Initialize vectors to store threshold values and the number of switched genes
  threshold_values <- seq(t_increment, 1, by = t_increment)  # Start after 0
  switched_genes_counts <- numeric(length(threshold_values))
  
  # Iterate over the range of threshold values
  for (i in seq_along(threshold_values)) {
    t <- threshold_values[i]
    
    # Process the data by applying the threshold
    processed_data <- data
    processed_data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data), with = FALSE], function(x) ifelse(as.numeric(x) >= t, 1, 0))
    
    # Count the number of genes that have "switched"
    num_switched_genes <- sum(sapply(processed_data[, 2:ncol(data), with = FALSE], has_switched))
    
    # Store the result
    switched_genes_counts[i] <- num_switched_genes
  }
  
  # Create a data frame for plotting
  plot_data <- data.frame(Threshold = threshold_values, SwitchedGenes = switched_genes_counts)
  
  # Plot the results
  plot <- ggplot(plot_data, aes(x = Threshold, y = SwitchedGenes)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    theme(plot.title = element_blank(),  # Remove title
          axis.title.x = element_blank(),  # Remove x-axis label
          axis.title.y = element_blank(),  # Remove y-axis label
          text = element_text(size = 15 * 1.5))  # Increase font size by 50%
  
  # Save the plot
  plot_file_path <- file.path(plot_output_dir, paste0("plot_", tissue_name, ".png"))
  ggsave(plot_file_path, plot)
  
  cat("Generated plot for:", tissue_name, "\n")
}

# Generate plots for each normalized file
for (file_path in normalized_files) {
  generate_plot(file_path)
}

cat("All plots generated.\n")
