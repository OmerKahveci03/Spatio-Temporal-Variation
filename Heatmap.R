# This script performs a Jaccard similarity for the genes between every combination of tissues
library(dplyr)
library(readxl)
library(writexl)
library(ComplexHeatmap)
library(circlize)  # For colorRamp2 function

# Define directories and output file paths
sets_dirs <- list(all = "all_sets", ad = "ad_sets")

# Define a "heatmaps" directory and create it if it doesn't exist
heatmap_dir <- "heatmaps"
if (!dir.exists(heatmap_dir)) {
  dir.create(heatmap_dir)
}

output_files <- list(
  all_age = file.path("output", "heatmap_all_age.xlsx"),
  all_life = file.path("output", "heatmap_all_life.xlsx"),
  ad_age = file.path("output", "heatmap_ad_age.xlsx"),
  ad_life = file.path("output", "heatmap_ad_life.xlsx")
)

heatmap_files <- list(
  all_age = file.path(heatmap_dir, "heatmap_all_age.png"),
  all_life = file.path(heatmap_dir, "heatmap_all_life.png"),
  ad_age = file.path(heatmap_dir, "heatmap_ad_age.png"),
  ad_life = file.path(heatmap_dir, "heatmap_ad_life.png")
)

# Function to read sets files and compute Jaccard similarities
process_sets <- function(sets_dir, output_age_file, heatmap_age_file, output_life_file, heatmap_life_file) {
  # List all the (tissue name)_sets.txt files in the sets directory
  set_files <- list.files(sets_dir, pattern = ".*_sets\\.txt$", full.names = TRUE)
  
  # Initialize a list to store data
  data_list <- list()
  
  # Read the data from each sets file
  for (file in set_files) {
    tissue_name <- gsub("_sets\\.txt$", "", basename(file))
    lines <- readLines(file)
    age_groups <- strsplit(lines, " ")
    data_list[[tissue_name]] <- lapply(age_groups, function(group) group[group != ""])
  }
  
  # Function to compute Jaccard similarity
  jaccard_similarity <- function(set1, set2) {
    intersection <- length(intersect(set1, set2))
    union <- length(union(set1, set2))
    if (union == 0) {
      return(0)
    }
    return(intersection / union)
  }
  
  # Initialize matrices to store the results for age and life
  tissue_names <- names(data_list)
  num_tissues <- length(tissue_names)
  results_matrix_age <- matrix(0, nrow = num_tissues, ncol = num_tissues, dimnames = list(tissue_names, tissue_names))
  results_matrix_life <- matrix(0, nrow = num_tissues, ncol = num_tissues, dimnames = list(tissue_names, tissue_names))
  
  # Compute the Jaccard similarity for each tissue pair for age groups
  for (i in 1:num_tissues) {
    for (j in 1:num_tissues) {
      set_i <- data_list[[i]]
      set_j <- data_list[[j]]
      similarities <- sapply(1:5, function(k) jaccard_similarity(set_i[[k]], set_j[[k]]))
      average_similarity <- mean(similarities)
      results_matrix_age[i, j] <- average_similarity
    }
  }
  
  # Compute the Jaccard similarity for each tissue pair for entire life set
  for (i in 1:num_tissues) {
    for (j in 1:num_tissues) {
      set_i <- unlist(data_list[[i]])
      set_j <- unlist(data_list[[j]])
      similarity <- jaccard_similarity(set_i, set_j)
      results_matrix_life[i, j] <- similarity
    }
  }
  
  # Convert the results matrices to data frames
  results_df_age <- as.data.frame(results_matrix_age)
  results_df_age <- cbind(Tissue = rownames(results_df_age), results_df_age)
  
  results_df_life <- as.data.frame(results_matrix_life)
  results_df_life <- cbind(Tissue = rownames(results_df_life), results_df_life)
  
  # Write the results to Excel files
  write_xlsx(results_df_age, output_age_file)
  write_xlsx(results_df_life, output_life_file)
  
  # Perform hierarchical clustering for age and life matrices
  hc_age <- hclust(dist(results_matrix_age))
  hc_life <- hclust(dist(results_matrix_life))
  
  # Method 1: Using colorRamp2 from circlize package
  if (exists("colorRamp2")) {
    color_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  } else {
    # Method 2: Using colorRampPalette from base R as a fallback
    color_fun <- colorRampPalette(c("blue", "white", "red"))
  }
  
  # Create the heatmap for age with dendrograms
  heatmap_age <- Heatmap(results_matrix_age, 
                         name = "Jaccard Similarity (Age)", 
                         cluster_rows = hc_age, 
                         cluster_columns = hc_age, 
                         show_row_dend = TRUE, 
                         show_column_dend = TRUE,
                         col = color_fun,
                         column_names_rot = 45)
  
  # Save the age heatmap to a file
  png(heatmap_age_file, width = 1000, height = 800)
  draw(heatmap_age, heatmap_legend_side = "right")
  dev.off()
  
  # Create the heatmap for life with dendrograms
  heatmap_life <- Heatmap(results_matrix_life, 
                          name = "Jaccard Similarity (Life)", 
                          cluster_rows = hc_life, 
                          cluster_columns = hc_life, 
                          show_row_dend = TRUE, 
                          show_column_dend = TRUE,
                          col = color_fun,
                          column_names_rot = 45)
  
  # Save the life heatmap to a file
  png(heatmap_life_file, width = 1000, height = 800)
  draw(heatmap_life, heatmap_legend_side = "right")
  dev.off()
  
  cat("Heatmap of tissue comparisons by age saved to:", output_age_file, "\n")
  cat("Heatmap image by age saved to:", heatmap_age_file, "\n")
  cat("Heatmap of tissue comparisons by life saved to:", output_life_file, "\n")
  cat("Heatmap image by life saved to:", heatmap_life_file, "\n")
}

# Process the all_sets and ad_sets directories
process_sets(sets_dirs$all, output_files$all_age, heatmap_files$all_age, output_files$all_life, heatmap_files$all_life)
process_sets(sets_dirs$ad, output_files$ad_age, heatmap_files$ad_age, output_files$ad_life, heatmap_files$ad_life)
