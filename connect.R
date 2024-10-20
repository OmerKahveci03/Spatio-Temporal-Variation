# Creates connectivity heatmaps
library(data.table)
library(dplyr)
library(igraph)
library(biomaRt)
library(ggplot2)
library(reshape2)
library(openxlsx)

# Define file paths
input_dir <- "ad_sets"
protein_links_file <- "data/9606.protein.links.v12.0.txt"
output_dir <- "connectivity"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Read protein links file
protein_links <- fread(protein_links_file)

# Remove the "9606." prefix to get the Ensembl protein IDs
protein_links$protein1 <- gsub("9606\\.", "", protein_links$protein1)
protein_links$protein2 <- gsub("9606\\.", "", protein_links$protein2)

# Filter the links to include only those with a combined score of 900 or more
protein_links_filtered <- protein_links %>% filter(combined_score >= 900)

# Extract unique Ensembl protein IDs
ensembl_ids <- unique(c(protein_links_filtered$protein1, protein_links_filtered$protein2))

# Function to connect to Ensembl with mirror site support
connect_to_ensembl <- function() {
  ensembl <- NULL
  try({
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }, silent = TRUE)
  
  if (is.null(ensembl)) {
    mirrors <- c("useast", "uswest", "asia")
    for (mirror in mirrors) {
      try({
        ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = mirror)
        if (!is.null(ensembl)) {
          break
        }
      }, silent = TRUE)
    }
  }
  
  if (is.null(ensembl)) {
    stop("All Ensembl services are currently unavailable.")
  }
  
  return(ensembl)
}

# Connect to Ensembl
ensembl <- connect_to_ensembl()

# Function to map Ensembl IDs to gene symbols
map_ensembl_to_gene <- function(ensembl_ids) {
  mapped <- getBM(attributes = c("ensembl_peptide_id", "external_gene_name"), 
                  filters = "ensembl_peptide_id", 
                  values = ensembl_ids, 
                  mart = ensembl)
  return(mapped)
}

# Map Ensembl IDs to gene symbols in chunks
chunk_size <- 10000
results <- lapply(seq(1, length(ensembl_ids), by = chunk_size), function(i) {
  map_ensembl_to_gene(ensembl_ids[i:min(i + chunk_size - 1, length(ensembl_ids))])
})
result_df <- do.call(rbind, results)

# Create a mapping from Ensembl protein IDs to gene symbols
ensembl_to_gene <- setNames(result_df$external_gene_name, result_df$ensembl_peptide_id)

# Map the Ensembl IDs in the protein links to gene symbols
protein_links_filtered$gene1 <- ensembl_to_gene[protein_links_filtered$protein1]
protein_links_filtered$gene2 <- ensembl_to_gene[protein_links_filtered$protein2]

# Remove rows where gene symbols could not be mapped
protein_links_filtered <- protein_links_filtered %>% 
  filter(!is.na(gene1) & !is.na(gene2))

# Create an igraph object using the mapped protein links data
g <- graph_from_data_frame(protein_links_filtered[, c("gene1", "gene2", "combined_score")], directed = FALSE)

# Print graph information: number of nodes (genes) and edges (interactions)
cat("Graph information:\n")
cat("Number of nodes (genes):", vcount(g), "\n")
cat("Number of edges (interactions):", ecount(g), "\n")

# Initialize a data frame to store results for all tissues
all_results <- data.frame(Tissue = character(), AgeGroup = character(), AverageInverseDistance = numeric(), stringsAsFactors = FALSE)

# Process each tissue file in the input directory
tissue_files <- list.files(input_dir, pattern = "_sets.txt$", full.names = TRUE)

for (tissue_file in tissue_files) {
  tissue_name <- gsub("_sets.txt$", "", basename(tissue_file))
  
  # Replace underscores with spaces in tissue names
  tissue_name <- gsub("_", " ", tissue_name)
  
  # Read the gene sets file
  gene_sets <- readLines(tissue_file)
  
  # Initialize union of switching genes
  switching_genes_union <- character()
  
  # Initialize results list
  age_groups <- c("30-39", "40-49", "50-59", "60-69", "70-79")
  results <- data.frame(AgeGroup = age_groups, AverageInverseDistance = numeric(length(age_groups)), stringsAsFactors = FALSE)
  
  # Calculate distances and average inverse distances for each age group
  for (i in seq_along(gene_sets)) {
    current_age_group_genes <- unlist(strsplit(gene_sets[i], " "))
    
    # Update the union of switching genes
    switching_genes_union <- unique(c(switching_genes_union, current_age_group_genes))
    
    # Initialize distance vector
    distances <- rep(Inf, vcount(g))
    names(distances) <- V(g)$name
    
    # Calculate distances from each gene in the graph to the genes in the union of switching genes
    for (gene in switching_genes_union) {
      if (gene %in% V(g)$name) {
        dist <- distances(g, v = gene)
        distances <- pmin(distances, dist, na.rm = TRUE)
      }
    }
    
    # Replace infinite distances with zero
    distances[is.infinite(distances)] <- 0
    
    # Calculate the inverse distances
    inverse_distances <- ifelse(distances == 0, 0, 1 / distances)
    
    # Remove the switching genes from the calculation
    inverse_distances <- inverse_distances[!names(inverse_distances) %in% switching_genes_union]
    
    # Calculate average value (inverse distance)
    average_inverse_distance <- mean(inverse_distances)
    
    # Store the result
    results$AverageInverseDistance[i] <- average_inverse_distance
  }
  
  # Add tissue name to results
  results$Tissue <- tissue_name
  
  # Append results to all_results
  all_results <- rbind(all_results, results)
}

# Transform data for heatmap
heatmap_data <- dcast(all_results, Tissue ~ AgeGroup, value.var = "AverageInverseDistance")

# Replace underscores with spaces in tissue names in heatmap data
heatmap_data$Tissue <- gsub("_", " ", heatmap_data$Tissue)

# Order tissues by the value in the first age group (30-39) in descending order
heatmap_data <- heatmap_data[order(-heatmap_data$`30-39`), ]

# Convert data to long format for ggplot2
heatmap_long <- melt(heatmap_data, id.vars = "Tissue")

# Set the levels of the Tissue factor according to the sorted order
heatmap_long$Tissue <- factor(heatmap_long$Tissue, levels = heatmap_data$Tissue)

# Adjust the plot size and layout
p <- ggplot(heatmap_long, aes(x = variable, y = Tissue, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue", guide = guide_colorbar(title = NULL)) +  # Remove color scale title
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_text(size = 15 * 1.1, face = "bold"),  # Increase tissue name size and make it bold
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    plot.title = element_blank()  # Remove plot title
  )

# Save the heatmap with adjusted aspect ratio
ggsave(filename = file.path(output_dir, "heatmap_connectivity.png"), plot = p, width = 6, height = 15)

# Save the heatmap
ggsave(filename = file.path(output_dir, "heatmap_connectivity.png"), plot = p, width = 12, height = 10)

# Save the data used in the heatmap to an Excel file
write.xlsx(heatmap_data, file = file.path(output_dir, "heatmap_connectivity.xlsx"))

print("Heatmap and Excel file generated and saved.")
