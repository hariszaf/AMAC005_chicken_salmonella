install.packages("BiocManager")
BiocManager::install("RCy3")

# library(tidyr)
library(RCX)
library(RCy3)
library(igraph)
library(dplyr)
library(tibble)
library(tidyverse)

setwd("github_repos/contribute-to/AMAC005_chicken_salmonella/")



library(dplyr)
library(scales)
library(igraph)


# Build adj table
adjacency_matrix_frox_cx2 <- function(cxFile) {

  # Read the CX file
  rcx <- readCX(cxFile)

  # Convert to igraph object
  ig <- RCX::toIgraph(rcx)

  # Get edge list and nodes
  df <- igraph::get.data.frame(ig)
  edgelist <- df[, c("from", "to", "microbetag::weight")]
  nodes <- igraph::as_data_frame(ig, what = "vertices")

  # Merge 'from' node names
  edgelist <- merge(edgelist, nodes[, c("id", "nodeName")], by.x = "from", by.y = "id", all.x = TRUE)
  colnames(edgelist)[ncol(edgelist)] <- "from_nodeName"

  # Merge 'to' node names
  edgelist <- merge(edgelist, nodes[, c("id", "nodeName")], by.x = "to", by.y = "id", all.x = TRUE)
  colnames(edgelist)[ncol(edgelist)] <- "to_nodeName"

  # Filter out NA values in weight column
  edgelist <- edgelist %>%
    filter(!is.na(`microbetag::weight`))

  # Create 3-column edge list
  edgelist_3col <- edgelist[, c("from_nodeName", "to_nodeName", "microbetag::weight")]
  colnames(edgelist_3col) <- c("from", "to", "weight")

  # Ensure symmetry by duplicating edges
  symmetric_edges <- edgelist_3col %>%
    bind_rows(edgelist_3col %>% rename(from = to, to = from))

  # Get unique nodes
  nodes <- unique(c(symmetric_edges$from, symmetric_edges$to))

  # Ensure weight column is numeric
  symmetric_edges <- symmetric_edges %>%
    mutate(weight = as.numeric(weight))

  # Create a data frame with all possible combinations of nodes
  grid <- expand_grid(from = nodes, to = nodes)

  # Join with existing edges and replace missing values with 0
  adj_matrix <- grid %>%
    left_join(symmetric_edges, by = c("from", "to")) %>%
    mutate(weight = replace_na(weight, 0)) %>%
    pivot_wider(names_from = to, values_from = weight) %>%
    column_to_rownames(var = "from") %>%
    as.matrix()

  # Set diagonal to 1
  diag(adj_matrix) <- 1

  # Return the adjacency matrix
  return(adj_matrix)
}




# Set size, colors etc
set_node_attributes <- function(graph_pos, genome_metadata, genome_counts_filt, order_colors, default_size = 8) {
  # Ensure genome_counts_filt and genome_metadata are compatible with graph_pos names
  nodes <- V(graph_pos)$name

  genome_metadata <- add_missing_records(nodes, genome_metadata)

  # Node color
  V(graph_pos)$color <- tibble(
    order = genome_metadata$order %>% unique() %>% sort(),
    color = order_colors
  ) %>%
    right_join(genome_metadata, by = "order") %>%
    filter(genome %in% nodes) %>%
    select(genome, color) %>%
    pull(color)

  # Nodes not present in genome_counts_filt
  missing_nodes <- setdiff(nodes, genome_counts_filt$genome)

  # Calculate sizes only for nodes that are in genome_counts_filt
  sizes <- genome_counts_filt %>%
    mutate_at(vars(-genome), ~ . / sum(.)) %>%
    rowwise() %>%
    mutate(average = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    ungroup() %>%
    select(genome, average) %>%
    filter(genome %in% nodes) %>%
    arrange(match(genome, nodes)) %>%
    mutate(size = scales::rescale(average, to = c(2, 10))) %>%
    pull(size)

  # Set a specific size for the missing nodes
  size_vector <- rep(default_size, length(nodes))  # Start with default sizes for all nodes
  size_vector[nodes %in% genome_counts_filt$genome] <- sizes  # Assign calculated sizes to present nodes

  V(graph_pos)$size <- size_vector

  # Set shape for missing nodes (square for missing, circle for others)
  V(graph_pos)$shape <- ifelse(nodes %in% missing_nodes, "square", "circle")  # Or "hexagon"

  return(graph_pos)
}



# Add metavars to genome_metadata
add_missing_records <- function(nodes, genome_metadata) {
  # Get the difference (nodes not in genome_metadata$genome)
  missing_genomes <- setdiff(nodes, genome_metadata$genome)

  # Create a new row for each missing genome
  new_rows <- tibble(
    genome = missing_genomes,
    domain = "metavar",
    phylum = "metavar",
    class = "metavar",
    order = "metavar",
    family = "metavar",
    genus = "metavar",
    species = "metavar",
    completeness = 100,
    contamination = 0.2,
    length = 0
  )

  # Append the new rows to the existing genome_metadata
  genome_metadata <- bind_rows(genome_metadata, new_rows)

  return(genome_metadata)
}







# ----------------------------------

#
# ig2 <- createIgraphFromNetwork("microbetag_analysis/microbetag_nets/microbetag_net_day_14.cx2")
# ig2 <- importNetworkFromFile("microbetag_analysis/microbetag_nets/microbetag_net_day_14.cx2")
#
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("RCyjs")
#
#
#
# #sif <- "/home/luna.kuleuven.be/u0156635/github_repos/contribute-to/AMAC005_chicken_salmonella/microbetag_analysis/microbetag_nets/day7.sif"
#
# sif <- "microbetag_analysis/microbetag_nets/day7.sif"
#
# oncytoscape <- importNetworkFromFile(sif)
#
# #  You need to keep cytoscape in the network of interest
# ig <- createIgraphFromNetwork(oncytoscape)
#
#
# edgelist <- igraph::get.edgelist(ig)
#
# edge_attributes <- igraph::edge_attr(ig)
