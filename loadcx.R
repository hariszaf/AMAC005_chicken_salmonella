# install.packages("BiocManager")
# BiocManager::install("RCy3")
# setwd("github_repos/contribute-to/AMAC005_chicken_salmonella/")


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


# Graph metrics
graph_metrics <- function(graph, col_name, name){
  return(
    tibble(
      !!col_name := name,  # Unquote the symbol to use it as a column name
      density=graph.density(graph),
      modularity=modularity(cluster_walktrap(graph)),
      assortability=assortativity_degree(graph),
      connectivity=mean(components(graph)$csize),
      centrality=eigen_centrality(graph)$value
    )
  )
}


# Process a file
process_cxfile <- function(cxFile, col_name, net_name, genome_metadata, genome_counts_filt, order_colors, threshold = 0.7) {

  # --------------
  # POSITIVE EDGES
  # --------------

  # Load the adjacency matrix (replace with the correct function to read CX files)
  mag_cor <- adjacency_matrix_frox_cx2(cxFile)  # Replace with actual function

  # Apply threshold to correlations (positive correlations only)
  # NOTE: for the case of overall:

  # if ( endsWith(cxFile, "overall.cx") | endsWith(cxFile, "microbetag_net_TG2.cx")  ) {
  w <- matrix(unlist(mag_cor), nrow = nrow(mag_cor), ncol = ncol(mag_cor))
  mag_cor_pos<- ifelse(w > threshold, 1, 0)
  dimnames(mag_cor_pos) <- dimnames(mag_cor)
  # }
  # mag_cor_pos <- ifelse(mag_cor > threshold, 1, 0)

  # Create an undirected graph from the adjacency matrix
  graph_pos <- igraph::graph_from_adjacency_matrix(mag_cor_pos, mode = "undirected", diag = FALSE)

  # Set node attributes (you should ensure these variables are defined)
  graph_pos <- set_node_attributes(graph_pos, genome_metadata, genome_counts_filt, order_colors)

  # Perform clustering using edge betweenness
  cluster_pos <- cluster_edge_betweenness(graph_pos)

  # Split communities based on clustering and filter for those with more than one member
  communities_pos <- split(V(graph_pos)$name, membership(cluster_pos)) %>% keep(~ length(.x) > 1)

  # --------------
  # NEGATIVE EDGES
  # --------------

  #Negative correlations
  mag_cor_neg <- ifelse(w < -threshold, 1, 0)
  dimnames(mag_cor_neg) <- dimnames(mag_cor)

  graph_neg <- graph_from_adjacency_matrix(mag_cor_neg, mode = "undirected", diag = FALSE)

  graph_neg <- set_node_attributes(graph_neg, genome_metadata, genome_counts_filt, order_colors)

  cluster_neg <- cluster_edge_betweenness(graph_neg)

  communities_neg <- split(V(graph_neg)$name, membership(cluster_neg) )%>% keep(~ length(.x) > 1)

  # --------------
  # METRICS
  # --------------
  col_name <- ensym(col_name)
  metrics_pos <- graph_metrics(graph_pos, col_name, net_name)
  metrics_neg <- graph_metrics(graph_neg, col_name, net_name)

  # Return both graph and communities
  return(
    list(
      graph_pos    = graph_pos,
      clusters_pos = cluster_pos,
      graph_neg    = graph_neg,
      clusters_neg = cluster_neg,
      metrics_pos  = metrics_pos,
      metrics_neg  = metrics_neg,
      communities_pos = communities_pos,
      communities_neg = communities_neg
    )
  )
}


# Plot a graph function along with its clusters
plot_graph <- function(graph, cluster, outputfile, min_community_size = 2, vertex_color = NULL, vertex_size = NULL, edge_width = 1, mark_col = "#f4f4f4") {

  png(outputfile, width = 2600, height = 2400, res  = 300)

  # Create the igraph plot
  p <- igraph::plot.igraph(
    graph,
    vertex.color = vertex_color %||% V(graph)$color,
    vertex.size = vertex_size %||% V(graph)$size,
    vertex.label = NA,
    vertex.frame.color = NA,
    edge.width = edge_width,
    mark.groups = communities(cluster) %>% keep(~ length(.x) >= min_community_size),
    mark.col = mark_col,
    mark.border = NA,
    layout = layout_with_fr
  )

  dev.off()
}



#
compute_cluster_distances <- function(clusters, genome_gifts) {
  table(map_chr(clusters, ~ paste(.x, collapse = " - "))) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rename(genomes = 1, count = 2) %>%
    mutate(distance = map_dbl(genomes, function(pair) {
      genomes <- str_split(pair, " - ", simplify = TRUE) %>% as.vector()
      genomes <- genomes[genomes %in% rownames(genome_gifts)]  # Keep only valid genomes

      if (length(genomes) < 2) return(NA_real_)  # Avoid errors if less than 2 genomes

      dist_value <- stats::dist(genome_gifts[genomes, , drop = FALSE], method = "manhattan") /
        ncol(genome_gifts[genomes, , drop = FALSE])
      return(as.numeric(dist_value))
    })) %>%
    arrange(-count)
}

