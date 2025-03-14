install.packages("BiocManager")
BiocManager::install("RCy3")
# library(RCyjs)
library(RCX)
library(RCy3)
library(dplyr)
library(igraph)


setwd("github_repos/contribute-to/AMAC005_chicken_salmonella/")
remotes::install_gitlab("artemklevtsov/uchardet@devel")

# asda
cxFile <- "microbetag_analysis/microbetag_nets/microbetag_net_day_7.cx"
rcx = readCX(cxFile)

ig = RCX::toIgraph(rcx)

#
df <- igraph::get.data.frame(ig)
edgelist <- df[, c("from", "to", "microbetag::weight")]
nodes <- as_data_frame(ig, what = "vertices")

#
edgelist <- merge(edgelist, nodes[, c("id", "nodeName")], by.x = "from", by.y = "id", all.x = TRUE)
colnames(edgelist)[ncol(edgelist)] <- "from_nodeName"  # Rename to avoid confusion

# Merge 'to' column with 'nodes' to get 'to' node names
edgelist <- merge(edgelist, nodes[, c("id", "nodeName")], by.x = "to", by.y = "id", all.x = TRUE)
colnames(edgelist)[ncol(edgelist)] <- "to_nodeName"  # Rename

#
edgelist <- edgelist %>%
  filter(!is.na(`microbetag::weight`))

edgelist_3col <- edgelist[, c("from_nodeName", "to_nodeName", "microbetag::weight")]

#
adj_matrix <- edgelist_3col %>%
    pivot_wider(names_from = to_nodeName, values_from = `microbetag::weight`) %>%
    column_to_rownames(var = "from_nodeName")

# Replace NA with 0 using base R
adj_matrix[is.na(adj_matrix)] <- 0






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
