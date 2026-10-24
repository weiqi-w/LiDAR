CreateNetworks <- function(lidardata, n) {
  nets = list()

  if (length(n) > 1) {
    for (i in seq_along(n)) {
      cat()
      nets[[i]] = NetworkFromMCA(lidardata, nn = n[[i]])
    }
  } else {
    nets[[1]] = NetworkFromMCA(lidardata, nn = n) #20260603
  }

  lidardata@Coembedding@nets = nets
  lidardata@Coembedding@n_list = n
  return(lidardata)
}

##Construct network----
NetworkFromMCA = function(lidardata, nn = 30){
  total <- nrow(lidardata@Coembedding@cell_embedding) + nrow(lidardata@Coembedding@gene_embedding)
  if (length(nn) != 1L || !is.finite(nn) || nn < 1 || nn != as.integer(nn) || nn >= total) stop("Network neighborhood size is invalid.")
  coembed = rbind(lidardata@Coembedding@cell_embedding,
                  lidardata@Coembedding@gene_embedding)

  # Compute the nearest neighbors for the coembed matrix, including self as the first neighbor
  nbrs = nn2(coembed, k = nn + 1)
  nbrs_idx = nbrs$nn.idx

  # Create an edge list where each row represents an edge between two cells (neighbors)
  edge_list = cbind(nbrs_idx[, 1], c(nbrs_idx[, -1]))

  # Create an edge list with cell namesr
  edge_list_name = cbind(rownames(coembed)[edge_list[, 1]],
                         rownames(coembed)[edge_list[, 2]])
  net = graph_from_edgelist(as.matrix(edge_list_name), directed = F)
  return(net)
}



#Run RWR----
