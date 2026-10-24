ImputeExp <- function(lidardata, nn = 30) {
  cell_embedding = lidardata@Coembedding@cell_embedding
  #names(cell_embedding) = colnames(lidardata)
  cell_cell_nbrs = FindCellNeighbors(cell_embedding, nn = nn)  # Compute nearest neighbors and weights
  ligand_exp_imp = ImputeByNeighbors(cell_cell_nbrs, lidardata@Data@ligand_exp)
  receptor_exp_imp = ImputeByNeighbors(cell_cell_nbrs, lidardata@Data@receptor_exp)
  lidardata@ImputationData@ligand_exp = ligand_exp_imp
  lidardata@ImputationData@receptor_exp = receptor_exp_imp
  return(lidardata)
}

# Function to compute nearest neighbors and their weights
FindCellNeighbors <- function(embedding, nn = 30) {
  if (length(nn) != 1L || !is.finite(nn) || nn < 1 || nn != as.integer(nn) || nn >= nrow(embedding)) stop("nn must be an integer between 1 and number of cells minus one.")
  # Compute nearest neighbors (including self as first neighbor)
  nbrs = nn2(embedding, k = nn + 1)

  # Remove self
  nbrs[["nn.idx"]] = nbrs[["nn.idx"]][, -1, drop = FALSE]

  # Compute weights for the neighbors based on the distance (inverse distance)
  nbrs[["nn.weights"]] = apply(nbrs$nn.dists[, -1, drop = FALSE], 1, function(x) {
    w = 1 / (x + 1e-10)
    w / sum(w)
  })

  return(nbrs)
}


# Function to perform imputation based on neighbors and weights
ImputeByNeighbors <- function(nbrs, gene_exp) {
  ncell = ncol(gene_exp)
  nn = ncol(nbrs[["nn.idx"]])

  # Construct sparse matrix representing the weighted adjacency matrix
  nbr_matrix = sparseMatrix(
    i = as.vector(t(nbrs[["nn.idx"]])),
    j = rep(1:ncell, each = nn),  # Repeated col indices (for each cell)
    x = as.vector(nbrs[["nn.weights"]]),
    dims = c(ncell, ncell)
  )

  # Perform the imputation: a weighted sum of gene exp based on neighbors
  all_exp_imp = gene_exp %*% nbr_matrix

  # Update where values are missing (assuming 0 represents missing data)
  gene_exp[gene_exp == 0] = all_exp_imp[gene_exp == 0]

  return(gene_exp)
}


#Construct network at multiplt n----
