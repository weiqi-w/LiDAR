FindNeighbors <- function(embedding, nn = 30) {
    .Need("RANN")
    embedding <- as.matrix(embedding)
    nn <- .Integer(nn, "nn", max = nrow(embedding) - 1L)
    if (!ncol(embedding) || any(!is.finite(embedding))) 
        stop("Invalid embedding")
    a <- RANN::nn2(embedding, k = nn + 1L)
    idx <- matrix(0L, nrow(embedding), nn)
    ds <- matrix(0, nrow(embedding), nn)
    for (i in seq_len(nrow(embedding))) {
        keep <- which(a$nn.idx[i, ] != i)[seq_len(nn)]
        idx[i, ] <- a$nn.idx[i, keep]
        ds[i, ] <- a$nn.dists[i, keep]
    }
    w <- 1/(ds + 1e-10)
    w <- w/rowSums(w)
    list(nn.idx = idx, nn.dists = ds, nn.weights = t(w), cells = rownames(embedding))
}
ImputeZeros <- function(nbrs, gene_exp) {
    gene_exp <- .AlignCols(gene_exp, nbrs$cells, "Expression for imputation")
    ncell <- ncol(gene_exp)
    nn <- ncol(nbrs$nn.idx)
    if (!identical(dim(nbrs$nn.weights), c(nn, ncell))) 
        stop("Neighbor weight dimensions mismatch")
    W <- Matrix::sparseMatrix(i = as.vector(t(nbrs$nn.idx)), 
        j = rep(seq_len(ncell), each = nn), x = as.vector(nbrs$nn.weights), 
        dims = c(ncell, ncell))
    good <- vapply(seq_len(nrow(gene_exp)), function(i) all(is.finite(gene_exp[i, 
        ])) && all(gene_exp[i, ] >= 0), logical(1))
    if (any(good)) {
        part <- gene_exp[good, , drop = FALSE]
        imputed <- part %*% W
        zero <- part == 0
        part[zero] <- imputed[zero]
        gene_exp[good, ] <- part
    }
    gene_exp
}
ImputeExpr <- function(lidardata, nn = 30) {
    .CheckObject(lidardata)
    ce <- .AlignRows(lidardata@Coembedding@cell_embedding, lidardata@cellname, 
        "Embedding")
    nbrs <- FindNeighbors(ce, nn)
    lidardata@ImputationData@ligand_exp <- ImputeZeros(nbrs, 
        lidardata@Data@ligand_exp)
    lidardata@ImputationData@receptor_exp <- ImputeZeros(nbrs, 
        lidardata@Data@receptor_exp)
    lidardata@Diffusion@ligand_dist <- array(0, c(0, 0, 0))
    lidardata@Diffusion@ligand_diff <- .EmptyFit()
    lidardata@Diffusion@ligand_conc <- .EmptyMatrix(lidardata@cellname)
    lidardata@Lidarscore <- .ScoresEmpty(lidardata@cellname)
    lidardata
}
