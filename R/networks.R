BuildNetwork <- function(lidardata, nn = 30) {
    .Need("igraph")
    embed <- rbind(lidardata@Coembedding@cell_embedding, lidardata@Coembedding@gene_embedding)
    .Names(rownames(embed), "Network nodes")
    a <- FindNeighbors(embed, nn)
    edges <- data.frame(from = rep(rownames(embed), each = nn), 
        to = rownames(embed)[as.vector(t(a$nn.idx))])
    net <- igraph::graph_from_data_frame(edges, directed = FALSE, 
        vertices = rownames(embed))
    igraph::simplify(net, remove.multiple = TRUE, remove.loops = TRUE)
}
BuildNetworks <- function(lidardata, n = c(500, 700)) {
    .CheckObject(lidardata)
    if (!length(n) || anyDuplicated(n)) 
        stop("n must contain distinct neighbor counts")
    nets <- lapply(n, function(k) BuildNetwork(lidardata, k))
    lidardata@Coembedding@nets <- nets
    lidardata@Coembedding@n_list <- as.numeric(n)
    lidardata@Diffusion@ligand_act <- .EmptyMatrix(lidardata@cellname)
    lidardata@Diffusion@ligand_diff <- .EmptyFit()
    lidardata@Diffusion@ligand_conc <- .EmptyMatrix(lidardata@cellname)
    lidardata@Lidarscore <- .ScoresEmpty(lidardata@cellname)
    lidardata
}
SelectSeeds <- function(lidardata) {
    .Need("igraph")
    if (!length(lidardata@Coembedding@nets)) 
        stop("Run BuildNetworks first")
    nodes <- igraph::V(lidardata@Coembedding@nets[[1]])$name
    ltm <- lidardata@DB@ltm
    if (!all(rownames(ltm) %in% nodes)) 
        stop("Target genes missing from network")
    seed <- .EmptyMatrix(nodes, colnames(ltm))
    invalid <- character()
    for (l in colnames(ltm)) {
        v <- as.numeric(ltm[, l])
        names(v) <- rownames(ltm)
        if (any(!is.finite(v)) || any(v < 0) || !any(v > 0)) {
            invalid <- c(invalid, l)
            next
        }
        positive <- v[v > 0]
        threshold <- stats::quantile(positive, 0.99, names = FALSE, 
            type = 7)
        selected <- positive[positive >= threshold]
        seed[names(selected), l] <- selected
    }
    seed <- seed[, setdiff(colnames(seed), invalid), drop = FALSE]
    attr(seed, "failed_ligand") <- invalid
    seed
}
