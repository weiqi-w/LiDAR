MatchResponse <- function(rwrs, ligand, lr_list, receptor_exp, 
    receptorname = NULL, i = 1L) {
    mat <- rwrs[[i]]
    if (!ligand %in% colnames(mat)) 
        stop("Ligand missing from candidate response")
    response <- mat[, ligand]
    .Nonnegative(response, "Response")
    receptors <- intersect(unique(lr_list$to[lr_list$from == 
        ligand]), rownames(receptor_exp))
    if (!length(receptors)) 
        stop("No available receptor")
    rex <- .AlignCols(receptor_exp[receptors, , drop = FALSE], 
        rownames(mat), "Receptors")
    .Nonnegative(rex, "Receptor expression")
    data.frame(cell = rownames(mat), rwr = as.numeric(response), 
        receptor = as.integer(Matrix::colSums(rex) > 0))
}
SelectNetwork <- function(rwrs, n_list, ligandname, cell_data, 
    lr_list, receptor_exp, receptorname = NULL, normalize = TRUE) {
    .Flag(normalize, "normalize")
    if (!length(rwrs) || length(rwrs) != length(n_list)) 
        stop("Candidate networks are missing")
    rwrs <- lapply(rwrs, function(z) .AlignRows(z, cell_data, 
        "Candidate response"))
    selected <- list()
    info <- list()
    failed <- list()
    candidates <- list()
    for (l in ligandname) {
        nonfinite <- any(vapply(rwrs, function(z) l %in% colnames(z) && 
            any(!is.finite(z[, l])), logical(1)))
        if (nonfinite) {
            failed[[l]] <- "Non-finite candidate response"
            next
        }
        scores <- rep(NA_real_, length(rwrs))
        reasons <- rep("", length(rwrs))
        fatal <- NULL
        for (j in seq_along(rwrs)) {
            a <- tryCatch(MatchResponse(rwrs, l, lr_list, receptor_exp, 
                i = j), error = function(e) e)
            if (inherits(a, "error")) {
                fatal <- conditionMessage(a)
                break
            }
            if (sum(a$rwr) == 0) {
                reasons[j] <- "Zero response total"
                next
            }
            scores[j] <- if (stats::sd(a$rwr) == 0 || stats::sd(a$receptor) == 
                0) 
                0
            else stats::cor(a$rwr, a$receptor)
        }
        candidates[[l]] <- data.frame(n = n_list, r = scores, 
            reason = reasons)
        if (!is.null(fatal) || !any(is.finite(scores))) {
            failed[[l]] <- if (!is.null(fatal)) 
                fatal
            else "All candidate networks unavailable"
            next
        }
        best <- which.max(replace(scores, !is.finite(scores), 
            -Inf))
        v <- as.numeric(rwrs[[best]][, l])
        if (normalize) 
            v <- log1p(v/sum(v) * 10000)
        selected[[l]] <- v
        info[[l]] <- data.frame(Ligand = l, n = n_list[best], 
            r = scores[best])
    }
    result <- if (length(selected)) 
        do.call(cbind, selected)
    else .EmptyMatrix(cell_data)
    rownames(result) <- cell_data
    attr(result, "network_selection") <- if (length(info)) 
        do.call(rbind, info)
    else data.frame()
    attr(result, "candidate_results") <- candidates
    attr(result, "failures") <- failed
    result
}
.AlignRWR <- function(result, nodes, ligands) {
    .Names(nodes, "RWR graph nodes")
    .Names(ligands, "RWR seed columns")
    result <- as.matrix(result)
    if (nrow(result) != length(nodes) || ncol(result) != length(ligands)) 
        stop("dRWR returned dimensions inconsistent with graph nodes/seed columns")
    if (is.null(rownames(result))) 
        rownames(result) <- nodes
    if (is.null(colnames(result))) 
        colnames(result) <- ligands
    .Names(rownames(result), "RWR result rows")
    .Names(colnames(result), "RWR result columns")
    if (!setequal(rownames(result), nodes) || !setequal(colnames(result), 
        ligands)) 
        stop("dRWR returned names inconsistent with graph nodes/seed columns")
    result[nodes, ligands, drop = FALSE]
}
InferResponse <- function(lidardata, normalize = TRUE, parallel = FALSE) {
    .CheckObject(lidardata)
    .Need("dnet")
    .Need("igraph")
    .Flag(parallel, "parallel")
    .Flag(normalize, "normalize")
    seeds <- SelectSeeds(lidardata)
    expr <- .Expression(lidardata, FALSE)
    valid <- character()
    for (l in colnames(seeds)) {
        err <- tryCatch({
            .Nonnegative(expr$ligand[l, ], "Ligand expression")
            .Receptors(lidardata, l, expr$receptor)
            NULL
        }, error = function(e) conditionMessage(e))
        if (is.null(err)) 
            valid <- c(valid, l)
        else lidardata <- .Failure(lidardata, l, "response", 
            err)
    }
    for (l in attr(seeds, "failed_ligand")) lidardata <- .Failure(lidardata, 
        l, "seeds", "No valid positive target seeds")
    seeds <- seeds[, valid, drop = FALSE]
    if (!length(valid)) 
        stop("No ligands available for response inference")
    rwrs <- lapply(lidardata@Coembedding@nets, function(net) {
        nodes <- igraph::V(net)$name
        .Names(nodes, "RWR graph nodes")
        if (!setequal(nodes, rownames(seeds))) 
            stop("Candidate network nodes do not match seed rows")
        net_seeds <- seeds[nodes, valid, drop = FALSE]
        r <- dnet::dRWR(net, setSeeds = net_seeds, parallel = parallel, 
            verbose = FALSE, restart = 0.75)
        r <- .AlignRWR(r, nodes, colnames(net_seeds))
        r <- .AlignRows(r, lidardata@cellname, "RWR result")
        .AlignCols(r, valid, "RWR result")
    })
    result <- SelectNetwork(rwrs, lidardata@Coembedding@n_list, 
        valid, lidardata@cellname, lidardata@DB@lr_list, expr$receptor, 
        normalize = normalize)
    for (l in names(attr(result, "failures"))) lidardata <- .Failure(lidardata, 
        l, "response", attr(result, "failures")[[l]])
    lidardata@Diffusion@ligand_act <- result
    lidardata@Diffusion@ligand_diff <- .EmptyFit()
    lidardata@Diffusion@ligand_conc <- .EmptyMatrix(lidardata@cellname)
    lidardata@Lidarscore <- .ScoresEmpty(lidardata@cellname)
    lidardata
}
