ScaleMax <- function(x) {
    .Nonnegative(x, "Score component")
    if (!length(x)) 
        return(x)
    maximum <- max(x)
    if (maximum == 0) 
        return(x * 0)
    x/maximum
}
ComputeScores <- function(lidardata, mode = c("aggregate", "pairwise"), 
    ligand_list = NULL, pval = 0.001, imp = TRUE, coupling = 0, Log = FALSE) {
    .CheckObject(lidardata)
    mode <- match.arg(mode)
    .Flag(Log, "Log")
    if (!is.null(coupling) && (!is.numeric(coupling) || length(coupling) != 1L ||
        !is.finite(coupling) || abs(coupling) > 1))
        stop("coupling must be NULL or a number between -1 and 1")
    expr <- .Expression(lidardata, imp)
    cells <- lidardata@cellname
    df <- lidardata@Diffusion@ligand_diff
    cfg <- attr(df, "fit_config")
    if (is.null(cfg)) 
        stop("Run FitPropagation first")
    if (!identical(cfg$imp, imp)) 
        stop("imp differs from fitted exposure; refit with matching imp")
    eligible <- GetDecay(lidardata, pval = pval, p_column = "p.value")
    if (!is.null(coupling))
        eligible <- eligible[is.finite(eligible$coupling) & eligible$coupling > coupling, , drop = FALSE]
    selected <- eligible$Ligand
    if (!is.null(ligand_list))
        selected <- intersect(as.character(ligand_list), selected)
    response <- .AlignRows(lidardata@Diffusion@ligand_act, cells, 
        "Response")
    aggregate <- list()
    pair <- list()
    send <- list()
    receive <- list()
    if (mode == "pairwise") {
        bins <- .AlignCols(.AlignRows(lidardata@Image@dist_bin_all, 
            cells, "Distance bins"), cells, "Distance bins")
        if (lidardata@Diffusion@bin != cfg$bin || lidardata@Diffusion@dist_threshold != 
            cfg$n_bins) 
            stop("Distance bins changed after fitting")
    }
    for (l in selected) {
        value <- tryCatch({
            if (!l %in% df$Ligand || !l %in% colnames(response) || 
                !l %in% rownames(expr$ligand)) 
                stop("Ligand lacks valid fitted inputs")
            availability <- .Receptors(lidardata, l, expr$receptor)
            ligand_expr <- as.numeric(expr$ligand[l, ])
            .Nonnegative(ligand_expr, "Ligand expression")
            act <- ScaleMax(as.numeric(response[, l])) * 10
            rec <- ScaleMax(availability) * 10
            if (mode == "aggregate") {
                exposure <- .AlignRows(lidardata@Diffusion@ligand_conc, 
                  cells, "Exposure")
                if (!l %in% colnames(exposure)) 
                  stop("Exposure missing; fit with CalcLigandConc=TRUE")
                ScaleMax(as.numeric(exposure[, l])) * 10 * act * 
                  rec
            }
            else {
                fitrow <- df[match(l, df$Ligand), , drop = FALSE]
                w <- .DecayWeights(c(fitrow$ki, fitrow$n), cfg$distance, 
                  cfg$B)
                weight <- bins
                weight@x <- w[as.integer(weight@x)]
                weight <- Matrix::drop0(weight)
                Matrix::Diagonal(x = ScaleMax(ligand_expr) * 
                  10) %*% weight %*% Matrix::Diagonal(x = rec * 
                  act)
            }
        }, error = function(e) e)
        if (inherits(value, "error")) {
            lidardata <- .Failure(lidardata, l, "score", conditionMessage(value))
            next
        }
        if (mode == "aggregate") 
            aggregate[[l]] <- if (Log) log1p(value) / log(10) else value
        else {
            dimnames(value) <- list(cells, cells)
            send[[l]] <- as.numeric(Matrix::rowSums(value))
            receive[[l]] <- as.numeric(Matrix::colSums(value))
            if (Log) {
                value@x <- log1p(value@x) / log(10)
                send[[l]] <- log1p(send[[l]]) / log(10)
                receive[[l]] <- log1p(receive[[l]]) / log(10)
            }
            pair[[l]] <- value
        }
    }
    if (mode == "aggregate") {
        m <- if (length(aggregate)) 
            do.call(cbind, aggregate)
        else .EmptyMatrix(cells)
        rownames(m) <- cells
        lidardata@Lidarscore@aggregated_score <- Matrix::Matrix(m, 
            sparse = TRUE)
    }
    else {
        sender <- if (length(send)) 
            do.call(rbind, send)
        else .EmptyMatrix(character(), cells)
        receiver <- if (length(receive)) 
            do.call(rbind, receive)
        else .EmptyMatrix(character(), cells)
        colnames(sender) <- colnames(receiver) <- cells
        lidardata@Lidarscore@pairwise_score <- pair
        lidardata@Lidarscore@sender_score <- Matrix::Matrix(sender, 
            sparse = TRUE)
        lidardata@Lidarscore@receiver_score <- Matrix::Matrix(receiver, 
            sparse = TRUE)
    }
    lidardata
}
LoadPriors <- function(ligand_target_path, lr_network_path) {
    if (!file.exists(ligand_target_path) || !file.exists(lr_network_path)) 
        stop("Prior data file does not exist")
    list(ligand_target_matrix = readRDS(ligand_target_path), 
        lr_network = readRDS(lr_network_path))
}
