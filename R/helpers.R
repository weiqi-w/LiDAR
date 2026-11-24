.Integer <- function(x, name, min = 1L, max = .Machine$integer.max) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || 
        x < min || x > max || x != floor(x)) 
        stop(name, " must be an integer in [", min, ", ", max, 
            "]")
    as.integer(x)
}
.Flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) 
        stop(name, " must be TRUE or FALSE")
    x
}
.Names <- function(x, label) {
    if (is.null(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) 
        stop(label, " must have unique, nonempty names")
    invisible(x)
}
.Nonnegative <- function(x, label) {
    if (any(!is.finite(x)) || any(x < 0)) 
        stop(label, " contains non-finite or negative values")
    invisible(x)
}
.AlignRows <- function(x, ids, label) {
    .Names(rownames(x), paste(label, "rows"))
    if (!all(ids %in% rownames(x))) 
        stop(label, " is missing required rows")
    x[ids, , drop = FALSE]
}
.AlignCols <- function(x, ids, label) {
    .Names(colnames(x), paste(label, "columns"))
    if (!all(ids %in% colnames(x))) 
        stop(label, " is missing required columns")
    x[, ids, drop = FALSE]
}
.EmptyMatrix <- function(rows = character(), cols = character()) matrix(0, 
    length(rows), length(cols), dimnames = list(rows, cols))
.EmptyFit <- function() data.frame(Ligand = character(), ki = numeric(), 
    n = numeric(), coupling = numeric(), p.value = numeric(), r50 = numeric(), 
    r01 = numeric(), convergence = integer(), n_receivers = integer(), 
    stringsAsFactors = FALSE)
.ScoresEmpty <- function(cells) methods::new("Lidarscore", aggregated_score = .EmptyMatrix(cells), 
    pairwise_score = list(), sender_score = .EmptyMatrix(character(), 
        cells), receiver_score = .EmptyMatrix(character(), cells))
.CheckObject <- function(x) {
    if (!methods::is(x, "LidarData")) 
        stop("Expected a LidarData object")
    methods::validObject(x)
    .Names(x@cellname, "Spatial units")
    invisible(x)
}
.Failure <- function(x, ligand, stage, reason) {
    entry <- data.frame(Ligand = ligand, stage = stage, reason = reason, 
        stringsAsFactors = FALSE)
    attr(x, "failure_log") <- rbind(attr(x, "failure_log"), entry)
    x@Diffusion@failed_ligand <- unique(c(x@Diffusion@failed_ligand, 
        ligand))
    message(ligand, " [", stage, "]: ", reason)
    x
}
.Expression <- function(x, imp) {
    .Flag(imp, "imp")
    layer <- if (imp) 
        x@ImputationData
    else x@Data
    if (!nrow(layer@ligand_exp) || !nrow(layer@receptor_exp)) 
        stop("Expression is empty; run ImputeExpr first or use imp=FALSE")
    list(ligand = .AlignCols(layer@ligand_exp, x@cellname, "Ligand expression"), 
        receptor = .AlignCols(layer@receptor_exp, x@cellname, 
            "Receptor expression"))
}
.Receptors <- function(x, ligand, expression) {
    ids <- unique(x@DB@lr_list$to[x@DB@lr_list$from == ligand])
    ids <- intersect(ids, rownames(expression))
    if (!length(ids)) 
        stop("No available receptor")
    values <- expression[ids, , drop = FALSE]
    .Nonnegative(values, "Receptor expression")
    setNames(as.numeric(Matrix::colSums(values)), colnames(values))
}
.Coords <- function(pos, cells, coord_cols = NULL, cell_col = NULL) {
    pos <- as.data.frame(pos)
    if (is.null(cell_col)) {
        candidate <- intersect(c("cell", "cells", "barcode", 
            "barcodes", "id", "ids"), names(pos))
        if (length(candidate)) 
            cell_col <- candidate[1L]
    }
    ids <- if (!is.null(cell_col)) {
        if (length(cell_col) != 1L || !cell_col %in% names(pos)) 
            stop("Invalid cell_col")
        as.character(pos[[cell_col]])
    }
    else rownames(pos)
    .Names(ids, "Coordinate rows")
    if (!all(cells %in% ids)) 
        stop("Coordinate cell names do not match expression; provide cell_col or row names")
    pairs <- list(c("X", "Y"), c("x", "y"), c("col", "row"), 
        c("imagecol", "imagerow"), c("imagecl", "imagerow"), 
        c("pxl_col_in_fullres", "pxl_row_in_fullres"), c("spatial_x", 
            "spatial_y"))
    if (is.null(coord_cols)) {
        found <- Filter(function(z) all(z %in% names(pos)), pairs)
        if (!length(found)) 
            stop("Cannot identify coordinates; supply coord_cols=c('x_column','y_column')")
        coord_cols <- found[[1L]]
        message("Coordinate columns: ", paste(coord_cols, collapse = ", "))
    }
    if (length(coord_cols) != 2L || anyDuplicated(coord_cols) || 
        !all(coord_cols %in% names(pos))) 
        stop("coord_cols must name two distinct columns, in x/y order")
    out <- pos[match(cells, ids), coord_cols, drop = FALSE]
    if (!all(vapply(out, is.numeric, logical(1))) || any(!is.finite(as.matrix(out)))) 
        stop("Coordinates must be finite numeric values")
    names(out) <- c("x", "y")
    rownames(out) <- cells
    attr(out, "source_columns") <- coord_cols
    out
}
