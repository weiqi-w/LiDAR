CreateLiDAR <- function(data, assay = "SCT", ligand_target_matrix, 
    lr_network, coord = NULL, coord_cols = NULL, cell_col = NULL, 
    layer = "data") {
    .Need("Seurat")
    .Need("SeuratObject")
    .Need("CelliD")
    if (!assay %in% names(data@assays)) 
        stop("Assay does not exist: ", assay)
    SeuratObject::DefaultAssay(data) <- assay
    exp_data <- if (utils::packageVersion("SeuratObject") >= 
        "5.0.0") 
        SeuratObject::GetAssayData(data, assay = assay, layer = layer)
    else SeuratObject::GetAssayData(data, assay = assay, slot = layer)
    .Names(rownames(exp_data), "Gene expression rows")
    .Names(colnames(exp_data), "Gene expression columns")
    if (!nrow(exp_data) || !ncol(exp_data)) 
        stop("Expression layer is empty")
    message("Run MCA on assay ", assay)
    # Use the selected expression layer directly for the joint MCA embedding.
    mca <- CelliD::RunMCA(as.matrix(exp_data), features = rownames(exp_data))
    ce <- as.matrix(mca$cellsCoordinates)
    ge <- as.matrix(mca$featuresCoordinates)
    cells <- colnames(exp_data)
    genes <- rownames(ge)
    ce <- .AlignRows(ce, cells, "Cell embedding")
    .Names(genes, "Gene embedding rows")
    if (length(intersect(cells, genes))) 
        stop("Gene and spatial-unit names overlap")
    if (ncol(ce) != ncol(ge) || any(!is.finite(ce)) || any(!is.finite(ge))) 
        stop("Invalid MCA embedding")
    .Names(rownames(ligand_target_matrix), "Target matrix rows")
    .Names(colnames(ligand_target_matrix), "Target matrix columns")
    required <- c("from", "to", "database")
    if (!all(required %in% names(lr_network))) 
        stop("lr_network needs from, to, database")
    net <- as.data.frame(lr_network, stringsAsFactors = FALSE)
    net$from <- as.character(net$from)
    net$to <- as.character(net$to)
    present <- intersect(genes, rownames(exp_data))
    keep <- !is.na(net$from) & !is.na(net$to) & !is.na(net$database) & 
        net$from %in% present & net$to %in% present & !grepl("ppi", 
        net$database)
    net <- unique(net[keep, , drop = FALSE])
    ligands <- intersect(colnames(ligand_target_matrix), unique(net$from))
    targets <- intersect(rownames(ligand_target_matrix), genes)
    ltm <- ligand_target_matrix[targets, ligands, drop = FALSE]
    good <- vapply(seq_len(ncol(ltm)), function(j) all(is.finite(ltm[, 
        j])) && all(ltm[, j] >= 0) && any(ltm[, j] > 0), logical(1))
    excluded <- ligands[!good]
    ligands <- ligands[good]
    ltm <- ltm[, ligands, drop = FALSE]
    if (!length(ligands)) 
        stop("No ligands with usable receptors and positive target weights")
    net <- net[net$from %in% ligands, , drop = FALSE]
    receptors <- unique(net$to)
    if (is.null(coord)) 
        coord <- Seurat::GetTissueCoordinates(data)
    pos <- .Coords(coord, cells, coord_cols, cell_col)
    z <- .EmptyMatrix()
    x <- methods::new("LidarData", Data = methods::new("Data", 
        data = exp_data, ligand_exp = exp_data[ligands, , drop = FALSE], 
        receptor_exp = exp_data[receptors, , drop = FALSE]), 
        ImputationData = methods::new("ImputationData", data = z, 
            ligand_exp = z, receptor_exp = z), Image = methods::new("Image", 
            pos = pos, dist = z, dist_bin_all = Matrix::Matrix(z, 
                sparse = TRUE)), Coembedding = methods::new("Coembedding", 
            cell_embedding = ce, gene_embedding = ge, nets = list(), 
            n_list = numeric()), Diffusion = methods::new("Diffusion", 
            ligand_act = .EmptyMatrix(cells), ligand_dist = array(0, 
                c(0, 0, 0)), ligand_diff = .EmptyFit(), ligand_conc = .EmptyMatrix(cells), 
            dist_threshold = numeric(), bin = numeric(), failed_ligand = character()), 
        Lidarscore = .ScoresEmpty(cells), DB = methods::new("DB", 
            ltm = ltm, lr_list = net, ligandname = ligands, receptorname = receptors, 
            nligand = as.integer(length(ligands)), nreceptor = as.integer(length(receptors))), 
        cellname = cells, genename = genes, ncell = as.numeric(length(cells)), 
        ngene = as.numeric(length(genes)))
    attr(x, "input_config") <- list(assay = assay, layer = layer, 
        coordinate_units = "input units")
    for (l in excluded) x <- .Failure(x, l, "prepare", "Invalid or absent positive target weights")
    x
}
