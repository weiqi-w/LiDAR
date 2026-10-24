PrepLidarData <- function(data, assay = "SCT", ligand_target_matrix, lr_network,
                          coord = NULL) {

  if (!assay %in% names(data@assays)) {
    stop(paste("Error: Assay named", assay, "does not exist in the data object."))
  }

  message("Run MCA")
  if (!requireNamespace("CelliD", quietly = TRUE)) stop("Install CelliD to run MCA.")
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Install Seurat to prepare spatial data.")
  DefaultAssay(data) <- assay
  data = CelliD::RunMCA(data, features = rownames(data))
  message(paste0("Coembed ", nrow(data@reductions$mca@cell.embeddings),
                 " cells with ", nrow(data@reductions$mca@feature.loadings), " features"))

  # Extract gene and cell names
  DefaultAssay(data) = assay
  gene_data <- rownames(data@reductions$mca@feature.loadings)
  cell_data <- colnames(data)
  ngene <- length(gene_data)
  ncell <- length(cell_data)

  # Process NicheNet to create ligand target matrix (ltm)
  ltm <- ligand_target_matrix[rownames(ligand_target_matrix) %in% gene_data,
                              colnames(ligand_target_matrix) %in% gene_data, drop = FALSE]
  ltm <- ltm[, colSums(ltm) > 0, drop = FALSE]
  ligandname <- colnames(ltm)

  # Filter the ligand-receptor network
  lr_list <- lr_network %>%
    filter((from %in% gene_data) & (from %in% ligandname) & (to %in% gene_data) & !str_detect(database, "ppi"))
  ligand_lr_list <- unique(lr_list$from)
  receptor_lr_list <- unique(lr_list$to)

  # Determine receptor names and the number of receptors
  receptorname <- receptor_lr_list
  nreceptor <- length(receptorname)

  # Adjust ligand names (remove ligands without receptors)
  ligandname <- ligandname[ligandname %in% ligand_lr_list]
  nligand <- length(ligandname)

  # Filter ltm
  ltm = ltm[, ligandname, drop = FALSE]


  # Extract ligand and receptor expression data from SCT assay
  exp_data <- SeuratObject::GetAssayData(data, assay = assay, layer = "data")
  ligand_exp <- exp_data[match(ligandname, rownames(exp_data)), , drop = FALSE]
  receptor_exp <- exp_data[match(receptorname, rownames(exp_data)), , drop = FALSE]

  if (!length(ligandname)) stop("No supported ligands remain after database filtering.")
  # Create layers of the appropriate classes
  data_layer <- methods::new("Data", data = exp_data,
                             ligand_exp = ligand_exp,
                             receptor_exp = receptor_exp)
  imputation_data_layer <- methods::new("ImputationData",
                                        data = matrix(NA, nrow = 0, ncol = 0),
                                        ligand_exp = matrix(NA, nrow = 0, ncol = 0),
                                        receptor_exp = matrix(NA, nrow = 0, ncol = 0))

  # position
  if (!is.null(coord)) {
    pos <- coord
  } else {
    pos <- Seurat::GetTissueCoordinates(data)
  }

  pos <- as.data.frame(pos)
  if (all(c("imagerow", "imagecol") %in% names(pos))) {
    pos <- pos[, c("imagerow", "imagecol"), drop = FALSE]
  } else if (all(c("x", "y") %in% names(pos))) {
    pos <- pos[, c("x", "y"), drop = FALSE]
  }
  if (ncol(pos) != 2L || !all(vapply(pos, is.numeric, logical(1)))) stop("Coordinates must have two numeric columns.")
  if (is.null(rownames(pos)) || !setequal(rownames(pos), cell_data)) stop("Coordinate row names must match cell names.")
  pos <- pos[cell_data, , drop = FALSE]
  if (any(!is.finite(as.matrix(pos)))) stop("Coordinates must be finite.")
  # dist_pos = as.matrix(dist(pos)) %>% round(3)
  # if (!is.null(distance_bin)) {
  #   bin <- distance_bin
  # } else {
  #   bin <- min(dist_pos[dist_pos > 0])
  # }
  # dist_pos_bin_all = dist_pos%/%bin + 1
  # #dist_pos_bin[dist_pos_bin>dist_threshold] = 0
  #
  # dist_pos_bin_all = Matrix::Matrix(dist_pos_bin_all, sparse = T)
  image_layer <- methods::new("Image",
                              pos = pos,
                              dist = new("dgCMatrix", i = integer(0), p = 0L, Dim = c(0L, 0L), x = numeric(0)),
                              dist_bin_all = new("dgCMatrix", i = integer(0), p = 0L, Dim = c(0L, 0L), x = numeric(0)))


  # Create a coembedding object (assuming cell_embedding and gene_embedding need to be generated here)
  cell_embed = data@reductions$mca@cell.embeddings
  cell_embedding <- matrix(cell_embed,
                           nrow = dim(cell_embed)[1],
                           ncol = dim(cell_embed)[2],
                           dimnames = list(rownames(cell_embed), colnames(cell_embed)))
  gene_embed = data@reductions$mca@feature.loadings
  gene_embedding <- matrix(gene_embed,
                           nrow = dim(gene_embed)[1],
                           ncol = dim(gene_embed)[2],
                           dimnames = list(rownames(gene_embed), colnames(gene_embed)))
  nets = list()
  n_list = numeric()

  coembedding_layer <- methods::new("Coembedding",
                                    cell_embedding = cell_embedding,
                                    gene_embedding = gene_embedding,
                                    nets = nets,
                                    n_list = n_list)

  # Create a distance object (assuming ligand_dist needs some computation)
  ligand_act <- matrix(NA, ncol = 0, nrow = 0)
  ligand_dist <- array(NA, dim = c(0,0,0))  # Placeholder for actual computation
  ligand_diff <- data.frame()
  ligand_conc <- matrix(NA, ncol = 0, nrow = 0)
  dist_threshold <- numeric()
  bin = numeric()
  include_self = NA
  failed_ligand <- character()

  diffusion_layer <- methods::new("Diffusion",
                                  ligand_act = ligand_act,
                                  ligand_dist = ligand_dist,
                                  ligand_diff = ligand_diff,
                                  ligand_conc = ligand_conc,
                                  dist_threshold = dist_threshold,
                                  bin = bin,
                                  include_self = include_self,
                                  failed_ligand = failed_ligand)

  aggregated_score <- matrix()
  pairwise_score <- list()
  sender_score <- matrix()
  receiver_score <- matrix()

  lidarscore_layer <- methods::new("Lidarscore",
                                   aggregated_score = aggregated_score,
                                   pairwise_score = pairwise_score,
                                   sender_score = sender_score,
                                   receiver_score = receiver_score)

  # Create a database object
  db_layer <- methods::new("DB",
                           ltm = ltm,
                           lr_list = lr_list,
                           ligandname = ligandname,
                           receptorname = receptorname,
                           nligand = nligand,
                           nreceptor = nreceptor)

  # Create the LidarData layer
  lidar_data <- methods::new("LidarData",
                             Data = data_layer,
                             ImputationData = imputation_data_layer,
                             Image = image_layer,
                             Coembedding = coembedding_layer,
                             Diffusion = diffusion_layer,
                             Lidarscore = lidarscore_layer,
                             DB = db_layer,
                             cellname = cell_data,
                             genename = gene_data,
                             ncell = ncell,
                             ngene = ngene)

  return(lidar_data)
}

#Impute ligand and receptor exp----
