# Load external NicheNet resources----
load_nichenet_data <- function(species = c("human", "mouse"), data_dir = NULL,
                              ligand_target_path = NULL, lr_network_path = NULL) {
  species <- match.arg(species)
  if (is.null(data_dir)) data_dir <- getOption(paste0("LiDAR.nichenet.", species))
  if (is.null(ligand_target_path) && !is.null(data_dir)) ligand_target_path <- file.path(data_dir, "ligand_target_matrix.rds")
  if (is.null(lr_network_path) && !is.null(data_dir)) lr_network_path <- file.path(data_dir, "lr_network.rds")
  if (is.null(ligand_target_path) || is.null(lr_network_path)) stop("Supply data_dir or both database file paths. Files must match the requested species.")
  if (!file.exists(ligand_target_path) || !file.exists(lr_network_path)) stop("Database file not found.")
  ltm <- readRDS(ligand_target_path)
  lr <- readRDS(lr_network_path)
  if (!(is.matrix(ltm) || methods::is(ltm, "Matrix")) || is.null(rownames(ltm)) || is.null(colnames(ltm))) stop("Ligand-target matrix must have named rows and columns.")
  if (!is.data.frame(lr) || !all(c("from", "to", "database") %in% names(lr))) stop("Network must contain from, to and database columns.")
  list(ligand_target_matrix = ltm, lr_network = lr, species = species)
}
