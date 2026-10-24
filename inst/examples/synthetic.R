# Synthetic downstream demonstration; no MCA or RWR is run here.
local({
  library(LiDAR)
make_fixture <- function(ncell = 40L) {
  cells <- paste0("cell", seq_len(ncell))
  ligands <- c("L1", "L2")
  receptors <- c("R1", "R2")
  genes <- c(ligands, receptors, "T1", "T2")
  j <- seq_len(ncell)
  expression <- rbind((sin(j / 4) + 1.2), (cos(j / 5) + 1.3), rep(1, ncell), rep(1, ncell), j / ncell, rev(j) / ncell)
  dimnames(expression) <- list(genes, cells)
  expression <- Matrix::Matrix(expression, sparse = TRUE)
  pos <- data.frame(x = j, y = rep(0, ncell), row.names = cells)
  empty <- matrix(numeric(), 0, 0)
  prior <- matrix(1, length(genes), 2, dimnames = list(genes, ligands))
  methods::new("LidarData",
    Data = methods::new("Data", data = expression, ligand_exp = expression[ligands, ], receptor_exp = expression[receptors, ]),
    ImputationData = methods::new("ImputationData", data = empty, ligand_exp = expression[ligands, ], receptor_exp = expression[receptors, ]),
    Image = methods::new("Image", pos = pos, dist = empty, dist_bin_all = Matrix::Matrix(empty, sparse = TRUE)),
    Coembedding = methods::new("Coembedding", cell_embedding = matrix(c(j, sin(j)), ncell, 2, dimnames = list(cells, c("a", "b"))), gene_embedding = matrix(seq_len(12), 6, 2, dimnames = list(genes, c("a", "b"))), nets = list(), n_list = numeric()),
    Diffusion = methods::new("Diffusion", ligand_act = empty, ligand_dist = array(numeric(), c(0, 0, 0)), ligand_diff = data.frame(), ligand_conc = empty, dist_threshold = numeric(), bin = numeric(), include_self = TRUE, failed_ligand = character()),
    Lidarscore = methods::new("Lidarscore", aggregated_score = empty, pairwise_score = list(), sender_score = empty, receiver_score = empty),
    DB = methods::new("DB", ltm = prior, lr_list = data.frame(from = ligands, to = receptors, database = "example"), ligandname = ligands, receptorname = receptors, nligand = 2L, nreceptor = 2L),
    cellname = cells, genename = genes, ncell = as.numeric(ncell), ngene = as.numeric(length(genes)))
}

fitted_fixture <- function() {
  x <- CalcLigandDist(make_fixture(), dist_threshold = 5, distance_bin = 1)
  d <- 0:5
  kernel <- 1 / (1 + ((d / 5) / 0.5)^2)
  kernel <- (kernel - min(kernel)) / (max(kernel) - min(kernel))
  activity <- sapply(x@DB@ligandname, function(lg) as.vector(x@Diffusion@ligand_dist[lg, , ] %*% kernel))
  rownames(activity) <- x@cellname
  x@Diffusion@ligand_act <- activity
  OptLigandDiff(x)
}

  x <- fitted_fixture()
  x <- CompLidarScore(x, mode = "aggregate", ligand_list = c("L1", "L2"))
  print(x@Diffusion@ligand_diff)
  print(dim(x@Lidarscore@aggregated_score))
})
