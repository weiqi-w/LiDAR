test_that("distance bins respect self contribution and cell order", {
  x <- make_fixture()
  expect_true(methods::validObject(x))
  y <- CalcLigandDist(x, dist_threshold = 5, distance_bin = 1)
  expect_equal(y@Diffusion@ligand_dist[, , 1], as.matrix(x@Data@ligand_exp), ignore_attr = TRUE)
  expect_identical(dimnames(y@Diffusion@ligand_dist)$Cell, x@cellname)
  z <- CalcLigandDist(x, dist_threshold = 5, distance_bin = 1, include_self = FALSE)
  expect_true(all(z@Diffusion@ligand_dist[, , 1] == 0))
  expect_error(CalcLigandDist(x, distance_bin = 0), "positive")
  expect_error(CalcLigandDist(x, RealDist = TRUE), "not implemented")
})

test_that("imputation preserves observed expression", {
  x <- make_fixture()
  x@Data@ligand_exp[1, 2] <- 0
  y <- ImputeExp(x, nn = 3)
  nz <- x@Data@ligand_exp != 0
  expect_equal(y@ImputationData@ligand_exp[nz], x@Data@ligand_exp[nz])
  expect_gt(y@ImputationData@ligand_exp[1, 2], 0)
  expect_error(ImputeExp(x, nn = 40), "nn must")
  expect_s4_class(ImputeExp(x, nn = 1), "LidarData")
})

test_that("graphs and optional RWR preserve cell identifiers", {
  x <- make_fixture()
  rownames(x@Coembedding@cell_embedding) <- x@cellname
  x <- CreateNetworks(x, c(3, 5))
  expect_length(x@Coembedding@nets, 2)
  expect_setequal(igraph::V(x@Coembedding@nets[[1]])$name, c(x@cellname, x@genename))
  expect_true(igraph::is_connected(x@Coembedding@nets[[1]]))
  skip_if_not_installed("dnet")
  y <- InferLigandAct(x, parallel = FALSE)
  expect_identical(rownames(y@Diffusion@ligand_act), x@cellname)
  expect_true(all(is.finite(y@Diffusion@ligand_act)))
  expect_true(all(colSums(y@Diffusion@ligand_act) > 0))
})

test_that("optimization recovers a known spatial response", {
  x <- fitted_fixture()
  expect_equal(nrow(x@Diffusion@ligand_diff), 2)
  expect_true(all(x@Diffusion@ligand_diff$r > 0.99))
  expect_true(all(is.finite(x@Diffusion@ligand_conc)))
  y <- CompLidarScore(x, mode = "aggregate", ligand_list = "L1")
  expect_equal(dim(y@Lidarscore@aggregated_score), c(40L, 1L))
  expect_true(all(is.finite(y@Lidarscore@aggregated_score)))
  z <- CompLidarScore(x, mode = "pairwise", ligand_list = "L1")
  expect_equal(as.numeric(z@Lidarscore@sender_score["L1", ]), as.numeric(Matrix::rowSums(z@Lidarscore@pairwise_score$L1)))
  expect_equal(as.numeric(z@Lidarscore@receiver_score["L1", ]), as.numeric(Matrix::colSums(z@Lidarscore@pairwise_score$L1)))
})

test_that("database loading returns a list without global assignments", {
  d <- tempfile(); dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  prior <- matrix(1, 2, 1, dimnames = list(c("T1", "T2"), "L1"))
  saveRDS(prior, file.path(d, "ligand_target_matrix.rds"))
  saveRDS(data.frame(from = "L1", to = "R1", database = "toy"), file.path(d, "lr_network.rds"))
  old <- exists("ligand_target_matrix", .GlobalEnv, inherits = FALSE)
  db <- load_nichenet_data("human", data_dir = d)
  expect_equal(db$ligand_target_matrix, prior)
  expect_identical(exists("ligand_target_matrix", .GlobalEnv, inherits = FALSE), old)
  expect_error(load_nichenet_data(data_dir = file.path(d, "missing")), "not found")
})

test_that("empty fit results remain valid S4 slot values", {
  x <- CalcLigandDist(make_fixture(), distance_bin = 1)
  x@Diffusion@ligand_act <- matrix(0, x@ncell, 2, dimnames = list(x@cellname, x@DB@ligandname))
  x <- suppressMessages(OptLigandDiff(x))
  expect_equal(nrow(x@Diffusion@ligand_diff), 0)
  expect_equal(ncol(x@Diffusion@ligand_conc), 0)
  expect_setequal(x@Diffusion@failed_ligand, c("L1", "L2"))
  expect_true(methods::validObject(x))
  expect_error(CompLidarScore(x, mode = "aggregate"), "No ligands")
})


test_that("zero cell activity is preserved without division by zero", {
  x <- make_fixture()
  raw <- matrix(0, x@ncell, 2, dimnames = list(x@cellname, x@DB@ligandname))
  call_select <- function(raw, normalize = TRUE) {
    LiDAR:::SelectBestNN(list(raw), 3, x@DB@ligandname, x@cellname,
      x@DB@lr_list, x@Data@receptor_exp, x@DB@receptorname, normalize)
  }
  expect_warning(y <- call_select(raw), "No RWR activity reached cells")
  expect_true(all(is.finite(y)))
  expect_true(all(y == 0))
  raw[,] <- seq_len(x@ncell)
  expect_equal(call_select(raw), log1p(sweep(raw, 2, colSums(raw), "/") * 10000))
  raw[1, 1] <- NaN
  expect_error(call_select(raw), "finite, nonnegative")
})

test_that("disconnected gene seeds yield finite zero cell activity", {
  skip_if_not_installed("dnet")
  x <- make_fixture()
  x@Coembedding@gene_embedding <- matrix(seq_len(12), 6, 2,
    dimnames = list(x@genename, c("a", "b")))
  x <- CreateNetworks(x, c(3, 5))
  expect_warning(y <- InferLigandAct(x, parallel = FALSE), "No RWR activity reached cells")
  expect_true(all(is.finite(y@Diffusion@ligand_act)))
  expect_true(all(y@Diffusion@ligand_act == 0))
})
