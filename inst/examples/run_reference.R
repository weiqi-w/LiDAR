# 1. Read paths----
# Usage: Rscript run_reference.R input.rds database_directory reference.rds output_directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Expected input.rds database_directory reference.rds output_directory")
library(LiDAR)
dir.create(args[4], recursive = TRUE, showWarnings = FALSE)
seu <- readRDS(args[1])
db <- load_nichenet_data("human", data_dir = args[2])

# 2. Run analysis----
set.seed(123)
x <- PrepLidarData(seu, assay = "SCT", ligand_target_matrix = db$ligand_target_matrix, lr_network = db$lr_network)
x <- ImputeExp(x, nn = 30)
x <- CreateNetworks(x, c(500, 700))
x <- InferLigandAct(x, normalize = TRUE, parallel = FALSE)
x <- CalcLigandDist(x, dist_threshold = 10, include_self = TRUE)
x <- OptLigandDiff(x)
x <- CompLidarScore(x, mode = "aggregate")
saveRDS(x, file.path(args[4], "lidardata.rds"))

# 3. Compare named outputs----
reference <- readRDS(args[3])
compare_matrix <- function(observed, expected, label) {
  if (!identical(dimnames(observed), dimnames(expected))) stop(paste(label, "dimensions or identifiers differ"))
  observed <- as.numeric(as.matrix(observed))
  expected <- as.numeric(as.matrix(expected))
  data.frame(output = label, correlation = cor(observed, expected), max_abs_difference = max(abs(observed - expected)), rmse = sqrt(mean((observed - expected)^2)))
}
summary <- rbind(
  compare_matrix(x@Diffusion@ligand_act, reference$ligand_act, "activity"),
  compare_matrix(x@Diffusion@ligand_conc, reference$ligand_conc, "exposure"),
  compare_matrix(x@Lidarscore@aggregated_score, reference$aggregated_score, "aggregate"))
write.csv(summary, file.path(args[4], "reference_comparison.csv"), row.names = FALSE)
fit <- x@Diffusion@ligand_diff
expected <- reference$ligand_diff
comparison <- merge(expected, fit, by = "Ligand", all = TRUE, suffixes = c("_reference", "_package"))
write.csv(comparison, file.path(args[4], "diffusion_comparison.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[4], "sessionInfo.txt"))
print(summary)
