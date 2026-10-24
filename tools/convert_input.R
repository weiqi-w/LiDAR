# Usage: Rscript tools/convert_input.R input.qs output.rds
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Expected input.qs output.rds")
if (!requireNamespace("qs", quietly = TRUE)) stop("Install qs for this one-time conversion.")
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Install Seurat to load the input classes.")
x <- qs::qread(args[1])
saveRDS(x, args[2])
