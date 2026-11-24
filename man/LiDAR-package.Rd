\name{LiDAR-package}
\alias{LiDAR}
\alias{LiDAR-package}
\docType{package}
\title{Spatial ligand propagation and communication}
\description{Models ligand exposure and downstream responses using user-supplied prior networks.}
\details{Use CreateLiDAR, ImputeExpr, BuildNetworks, InferResponse, BinLigandExpr, FitPropagation and ComputeScores in sequence. Distances retain input units. The first bin has weight one; the final bin has weight zero. Continuous ranges are derived from the fitted curve. Optional permutation tests reoptimize decay after each response permutation. Ordinary Pearson P values after fitting are descriptive and do not account for parameter selection. Empirical P values are reported separately and do not change default score filtering.}
\examples{
CalcRanges(ki = 0.2, n = 2, B = 500)
ScaleMax(c(0, 0, 0))
}
