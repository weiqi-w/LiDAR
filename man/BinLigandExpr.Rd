\name{BinLigandExpr}
\alias{BinLigandExpr}
\title{Sum source expression into distance bins for each receiver.}
\description{Sum source expression into distance bins for each receiver.}
\usage{
BinLigandExpr(lidardata, bin = NULL, n_bins = 10L, imp = TRUE)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{bin}{Positive bin width in unchanged input coordinate units; NULL uses minimum positive distance.}
\item{n_bins}{Number of distance bins before the terminal zero-weight bin. A total of n_bins + 1 bins are constructed.}
\item{imp}{Whether to use imputed expression.}
}
\value{An updated LidarData object.}
\seealso{\code{\link{LiDAR-package}}}
