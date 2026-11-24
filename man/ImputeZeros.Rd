\name{ImputeZeros}
\alias{ImputeZeros}
\title{Replace zero entries using inverse-distance neighbor estimates.}
\description{Replace zero entries using inverse-distance neighbor estimates.}
\usage{
ImputeZeros(nbrs, gene_exp)
}
\arguments{
\item{nbrs}{Neighbor indices and distances returned by FindNeighbors.}
\item{gene_exp}{Gene-by-cell expression matrix.}
}
\value{A vector or matrix of calculated values with identifiers retained where applicable.}
\seealso{\code{\link{LiDAR-package}}}
