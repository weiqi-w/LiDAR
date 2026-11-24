\name{CalcExposure}
\alias{CalcExposure}
\title{Multiply binned ligand expression by decay weights.}
\description{Multiply binned ligand expression by decay weights.}
\usage{
CalcExposure(ligand_dist, decay_weights)
}
\arguments{
\item{ligand_dist}{Receiver-by-bin expression matrix.}
\item{decay_weights}{Numeric vector of weights, one per bin.}
}
\value{A vector or matrix of calculated values with identifiers retained where applicable.}
\seealso{\code{\link{LiDAR-package}}}
