\name{FitDecay}
\alias{FitDecay}
\title{Optimize decay parameters against receiver response.}
\description{Optimize decay parameters against receiver response.}
\usage{
FitDecay(response, ligand_dist, distance, B, fit_timeout = 300)
}
\arguments{
\item{response}{Numeric response vector aligned to potential receivers.}
\item{ligand_dist}{Receiver-by-bin expression matrix.}
\item{distance}{Bin left-edge distances aligned to columns of ligand_dist.}
\item{B}{Maximum decay distance in input coordinate units.}
\item{fit_timeout}{Cooperative time limit in seconds per fit; defaults to 300.}
}
\value{A named list of results.}
\seealso{\code{\link{LiDAR-package}}}
