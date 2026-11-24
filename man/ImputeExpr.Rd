\name{ImputeExpr}
\alias{ImputeExpr}
\title{Impute only zero ligand and receptor values using embedding neighbors.}
\description{Impute only zero ligand and receptor values using embedding neighbors.}
\usage{
ImputeExpr(lidardata, nn = 30)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{nn}{Number of neighbors, excluding self.}
}
\value{An updated LidarData object.}
\seealso{\code{\link{LiDAR-package}}}
