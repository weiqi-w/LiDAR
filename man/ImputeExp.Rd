\name{ImputeExp}
\alias{ImputeExp}
\title{Impute zero ligand and receptor expression}
\description{Impute zero ligand and receptor expression.}
\usage{
ImputeExp(lidardata, nn = 30)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{nn}{Number of nearest neighbors.}
}
\value{An updated LidarData object.}
\details{The core propagation and score formulas follow the supplied method script. Database loading returns a list and does not modify the global environment.}
