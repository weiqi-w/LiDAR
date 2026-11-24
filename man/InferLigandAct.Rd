\name{InferLigandAct}
\alias{InferLigandAct}
\title{Infer ligand responses by random walk and select a candidate network.}
\description{Infer ligand responses by random walk and select a candidate network.}
\usage{
InferLigandAct(lidardata, normalize = TRUE, parallel = FALSE)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{normalize}{Whether to scale selected response vectors by their maximum.}
\item{parallel}{Whether dnet may use parallel computation.}
}
\value{An updated LidarData object.}
\seealso{\code{\link{LiDAR-package}}}
