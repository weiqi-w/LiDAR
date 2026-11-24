\name{BuildNetwork}
\alias{BuildNetwork}
\title{Build one simple undirected cell-gene neighbor graph.}
\description{Build one simple undirected cell-gene neighbor graph.}
\usage{
BuildNetwork(lidardata, nn = 30)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{nn}{Number of neighbors, excluding self.}
}
\value{An igraph graph.}
\seealso{\code{\link{LiDAR-package}}}
