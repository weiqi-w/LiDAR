\name{NetworkFromMCA}
\alias{NetworkFromMCA}
\title{Build one simple undirected cell-gene neighbor graph.}
\description{Build one simple undirected cell-gene neighbor graph.}
\usage{
NetworkFromMCA(lidardata, nn = 30)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{nn}{Number of neighbors, excluding self.}
}
\value{An igraph graph.}
\seealso{\code{\link{LiDAR-package}}}
