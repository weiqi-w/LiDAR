\name{FindNeighbors}
\alias{FindNeighbors}
\title{Find nearest neighbors excluding each observation itself.}
\description{Find nearest neighbors excluding each observation itself.}
\usage{
FindNeighbors(embedding, nn = 30)
}
\arguments{
\item{embedding}{Numeric embedding matrix with named observation rows.}
\item{nn}{Number of neighbors, excluding self.}
}
\value{A named list of results.}
\seealso{\code{\link{LiDAR-package}}}
