\name{GetReceptorExpr}
\alias{GetReceptorExpr}
\title{Extract receptor expression.}
\description{Extract receptor expression.}
\usage{
GetReceptorExpr(lidardata, receptor_list = NULL, cells = NULL, imp = TRUE)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{receptor_list}{Optional receptor identifiers.}
\item{cells}{Optional spatial-unit identifiers in requested output order.}
\item{imp}{Whether to use imputed expression.}
}
\value{A vector or matrix of calculated values with identifiers retained where applicable.}
\seealso{\code{\link{LiDAR-package}}}
