\name{GetLigandExpr}
\alias{GetLigandExpr}
\title{Extract ligand expression.}
\description{Extract ligand expression.}
\usage{
GetLigandExpr(lidardata, ligand_list = NULL, cells = NULL, imp = TRUE)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{ligand_list}{Optional ligand identifiers; NULL uses available or filtered ligands as documented.}
\item{cells}{Optional spatial-unit identifiers in requested output order.}
\item{imp}{Whether to use imputed expression.}
}
\value{A vector or matrix of calculated values with identifiers retained where applicable.}
\seealso{\code{\link{LiDAR-package}}}
