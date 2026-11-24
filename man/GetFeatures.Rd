\name{GetFeatures}
\alias{GetFeatures}
\title{Assemble matrices for requested communication features.}
\description{Assemble matrices for requested communication features.}
\usage{
GetFeatures(lidardata, ligand_list = NULL, receptor_list = NULL, cells = NULL, features = c("LigandExpr", "LigandExprImp", "Response", "Exposure", "Availability",      "AvailabilityImp"))
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{ligand_list}{Optional ligand identifiers; NULL uses available or filtered ligands as documented.}
\item{receptor_list}{Optional receptor identifiers.}
\item{cells}{Optional spatial-unit identifiers in requested output order.}
\item{features}{Feature names to extract or plot.}
}
\value{A named list of results.}
\seealso{\code{\link{LiDAR-package}}}
