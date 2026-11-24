\name{GetExposure}
\alias{GetExposure}
\title{Extract fitted ligand exposures.}
\description{Extract fitted ligand exposures.}
\usage{
GetExposure(lidardata, ligand_list = NULL, cells = NULL)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{ligand_list}{Optional ligand identifiers; NULL uses available or filtered ligands as documented.}
\item{cells}{Optional spatial-unit identifiers in requested output order.}
}
\value{A vector or matrix of calculated values with identifiers retained where applicable.}
\seealso{\code{\link{LiDAR-package}}}
