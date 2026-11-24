\name{PlotCoupling}
\alias{PlotCoupling}
\title{Plot exposure-response correspondence across available receivers.}
\description{Plot exposure-response correspondence across available receivers.}
\usage{
PlotCoupling(lidardata, ligand, cells = NULL, point_size = 0.2)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{ligand}{Ligand identifier; NULL selects available ligands where supported.}
\item{cells}{Optional spatial-unit identifiers in requested output order.}
\item{point_size}{Point size.}
}
\value{A ggplot object.}
\seealso{\code{\link{LiDAR-package}}}
