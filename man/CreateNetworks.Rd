\name{CreateNetworks}
\alias{CreateNetworks}
\title{Build candidate graphs for several neighbor counts.}
\description{Build candidate graphs for several neighbor counts.}
\usage{
CreateNetworks(lidardata, n = c(500, 700))
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{n}{Positive decay exponent, candidate neighbor counts for BuildNetworks, a filter range for GetDecay, or number of colors for select_distinct_colors.}
}
\value{An updated LidarData object.}
\seealso{\code{\link{LiDAR-package}}}
