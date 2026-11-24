\name{select_distinct_colors}
\alias{select_distinct_colors}
\title{Choose separated colors from the LiDAR palette.}
\description{Choose separated colors from the LiDAR palette.}
\usage{
select_distinct_colors(n, palette = lidar_colors)
}
\arguments{
\item{n}{Positive decay exponent, candidate neighbor counts for BuildNetworks, a filter range for GetDecay, or number of colors for select_distinct_colors.}
\item{palette}{Character vector of colors.}
}
\value{A vector or matrix of calculated values with identifiers retained where applicable.}
\seealso{\code{\link{LiDAR-package}}}
