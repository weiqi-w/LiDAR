\name{PlotDecay}
\alias{PlotDecay}
\title{Plot a continuous fitted decay profile and optional characteristic distances.}
\description{Plot a continuous fitted decay profile and optional characteristic distances.}
\usage{
PlotDecay(lidardata, ligand, show_ranges = TRUE, distance_unit = NULL, n_points = 501L)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{ligand}{Ligand identifier; NULL selects available ligands where supported.}
\item{show_ranges}{Whether to mark half-decay and one-percent distances.}
\item{distance_unit}{Optional display label for unchanged coordinate units.}
\item{n_points}{Number of positions used to draw the continuous decay curve.}
}
\value{A ggplot object.}
\seealso{\code{\link{LiDAR-package}}}
