\name{CreateNetworks}
\alias{CreateNetworks}
\title{Construct joint cell-gene graphs}
\description{Construct joint cell-gene graphs.}
\usage{
CreateNetworks(lidardata, n)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{n}{Vector of cell-gene graph neighborhood sizes.}
}
\value{An updated LidarData object.}
\details{The core propagation and score formulas follow the supplied method script. Database loading returns a list and does not modify the global environment.}
