\name{CalcLigandDist}
\alias{CalcLigandDist}
\title{Summarize ligand expression by distance bin}
\description{Summarize ligand expression by distance bin.}
\usage{
CalcLigandDist(lidardata, dist_threshold = 10, distance_bin = NULL, RealDist = F, imp = T, include_self = TRUE)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{dist_threshold}{Number of nonzero distance-bin steps; internally one self bin is added.}
\item{distance_bin}{Positive bin width, or NULL to use minimum positive rounded distance.}
\item{RealDist}{Legacy argument. TRUE is unsupported.}
\item{imp}{Use imputed expression.}
\item{include_self}{Include the first distance bin.}
}
\value{An updated LidarData object.}
\details{The core propagation and score formulas follow the supplied method script. Database loading returns a list and does not modify the global environment.}
