\name{CompLidarScore}
\alias{CompLidarScore}
\title{Compute LiDAR communication scores}
\description{Compute LiDAR communication scores.}
\usage{
CompLidarScore(lidardata, mode = c("aggregate","pairwise"),
                          ligand_list = NULL,pval = 0.01, imp = T, include_self = T)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{mode}{Either aggregate or pairwise. Pairwise also calculates sender and receiver scores.}
\item{ligand_list}{Optional fitted ligands to score.}
\item{pval}{Unadjusted correlation p-value cutoff for score selection.}
\item{imp}{Use imputed expression.}
\item{include_self}{Include the first distance bin.}
}
\value{An updated LidarData object.}
\details{The core propagation and score formulas follow the supplied method script. Database loading returns a list and does not modify the global environment.}
