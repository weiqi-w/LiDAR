\name{ComputeScores}
\alias{ComputeScores}
\title{Compute aggregate or sender-receiver pairwise communication scores.}
\description{Compute aggregate or sender-receiver pairwise communication scores.}
\usage{
ComputeScores(lidardata, mode = c("aggregate", "pairwise"), ligand_list = NULL, pval = 0.001, imp = TRUE, coupling = 0, Log = FALSE)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{mode}{aggregate or pairwise scoring.}
\item{ligand_list}{Optional subset of fitted ligands. Active thresholds apply within this subset.}
\item{pval}{Strict upper threshold for the ordinary Pearson P value. NULL disables this filter.}
\item{coupling}{Strict lower threshold for exposure-response coupling. NULL disables this filter.}
\item{Log}{Store log10(1 + score) for each output. Sender and receiver totals are summed before transformation.}
\item{imp}{Whether to use imputed expression.}
}
\value{An updated LidarData object.}
\seealso{\code{\link{LiDAR-package}}}
