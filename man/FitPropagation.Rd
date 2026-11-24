\name{FitPropagation}
\alias{FitPropagation}
\title{Fit decay profiles and optionally compute reoptimized permutation P values.}
\description{Fit decay profiles and optionally compute reoptimized permutation P values.}
\usage{
FitPropagation(lidardata, CalcLigandConc = TRUE, imp = TRUE, ligand = NULL, permutation_p = FALSE, n_perm = 999L, perm_seed = 123L, fit_timeout = 300)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{CalcLigandConc}{Whether to retain fitted exposure in the returned object.}
\item{imp}{Whether to use imputed expression.}
\item{ligand}{Ligand identifier; NULL selects available ligands where supported.}
\item{permutation_p}{Whether to run and report the empirical permutation test; defaults to FALSE.}
\item{n_perm}{Number of permutations, default 999; use 9999 when finer resolution is needed.}
\item{perm_seed}{Nonnegative seed for reproducible permutations; caller RNG state is restored.}
\item{fit_timeout}{Cooperative time limit in seconds per fit; defaults to 300.}
}
\value{An updated LidarData object.}
\seealso{\code{\link{LiDAR-package}}}
