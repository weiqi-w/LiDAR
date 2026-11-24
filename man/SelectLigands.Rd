\name{SelectLigands}
\alias{SelectLigands}
\title{Select fitted ligands by correlation and P value without expression-percentage filtering.}
\description{Select fitted ligands by correlation and P value without expression-percentage filtering.}
\usage{
SelectLigands(lidardata, coupling = 0, pval = 0.01, ligandlist = NULL, p_column = "p.value", r = NULL)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{r}{Legacy argument name for coupling.}
\item{coupling}{Exposure-response coupling threshold for SelectLigands or two-element inclusive range for GetDecay.}
\item{pval}{Strict P-value upper threshold.}
\item{ligandlist}{Optional ligand identifiers.}
\item{p_column}{P-value column: p.value or p.empirical.}
}
\value{A character vector of ligand identifiers.}
\seealso{\code{\link{LiDAR-package}}}
