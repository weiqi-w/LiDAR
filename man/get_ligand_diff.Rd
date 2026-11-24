\name{get_ligand_diff}
\alias{get_ligand_diff}
\title{Extract decay fits and optionally filter result ranges.}
\description{Extract decay fits and optionally filter result ranges.}
\usage{
get_ligand_diff(lidardata, ligand_list = NULL, coupling = NULL, pval = NULL, n = NULL, ki = NULL, p_column = "p.value", r50 = NULL, r01 = NULL, r = NULL)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{ligand_list}{Optional ligand identifiers; NULL uses available or filtered ligands as documented.}
\item{r}{Legacy argument name for coupling.}
\item{coupling}{Exposure-response coupling threshold for SelectLigands or two-element inclusive range for GetDecay.}
\item{pval}{Strict P-value upper threshold.}
\item{n}{Positive decay exponent, candidate neighbor counts for BuildNetworks, a filter range for GetDecay, or number of colors for select_distinct_colors.}
\item{ki}{Positive relative decay parameter, or a two-element filter range for GetDecay.}
\item{p_column}{P-value column: p.value or p.empirical.}
\item{r50}{Optional inclusive two-element half-decay distance range.}
\item{r01}{Optional inclusive two-element one-percent distance range.}
}
\value{A data frame.}
\seealso{\code{\link{LiDAR-package}}}
