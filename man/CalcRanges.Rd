\name{CalcRanges}
\alias{CalcRanges}
\title{Calculate continuous half-decay and one-percent distances.}
\description{Calculate continuous half-decay and one-percent distances.}
\usage{
CalcRanges(ki, n, B)
}
\arguments{
\item{ki}{Positive relative decay parameter, or a two-element filter range for GetDecay.}
\item{n}{Positive decay exponent, candidate neighbor counts for BuildNetworks, a filter range for GetDecay, or number of colors for select_distinct_colors.}
\item{B}{Maximum decay distance in input coordinate units.}
}
\value{A named list of results.}
\seealso{\code{\link{LiDAR-package}}}
