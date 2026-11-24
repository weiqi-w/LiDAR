\name{SummarizeDecay}
\alias{SummarizeDecay}
\title{Summarize coupling, continuous ranges and optimization status.}
\description{Summarize coupling, continuous ranges and optimization status.}
\usage{
SummarizeDecay(ligand, response, exposure, fit, B)
}
\arguments{
\item{ligand}{Ligand identifier; NULL selects available ligands where supported.}
\item{response}{Numeric response vector aligned to potential receivers.}
\item{exposure}{Numeric exposure vector aligned to response.}
\item{fit}{Optimization result returned by FitDecay.}
\item{B}{Maximum decay distance in input coordinate units.}
}
\value{A data frame.}
\seealso{\code{\link{LiDAR-package}}}
