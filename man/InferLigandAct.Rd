\name{InferLigandAct}
\alias{InferLigandAct}
\title{Infer ligand-associated activity}
\description{Infer ligand-associated activity.}
\usage{
InferLigandAct(lidardata, normalize = TRUE, parallel = FALSE)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{normalize}{Normalize activity to log1p of proportions times 10000.}
\item{parallel}{Use parallel execution in dnet. Default FALSE.}
}
\value{An updated LidarData object.}
\details{The core propagation and score formulas follow the supplied method script. Database loading returns a list and does not modify the global environment.}
