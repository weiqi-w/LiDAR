\name{OptLigandDiff}
\alias{OptLigandDiff}
\title{Fit ligand diffusion and spatial coupling}
\description{Fit ligand diffusion and spatial coupling.}
\usage{
OptLigandDiff(lidardata, CalcLigandConc = TRUE, imp = TRUE, ligand = NULL)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{CalcLigandConc}{Calculate ligand exposure as well as fitted parameters.}
\item{imp}{Use imputed expression.}
\item{ligand}{Optional ligands to fit.}
}
\value{An updated LidarData object.}
\details{The core propagation and score formulas follow the supplied method script. Database loading returns a list and does not modify the global environment.}
