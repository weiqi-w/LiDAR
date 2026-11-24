\name{split_rwr_by_receptor}
\alias{split_rwr_by_receptor}
\title{Align one candidate response with binary receptor availability.}
\description{Align one candidate response with binary receptor availability.}
\usage{
split_rwr_by_receptor(rwrs, ligand, lr_list, receptor_exp, receptorname = NULL, i = 1L)
}
\arguments{
\item{rwrs}{List of candidate node-by-ligand random-walk matrices.}
\item{ligand}{Ligand identifier; NULL selects available ligands where supported.}
\item{lr_list}{Ligand-receptor data frame.}
\item{receptor_exp}{Receptor-by-cell expression matrix.}
\item{receptorname}{Compatibility argument; receptor matches are derived from lr_list.}
\item{i}{Candidate matrix index.}
}
\value{A data frame.}
\seealso{\code{\link{LiDAR-package}}}
