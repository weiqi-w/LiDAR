\name{SelectNetwork}
\alias{SelectNetwork}
\title{Select the candidate using correlation with binary receptor availability.}
\description{Select the candidate using correlation with binary receptor availability.}
\usage{
SelectNetwork(rwrs, n_list, ligandname, cell_data, lr_list, receptor_exp, receptorname = NULL, normalize = TRUE)
}
\arguments{
\item{rwrs}{List of candidate node-by-ligand random-walk matrices.}
\item{n_list}{Neighbor counts corresponding to candidate matrices.}
\item{ligandname}{Ligand names to evaluate.}
\item{cell_data}{Spatial-unit names in output order.}
\item{lr_list}{Ligand-receptor data frame.}
\item{receptor_exp}{Receptor-by-cell expression matrix.}
\item{receptorname}{Compatibility argument; receptor matches are derived from lr_list.}
\item{normalize}{Whether to scale selected response vectors by their maximum.}
}
\value{A vector or matrix of calculated values with identifiers retained where applicable.}
\seealso{\code{\link{LiDAR-package}}}
