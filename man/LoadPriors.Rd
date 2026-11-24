\name{LoadPriors}
\alias{LoadPriors}
\title{Read user-provided ligand-target and ligand-receptor RDS files.}
\description{Read user-provided ligand-target and ligand-receptor RDS files.}
\usage{
LoadPriors(ligand_target_path, lr_network_path)
}
\arguments{
\item{ligand_target_path}{Path to a target-by-ligand matrix saved with saveRDS.}
\item{lr_network_path}{Path to a receptor-network data frame saved with saveRDS.}
}
\value{A named list of results.}
\seealso{\code{\link{LiDAR-package}}}
