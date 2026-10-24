\name{load_nichenet_data}
\alias{load_nichenet_data}
\title{Load local NicheNet resources}
\description{Load local NicheNet resources.}
\usage{
load_nichenet_data(species = c("human", "mouse"), data_dir = NULL,
                              ligand_target_path = NULL, lr_network_path = NULL)
}
\arguments{
\item{species}{Species label; database species must be supplied correctly by the user.}
\item{data_dir}{Directory containing the two database RDS files.}
\item{ligand_target_path}{Explicit ligand-target RDS path.}
\item{lr_network_path}{Explicit ligand-receptor RDS path.}
}
\value{A named list of database objects and the species label.}
\details{The core propagation and score formulas follow the supplied method script. Database loading returns a list and does not modify the global environment.}
