\name{PrepLidarData}
\alias{PrepLidarData}
\title{Prepare spatial data and MCA embeddings}
\description{Prepare spatial data and MCA embeddings.}
\usage{
PrepLidarData(data, assay = "SCT", ligand_target_matrix, lr_network,
                          coord = NULL)
}
\arguments{
\item{data}{A Seurat object.}
\item{assay}{Assay used for MCA and expression extraction.}
\item{ligand_target_matrix}{Named target-by-ligand prior matrix.}
\item{lr_network}{Data frame with from, to and database columns.}
\item{coord}{Optional two-column coordinates with cell names as row names.}
}
\value{An updated LidarData object.}
\details{The core propagation and score formulas follow the supplied method script. Database loading returns a list and does not modify the global environment.}
