\name{PrepLidarData}
\alias{PrepLidarData}
\title{Create a LiDAR object and compute MCA embeddings on the selected assay.}
\description{Create a LiDAR object and compute MCA embeddings on the selected assay.}
\usage{
PrepLidarData(data, assay = "SCT", ligand_target_matrix, lr_network, coord = NULL, coord_cols = NULL, cell_col = NULL, layer = "data")
}
\arguments{
\item{data}{A Seurat object containing the selected expression assay.}
\item{assay}{Name of the assay used for expression and MCA.}
\item{ligand_target_matrix}{Numeric target-by-ligand matrix with unique row and column names and nonnegative weights.}
\item{lr_network}{Data frame with from (ligand), to (receptor) and database columns.}
\item{coord}{Optional coordinate data frame indexed by spatial-unit names.}
\item{coord_cols}{Explicit x and y column names, in that order.}
\item{cell_col}{Optional column containing spatial-unit identifiers.}
\item{layer}{Expression layer or legacy assay slot, usually data.}
}
\value{An updated LidarData object.}
\seealso{\code{\link{LiDAR-package}}}
