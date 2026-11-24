\name{ExportPlots}
\alias{ExportPlots}
\title{Export ligand-specific panels and return a file manifest and diagnostic log.}
\description{Export ligand-specific panels and return a file manifest and diagnostic log.}
\usage{
ExportPlots(lidardata, cells = NULL, ligand_list = NULL, receptor_list = NULL, features = c("LigandExpr", "LigandExprImp", "Response", "Exposure", "Availability",      "AvailabilityImp", "Decay", "Coupling"), out_dir = NULL, scale = c("zscore", "none"), show_border = FALSE, col_palette = NULL, reverse_y = FALSE, shared_scales = TRUE, distance_unit = NULL, format = c("jpg", "png", "pdf"), dpi = 300, point_size = 0.05, overwrite = FALSE, zscore = TRUE, ratio = 0.8)
}
\arguments{
\item{lidardata}{A LidarData object.}
\item{cells}{Optional spatial-unit identifiers in requested output order.}
\item{ligand_list}{Optional ligand identifiers; NULL uses available or filtered ligands as documented.}
\item{receptor_list}{Optional receptor identifiers.}
\item{features}{Feature names to extract or plot.}
\item{out_dir}{Output directory; NULL creates a timestamped directory.}
\item{ratio}{Positive ratio of the displayed length of one y unit to one x unit in spatial panels. The default is 0.8.}
\item{zscore}{Standardize each spatial feature across the plotted units by default. An explicit scale argument takes precedence.}
\item{scale}{none preserves raw values; zscore standardizes each feature.}
\item{show_border}{Plot appearance setting for show border ; see the default value in Usage.}
\item{col_palette}{Optional viridis palette override. NULL uses the feature-specific palette.}
\item{reverse_y}{Whether to reverse the vertical coordinate axis; default FALSE.}
\item{shared_scales}{Share limits across ligand panels for the same feature; different features remain separate.}
\item{distance_unit}{Optional display label for unchanged coordinate units.}
\item{format}{Image format: jpg, png or pdf.}
\item{dpi}{Raster export resolution.}
\item{point_size}{Point size.}
\item{overwrite}{Whether existing output files may be replaced.}
}
\value{Invisibly, a list containing files, log, color_limits, scale, shared_scales and out_dir.}
\seealso{\code{\link{LiDAR-package}}}
