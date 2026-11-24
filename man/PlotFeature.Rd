\name{PlotFeature}
\alias{PlotFeature}
\title{Plot a spatial feature with optional z-score scaling.}
\description{Plot a spatial feature with optional z-score scaling.}
\usage{
PlotFeature(values, cell_positions, show_border = FALSE, scale = c("zscore", "none"), title = NULL, labs_params = list(color = "Feature"), point_size = 0.05, col_palette = NULL, reverse_y = FALSE, color_limits = NULL, feature = NULL, zscore = TRUE, ratio = 0.8)
}
\arguments{
\item{feature}{Optional spatial feature name such as Exposure or Score for automatic palette selection.}
\item{values}{Named numeric spatial feature vector.}
\item{cell_positions}{Coordinate table with cell, x and y columns.}
\item{show_border}{Plot appearance setting for show border ; see the default value in Usage.}
\item{ratio}{Positive ratio of the displayed length of one y unit to one x unit in spatial panels. The default is 0.8.}
\item{zscore}{Standardize each spatial feature across the plotted units by default. An explicit scale argument takes precedence.}
\item{scale}{none preserves raw values; zscore standardizes each feature.}
\item{title}{Optional plot title.}
\item{labs_params}{Named list of arguments passed to ggplot2::labs.}
\item{point_size}{Point size.}
\item{col_palette}{Optional viridis palette override. NULL uses the feature-specific palette.}
\item{reverse_y}{Whether to reverse the vertical coordinate axis; default FALSE.}
\item{color_limits}{Optional two-element color range.}
}
\value{A ggplot object.}
\seealso{\code{\link{LiDAR-package}}}
