\name{LiDAR-classes}
\alias{LiDARMatrix-class}
\alias{Data-class}
\alias{ImputationData-class}
\alias{Image-class}
\alias{Coembedding-class}
\alias{Diffusion-class}
\alias{Lidarscore-class}
\alias{DB-class}
\alias{LidarData-class}
\title{LiDAR result containers}
\description{S4 containers for expression, coordinates, embeddings, networks, fitted propagation and scores. Create objects with CreateLiDAR rather than constructing slots manually.}
\details{LidarData contains Data, ImputationData, Image, Coembedding, Diffusion, Lidarscore and DB components. Coordinates retain their original units. Failure records are stored in the failure_log attribute.}
\seealso{\code{\link{CreateLiDAR}}}
