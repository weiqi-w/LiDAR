\name{LidarData-class}
\alias{LidarData-class}
\alias{Data-class}
\alias{ImputationData-class}
\alias{Image-class}
\alias{Coembedding-class}
\alias{Diffusion-class}
\alias{Lidarscore-class}
\alias{DB-class}
\alias{LidarMatrix-class}
\title{LiDAR S4 data classes}
\description{Container classes for expression, spatial coordinates, embeddings, diffusion results and communication scores. Matrix slots accept base matrices and Matrix objects.}
\details{Create objects using PrepLidarData. Cells are columns of expression matrices and rows of activity and exposure matrices. Sender and receiver score matrices have ligands in rows. Generic layer class names are retained for compatibility; load the package in a clean session without sourcing the old class definitions.}
