\name{bbRaster}
\docType{methods}
\alias{bbRaster}
\title{Get the bounding box information of a \code{RasterLayer} from its extents.}
\description{This is a simple function that returns a 4 by 2 matrix of a \code{RasterLayer} bounding box}

\usage{bbRaster(obj)
}
\arguments{
  \item{obj}{object of class \code{"RasterLayer"}} 
}
\value{Returns a 4  x 2 matrix with each row indicating a coordinate pair of the bounding box
}

\author{Brendan Malone}


\examples{
library(raster)
data(edgeTarget_C)
target
bbRaster(target)
}
\keyword{methods}
