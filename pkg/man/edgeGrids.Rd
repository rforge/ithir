\name{edgeroi covariates (subset)}
\docType{data}
\alias{edgeroi covariates}
\alias{edgeGrids}
\title{Selected subset of environmental covariates for the Edgeroi District, NSW}
\description{A \code{RasterStack} of selected environmental covariates for a small area of the Edgeroi district in NSW, Australia.}
\usage{data(edgeGrids)}
\format{
\code{edgeGrids} is an 87 row, 117 column, 5 layer \code{RasterStack} of selected environmental covariates from a small area in the Edgeroi district NSW, Australia. The grids have a pixel resolution of 90m x 90m. It contains the following layers:
  \describe{
  \item{\code{Doserate}}{numeric; gamma radiometric data}
  \item{\code{Elevation}}{numeric; topographic variable of bare earth ground elevation. Derived from digital elevation model}
  \item{\code{Panchromat}}{numeric; panchromatic band of the Landsat 7 satelite}
  \item{\code{Slope}}{numeric; Slope gradient of the land surface. Derived from digital elevation model }  
  \item{\code{TWI}}{numeric; topographic wetness index. Secondary derivative of the digital elevation model }  
}
}
\details{The Edgeroi District, NSW is an intensive cropping area upon the fertile alluvial Namoi River plain. The District has been the subject of many soil invetigations, namely McGarry et al. (1989) whom describe an extensive soil data set collected from the area. More recently, digital soil mapping studies of the area have been conducted, for example, Malone et al. (2009).}
\note{The raw spatial data that contributed to the creation of \code{edgeGrids} were sourced from publically accessable repositories hosted by various Australian Government and international agencies including CSIRO (for the DEM), Geosciences Australia (for the radiometric data) and NASA (for the Landsat 7 ETM+ data). The projection for the \code{RasterStack} is WGS 84 Zone 55.}

\references{
\itemize{
\item Malone, B.P., McBratney, A.B., Minasny, B. (2009) \href{http://dx.doi.org/10.1016/j.geoderma.2009.10.007}{Mapping continuous depth functions of soil carbon storage and available water capacity}. Geoderma 154, 138-152.
\item McGarry, D., Ward, W.T., McBratney, A.B. (1989) Soil Studies in the Lower Namoi Valley: Methods and Data. The Edgeroi Data Set. (2 vols) (CSIRO Division of Soils: Adelaide).
}
}
\examples{
library(ithir)
library(raster)

data(edgeGrids)
edgeGrids

#plot 1 of the \code{RasterStack} layers
plot(edgeGrids[[1]], main= names(edgeGrids[[1]]))

}
\keyword{datasets}
