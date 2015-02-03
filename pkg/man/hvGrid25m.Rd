\name{Hunter Valley Grid (S1)}
\docType{data}
\alias{Hunter Valley Grid}
\title{Raster grid of the Lower Hunter Valley, NSW, Australia}
\description{A \code{RasterLayer} of pixel numbers for the area covering the Lower Hunter Valley, NSW, Australia. The pixel resolution is 25m}
\usage{data(hvGrid25m)}
\format{
\code{hvGrid25m} is an 860 row, 676 column, 1 layer \code{Rasterlayer} of pixel numbers for the Lower Hunter Valley, NSW, Australia. The grid has a pixel resolution of 25m x 25m. It contains the following layers:
}
\details{The area in question is the Hunter Wine Country Private Irrigation District (HWCPID), situated in the Lower Hunter Valley, NSW (32.83S 151.35E), and covers an area of approximately 220 km2. The HWCPID is approximately 140 km north of Sydney, NSW, Australia. Climatically, the HWCPID is situated in a temperate climatic zone, and experiences warm humid summers, and relatively cooler yet also humid winters. Rainfall is mostly uniformly distributed throughout the year. On average the HWCPID receives just over 750 mm of rainfall annually. In terms of landuse, an expansive viticultural industry is situated in the area and is possibly most widespread
of rural industries, followed by dryland agricultural grazing systems.The projection for the \code{RasterLayer} is WGS 84 Zone 56. See Malone et al. (2014) for more detailed information about the area and research that has been conducted there by the University of Sydney.}

\references{
\itemize{
\item Malone, B.P., Hughes, P.,  McBratney, A.B., Minasny, B. (2009) \href{http://dx.doi.org/10.1016/j.geodrs.2014.08.001}{A model for the identification of terrons in the Lower Hunter
Valley, Australia}. Geoderma Regional 1, 31-47.
}
}

\examples{
library(ithir)
library(raster)

data(hvGrid25m)
hvGrid25m
plot(hvGrid25m)

}
\keyword{datasets}
