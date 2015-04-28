\name{Hunter Valley Points (S1)}
\docType{data}
\alias{Hunter Valley Points}
\title{Random selection of point locations: Hunter Valley}
\description{A \code{dataframe} carrying point coordinates for 250 randomly selected locations in the Lower Hunter Valley, NSW, Australia.}
\usage{data(hvPoints250)}
\format{
\code{hvPoints250} is 250 row \code{dataframe} with the first 2 columns indicating the projected coordinate locations. The coordinate reference system is WGS 84 Zone 56. 
}
\details{The area in question is the Hunter Wine Country Private Irrigation District (HWCPID), situated in the Lower Hunter Valley, NSW (32.83S 151.35E), and covers an area of approximately 220 km2. The HWCPID is approximately 140 km north of Sydney, NSW, Australia. Climatically, the HWCPID is situated in a temperate climatic zone, and experiences warm humid summers, and relatively cooler yet also humid winters. Rainfall is mostly uniformly distributed throughout the year. On average the HWCPID receives just over 750 mm of rainfall annually. In terms of landuse, an expansive viticultural industry is situated in the area and is possibly most widespread
of rural industries, followed by dryland agricultural grazing systems.}

\references{
\itemize{
\item Malone, B.P., Hughes, P.,  McBratney, A.B., Minasny, B. (2009) \href{http://dx.doi.org/10.1016/j.geodrs.2014.08.001}{A model for the identification of terrons in the Lower Hunter
Valley, Australia}. Geoderma Regional 1, 31-47.
}
}

\examples{
library(ithir)


data(hvPoints250)
summary(hvPoints250)

}
\keyword{datasets}
