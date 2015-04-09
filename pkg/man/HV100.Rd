\name{Hunter Valley soil data}
\docType{data}
\alias{Hunter Valley soil data}
\title{Soil point data from the Hunter Valley, NSW, Australia}
\description{A soil information \code{dataframe} of 100 obsevations from various locations of the Hunter Valley, NSW, Australia (32.83S 151.35E). The data were collected, using a stratified random sampling design in 2010 and is described in Malone et al. (2011). Each row represents an observation at the 0-5cm depth interval. Various soil attribute information is attibuted to each observation, which includes soil organic carbon, soil pH (1:5 soil:water), and electrical conductivity. Locations are attributed to a site name and have a recorded coordinate location. The CRS of the points is WGS84 UTM Zone 56.}
\usage{data(HV100)}
\format{
\code{HV100} is a 100 row \code{dataframe} with identifier, locational and soil attribute information labelled in respective columns. 
}
\details{This data frame is a typical soil information table}

\references{
\itemize{
\item Malone, B.P., de Gruijter, J.J., McBratney, A.B., Minasny, B., Brus, D.J. (2011) \href{http://dx.doi.org/10.2136/sssaj2010.0280}{Using Additional Criteria for Measuring the Quality of Predictions and Their Uncertainties in a Digital Soil Mapping Framework}. Soil Science Society of America Journal, 75(3): 1032-1043. 
}
}

\examples{
library(ithir)


data(HV100)
head(HV100)

}
\keyword{datasets}
