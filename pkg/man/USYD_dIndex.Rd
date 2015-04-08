\name{USYD drainage index data}
\docType{data}
\alias{USYD soil drainage data}
\title{Soil drainage index data}
\description{A \code{dataframe} of 446 rows. Each row corresponds to a location. Actual location not provided. The \code{dataframe} has 2 variables: \code{DI_observed} and \code{DI_predicted}.}
\usage{data(USYD_dIndex)}
\format{
\code{USYD_dIndex} is 446 row, 2 column \code{dataframe}. 
}
\details{This data frame is a typical soil information table where observations are matched prediction predictions from a specified model.}

\references{
\itemize{
\item This data is copyright of the Soil Security Lab, The University of Sydney. If any part of this data is to be used in any publication or report, please provide a citation:
Soil Security Laboratory, 2015. Use R for Digital Soil Mapping Manual. The University of Sydney, Sydney, Australia
}
}

\examples{
library(ithir)


data(USYD_dIndex)
summary(USYD_dIndex)

}
\keyword{datasets}
