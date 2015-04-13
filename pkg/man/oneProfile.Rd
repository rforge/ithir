\name{One soil profile}
\docType{data}
\alias{Soil carbon soil profile}
\title{One soil profile}
\description{A soil information \code{dataframe} of one soil profile for one soil attribute; soil carbon density. Each row is a specified depth interval.}
\usage{data(oneProfile)}
\format{
\code{oneProfile} is an 8 row \code{dataframe} with 4 columns. Columns correspond to identifier, upper soil depth interval, lower depth interval and soil carbon density value. 
}
\details{This data frame is a typical soil data table}

\references{
\itemize{
\item This data is copyright of the Soil Security Lab, The University of Sydney. If any part of this data is to be used in any publication or report, please provide a citation:
Soil Security Laboratory, 2015. Use R for Digital Soil Mapping Manual. The University of Sydney, Sydney, Australia
}
}

\examples{
library(ithir)


data(oneProfile)
str(oneProfile)

}
\keyword{datasets}
