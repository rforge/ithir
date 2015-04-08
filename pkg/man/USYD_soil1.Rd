\name{USYD soil data}
\docType{data}
\alias{USYD soil data}
\title{Random selection of soil point data}
\description{A soil information \code{dataframe} for 29 soil profiles. Each row is a horizon or depth interval observation for a given soil profile. Various soil attribute information is attibuted to each observation }
\usage{data(USYD_soil1)}
\format{
\code{hvPoints250} is 166 row \code{dataframe} with identifier, soil depth and soil attribute information labelled in respective columns. 
}
\details{This data frame is a typical soil information table}

\references{
\itemize{
\item This data is copyright of the Soil Security Lab, The University of Sydney. If any part of this data is to be used in any publication or report, please provide a citation:
Soil Security Laboratory, 2015. Use R for Digital Soil Mapping Manual. The University of Sydney, Sydney, Australia
}
}

\examples{
library(ithir)


data(USYD_soil1)
summary(USYD_soil1)

}
\keyword{datasets}
