\name{plot_soilProfile}
\docType{methods}
\alias{plot_soilProfile}
\alias{plot soil profile data}
\title{plot soil profile data}
\description{This is a simple function for plotting soil profile data that is supplied in \code{dataframe} form.}

\usage{plot_soilProfile(data, vals, depths, label="")}
\arguments{
  \item{data}{dataframe; Usually a typical soil information data frame with identifier and soil depth labelled columns, followed by one or more columns of soil attribute or even soil categorical data.} 
  \item{vals}{vector; Would generally point to a column of \code{data}, and is the soil attribute value that needs to be ploted.}     
  \item{depths}{dataframe; Would generally point to the two columns of \code{data} which specify the upper and lower depths of the soil profile observations.}
   \item{label}{character; An optional input the put a label on the x-axis, otherwise it will remain blank.}
}
\value{Function returns a labeled plot of the soil profile information. The length of the bars corresponds to the value of the soil variable being plotted.}

\note{This function is really designed to plot one soil profile at a time. Would need to be inserted into another function if plotting of a collection of soil profiles is required}

\author{Brendan Malone}


\examples{
library(ithir)

data(oneProfile)
str(oneProfile)

#do plot
plot_soilProfile(data = oneProfile, vals = oneProfile$C.kg.m3., depths = oneProfile[,2:3], label= names(oneProfile)[4])
}
\keyword{methods}
