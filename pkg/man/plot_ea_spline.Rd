\name{plot_ea_spline}
\docType{methods}
\alias{plot_ea_spline}
\title{plot soil profile outputs from \code{ea_spline}}
\description{This is a simple function for plotting outputs and soil profile data from the \code{ea_spline} function.}

\usage{plot_ea_spline(splineOuts, d = t(c(0,5,15,30,60,100,200)), maxd, type = 1 , label = "", plot.which = 1)}
\arguments{
  \item{splineOuts}{list; This is a returned object from the \code{ea_spline} function} 
  \item{d}{numeric; standard depths that were used during the fitting from \code{ea_spline}.} 
  \item{maxd}{numeric; Maximum soil depth for generating the outputs plots/s}
   \item{type}{numeric; Different themes of plot may be exported. Type 1 is to return the observed soil data plus the continuous spline (default). Type 2 is to return the observed data plus the averages of the spline at the specified depth intervals. Type 3 is to return the observed data, spline averages and continuous spline.}
   \item{label}{Character; Optional label to put on x-axis of plot.}
   \item{plot.which}{numeric; Integer selction of which plot to produce if multiple soil profiles have been fitted using the \code{ea_spline} function.}
}
\value{Function returns a labeled plot of the soil profile information and user selected \code{ea_spline} outputs.}

\note{This is a companion function for \code{ea_spline}. This function is really designed to plot one soil profile at a time. Would need to be inserted into another function if plotting of a collection of soil profiles is required.}

\author{Brendan Malone}


\examples{
library(ithir)
library(aqp)

data(oneProfile)
str(oneProfile)

#convert to SoilProfileCollection object
depths(oneProfile)<- Soil.ID ~ Upper.Boundary + Lower.Boundary


#fit spline
eaFit <- ea_spline(oneProfile, var.name="C.kg.m3.",d= t(c(0,5,15,30,60,100,200)),lam = 0.1, vlow=0, show.progress=FALSE )

#do plot
plot_ea_spline(splineOuts=eaFit, d= t(c(0,5,15,30,60,100,200)), maxd=200, type=1, label="carbon density") 
}
\keyword{methods}
