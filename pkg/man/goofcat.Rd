\name{Goodness of fit measures for categorial variable models}
\docType{methods}
\alias{goofcat}
\title{Goodness of fit measures for categorical variable models}
\description{This function performs a number of model diagnostic goodness of fit measures for categorical variables which include: Overall accuracy, Producer's Accuracy, User's Accuracy, and Kappa Statistic. See Congalton (1991) for details on these statistics. They are commonly used to assess the performance of models where the target variable is a class or category.}

\usage{goofcat(observed = NULL, predicted = NULL, conf.mat, imp=FALSE)}
\arguments{
  \item{observed}{a vector values that could either be integer or character class that are actual observations of some categorical phenomenon.} 
  \item{predicted}{a vector values that could either be integer or character class that are predictions of the phenomenon that was observed.}     
  \item{conf.mat}{matrix; Optional input. A square matrix called a confusion matrix that summarizes observation and subsequent prediction from a particular model.}
  \item{imp}{logical; If a confusion matrix is being entered directly into the function through the \code{conf.mat} parametisation, this is set to \code{TRUE}. Default is \code{FALSE}.}
}

\value{Returns an element labelled \code{list} containing the goodness of fit statistics.}

\note{If for whatever reason the \code{conf.mat} input is not of \code{matrix} class, and is not square in dimension the function will halt and provide an error message.}

\author{Brendan Malone}
\references{
\itemize{
\item Congalton, R. G., (1991) \href{http://dx.doi.org/10.1016/0034-4257(91)90048-B}{A review of assessing the accuracy of classifications of remotely sensed data.}. Remote Sensing of the Environment, 37:35-46.
}}


\examples{
library(ithir)

#Using a pre-constructed confusion matrix
con.mat <- matrix(c(5, 0, 1, 2, 0, 15, 0, 5, 0, 1, 31, 0, 0, 10, 2, 11),
nrow = 4, ncol = 4)
rownames(con.mat) <- c("DE", "VE", "CH", "KU")
colnames(con.mat) <- c("DE", "VE", "CH", "KU")
con.mat

goofcat(conf.mat = con.mat, imp=TRUE)

#Using observations and corresponding predictions
#Using random intgers
set.seed(123)
observed<- sample(1:5,1000,replace=TRUE)
set.seed(321)
predicted<- sample(1:5,1000,replace=TRUE)

goofcat(observed = observed, predicted = predicted)

}
\keyword{methods}
