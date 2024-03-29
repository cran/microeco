% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{clone}
\alias{clone}
\title{Copy an R6 class object}
\usage{
clone(x, deep = TRUE)
}
\arguments{
\item{x}{R6 class object}

\item{deep}{default TRUE; TRUE means deep copy, i.e. copied object is unlinked with the original one.}
}
\value{
identical but unlinked R6 object
}
\description{
Copy an R6 class object
}
\examples{
data("dataset")
clone(dataset)
}
