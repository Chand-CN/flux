% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clus.hop.R
\name{clus.hop}
\alias{clus.hop}
\title{Calculate Hopkins statistic}
\usage{
clus.hop(df, rept = 30)
}
\arguments{
\item{df}{a data frame.}

\item{rept}{number of the repeat times.}
}
\description{
Calculate Hopkins statistic
}
\examples{
# example 1
df <- iris[,-5]; clus.hop(df, 10)

# example 2
data(wine, package = 'rattle'); df <- wine[,-1]; clus.hop(df, 10)

}
\author{
Chand Zhang
}
