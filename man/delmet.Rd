\name{delmet}

\alias{delmet}

\title{
Delta method
}

\description{
\code{delmet} uses the delta method to estimate the variance of a 
scalar function of a parameter vector of length \code{n}. 
}

\usage{
delmet(fun, est, vcov)
}

\arguments{
  \item{fun}{
the scalar function.  
} 
  \item{est}{
either a vector of length \code{n} containing the estimate of the parameter 
vector, or a matrix with \code{p} rows and \code{n} columns, where the \code{p} 
rows are a sequence of estimates of the parameter vector. 
}
  \item{vcov}{
(an estimate of) the variance-covariance matrix of the parameter estimates. If 
\code{est} is a vector of length \code{n}, then \code{vcov} should be a 
\code{nxn} matrix. If \code{est} is a matrix with \code{p} rows and \code{n}
columns, then \code{vcov} should be a list of length \code{p}, where each 
element is a \code{nxn} matrix. 
}
}

\value{
An object of class \code{"delmet"} is a vector containing the variance of
the estimated function.
}

\author{
Arvid Sjolander. 
}

