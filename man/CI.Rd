\name{CI}

\alias{CI}

\title{
Confidence interval
}

\description{
\code{CI} computes Wald-type confidence intervals. 
}

\usage{
CI(est, var, CI.type="plain", CI.level=0.95)
}

\arguments{
  \item{est}{
a vector of length \code{p} containing the estimated parameters.  
} 
  \item{var}{
a vector of length \code{p} containing the variances of the estimated parameters. 
}
  \item{CI.type}{
a string specifying the type of confidence interval; \code{plain} or
\code{log} (for log-transformed). 
}
  \item{CI.level}{
the coverage probability of the confidence intervals. 
}

}

\value{
An object of class \code{"CI"} is a matrix with \code{p} rows and 2 columns,
containing the computed confidence intervals. 
}

\author{
Arvid Sjolander. 
}

