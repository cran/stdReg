\name{confint.stdGee}

\alias{confint.stdGee}

\title{
Confidence interval
}

\description{
This is a \code{confint} method for class \code{"stdGee"}.  
}

\usage{
\method{confint}{stdGee}(object, parm, level = 0.95, fun, type="plain", ...)
}

\arguments{
  \item{object}{
an object of class \code{"stdGee"}.  
} 
  \item{parm}{
not used.
}
  \item{level}{
the coverage probability of the confidence intervals. 
}
  \item{fun}{
a function of one vector argument of length \code{p}, which returns a scalar. 
}
  \item{type}{
a string specifying the type of confidence interval; \code{plain} 
(for untransformed) or \code{log} (for log-transformed). 
}
  \item{\dots}{
not used.
}

}

\details{
\code{confint.stdGee} extracts the \code{est} element from \code{object}, and 
inputs this to \code{fun}. It then uses the delta method to compute a confidence 
interval for the output of \code{fun}.   
}

\value{
a matrix with 1 row and 2 columns, containing the computed confidence interval. 
}

\author{
Arvid Sjolander. 
}

