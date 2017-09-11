\name{summary.parfrailty}

\alias{summary.parfrailty}


\title{
Summarizes parfrailty fit
}

\description{
This is a \code{summary} method for class \code{"parfrailty"}.
}

\usage{
\method{summary}{parfrailty}(object, CI.type = "plain", CI.level = 0.95, 
  digits=max(3L, getOption("digits") - 3L), ...)
}

\arguments{
  \item{object}{
an object of class \code{"parfrailty"}.
}
  \item{CI.type}{
string, indicating the type of confidence intervals. Either "plain", which
gives untransformed intervals, or "log", which gives log-transformed intervals.
}
  \item{CI.level}{
desired coverage probability of confidence intervals, in decimal form. 
}
  \item{digits}{
the number of significant digits to use when printing..  
}
  \item{\dots}{
not used.
}
}

\author{
Arvid Sjolander and Elisabeth Dahlqwist.
}

\seealso{
\code{\link{parfrailty}}
}

\examples{

##See documentation for frailty

}
