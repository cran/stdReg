\name{plot.stdCoxph}

\alias{plot.stdCoxph}

\title{
Plots Cox regression standardization fit
}

\description{
This is a \code{plot} method for class \code{"stdCoxph"}. 
}

\usage{
\method{plot}{stdCoxph}(x, plot.CI = TRUE, CI.type = "plain", CI.level = 0.95,
  transform = NULL, contrast = NULL, reference = NULL, legendpos="bottomleft", ...)
}

\arguments{
  \item{x}{
an object of class \code{"stdCoxph"}.
}
  \item{plot.CI}{
logical, indicating whether confidence intervals should be added to the plot. 
}
  \item{CI.type}{
string, indicating the type of confidence intervals. Either "plain", which
gives untransformed intervals, or "log", which gives log-transformed intervals. 
}
  \item{CI.level}{
desired coverage probability of confidence intervals, on decimal form. 
}
  \item{transform}{
a string. If set to \code{"log"}, \code{"logit"}, or \code{"odds"}, the standardized survival 
function \eqn{\theta(t,x)} is transformed into \eqn{\psi(t,x)=log\{\theta(t,x)\}}, 
\eqn{\psi(t,x)=log[\theta(t,x)/\{1-\theta(t,x)\}]}, or 
\eqn{\psi(t,x)=\theta(t,x)/\{1-\theta(t,x)\}}, respectively. If left unspecified, 
\eqn{\psi(t,x)=\theta(t,x)}.  
}
  \item{contrast}{
a string. If set to \code{"difference"} or \code{"ratio"}, then \eqn{\psi(t,x)-\psi(t,x_0)}
or \eqn{\psi(t,x) / \psi(t,x_0)} are constructed, where \eqn{x_0} is a reference 
level specified by the \code{reference} argument. 
}
  \item{reference}{
must be specified if \code{contrast} is specified. 
}
  \item{legendpos}{
position of the legend; see help for \code{legend}. 
}
  \item{\dots}{
further arguments passed on to plot.default. 
} 
}

\author{
Arvid Sjolander
}

\seealso{
\code{\link{stdCoxph}}
}

\examples{

##See documentation for stdCoxph

}

