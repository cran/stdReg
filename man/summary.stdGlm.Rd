\name{summary.stdGlm}

\alias{summary.stdGlm}

\title{
Summarizes GLM regression standardization fit
}

\description{
This is a \code{summary} method for class \code{"stdGlm"}.
}

\usage{
\method{summary}{stdGlm}(object, CI.type = "plain", CI.level = 0.95,
  transform = NULL, contrast = NULL, reference = NULL, ...)
}

\arguments{
  \item{object}{
an object of class \code{"stdGlm"}.
} 
 \item{CI.type}{
string, indicating the type of confidence intervals. Either "plain", which
gives untransformed intervals, or "log", which gives log-transformed intervals.
}
  \item{CI.level}{
desired coverage probability of confidence intervals, on decimal form. 
}
  \item{transform}{
a string. If set to \code{"log"}, \code{"logit"}, or \code{"odds"}, the standardized 
mean \eqn{\theta(x)} is transformed into \eqn{\psi(x)=log\{\theta(x)\}}, 
\eqn{\psi(x)=log[\theta(x)/\{1-\theta(x)\}]}, or 
\eqn{\psi(x)=\theta(x)/\{1-\theta(x)\}}, respectively. If left unspecified, 
\eqn{\psi(x)=\theta(x)}.  
}
  \item{contrast}{
a string. If set to \code{"difference"} or \code{"ratio"}, then \eqn{\psi(x)-\psi(x_0)}
or \eqn{\psi(x) / \psi(x_0)} are constructed, where \eqn{x_0} is a reference 
level specified by the \code{reference} argument. 
}
  \item{reference}{
must be specified if \code{contrast} is specified. 
}
 \item{\dots}{
not used.
}
}

\author{
Arvid Sjolander
}

\seealso{
\code{\link{stdGlm}}
}

\examples{

##See documentation for stdGlm

}
