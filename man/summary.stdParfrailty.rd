\name{summary.stdParfrailty}

\alias{summary.stdParfrailty}

\title{
Summarizes Frailty standardization fit
}

\description{
This is a \code{summary} method for class \code{"stdParfrailty"}.
}

\usage{
\method{summary}{stdParfrailty}(object, t, CI.type = "plain", CI.level = 0.95,
  transform = NULL, contrast = NULL, reference = NULL, ...)
}

\arguments{
  \item{object}{
an object of class \code{"stdParfrailty"}.
} 
  \item{t}{
numeric, indicating the times at which to summarize. It defaults to the specified 
value(s) of the argument \code{t} in the \code{stdCox} function.   
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
  \item{\dots}{
not used.
}
}

\author{
Arvid Sjolander
}

\seealso{
\code{\link{stdParfrailty}}
}

\examples{

##See documentation for stdParfrailty

}
