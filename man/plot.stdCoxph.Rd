\name{plot.stdCoxph}
\alias{plot.stdCoxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plots Cox regression standardization fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a \code{plot} method for class \code{"stdCoxph"}. 
}

\usage{
%%plot.stdCoxph(x, ..., plot.contrasts = FALSE, plot.CI = FALSE, CI.type = "log")
\method{plot}{stdCoxph}(x, plot.CI = TRUE, CI.type = "plain", CI.level = 95,
  transform = NULL, contrast = NULL, reference = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
an object of class \code{"stdCoxph"}.
}
 
  \item{plot.CI}{
%%     ~~Describe \code{plot.CI} here~~
logical, indicating whether confidence intervals should be added to the plot. 
}
  \item{CI.type}{
%%     ~~Describe \code{CI.type} here~~
string, indicating the type of confidence intervals. Either "plain", which
gives untransformed intervals, or "log", which gives log-transformed intervals. 
}
  \item{CI.level}{
%%     ~~Describe \code{CI.type} here~~
desired coverage probability of confidence intervals, in percent. 
}

\item{transform}{
%%     ~~Describe \code{CI.type} here~~
a string. If set to \code{"log"}, \code{"logit"}, or \code{"odds"}, the standardized survival 
function \eqn{\theta(t,x)} is transformed into \eqn{\psi(t,x)=log\{\theta(t,x)\}}, 
\eqn{\psi(t,x)=log[\theta(t,x)/\{1-\theta(t,x)\}]}, or 
\eqn{\psi(t,x)=\theta(t,x)/\{1-\theta(t,x)\}}, respectively. If left unspecified, 
\eqn{\psi(t,x)=\theta(t,x)}.  
}

\item{contrast}{
%%     ~~Describe \code{CI.type} here~~
a string. If set to \code{"difference"} or \code{"ratio"}, then \eqn{\psi(t,x)-\psi(t,x_0)}
or \eqn{\psi(t,x) / \psi(t,x_0)} are constructed, where \eqn{x_0} is a reference 
level specified by the \code{reference} argument. 
}

\item{reference}{
%%     ~~Describe \code{CI.type} here~~
must be specified if \code{contrast} is specified. 
}

  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
further arguments passed on to plot.default. 
} 

}
\author{
%%  ~~who you are~~
Arvid Sjolander
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{stdCoxph}}
}
\examples{

##See documentation for stdCoxph

}

