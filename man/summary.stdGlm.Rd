\name{summary.stdGlm}
\alias{summary.stdGlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summarizes GLM regression standardization fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a \code{summary} method for class \code{"stdGlm"}.
}
\usage{
\method{summary}{stdGlm}(object, CI.type = "plain", CI.level = 95,
  transform = NULL, contrast = NULL, reference = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
an object of class \code{"stdGlm"}.
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
a string. If set to \code{"log"}, \code{"logit"}, or \code{"odds"}, the standardized 
mean \eqn{\theta(x)} is transformed into \eqn{\psi(x)=log\{\theta(x)\}}, 
\eqn{\psi(x)=log[\theta(x)/\{1-\theta(x)\}]}, or 
\eqn{\psi(x)=\theta(x)/\{1-\theta(x)\}}, respectively. If left unspecified, 
\eqn{\psi(x)=\theta(x)}.  
}

\item{contrast}{
%%     ~~Describe \code{CI.type} here~~
a string. If set to \code{"difference"} or \code{"ratio"}, then \eqn{\psi(x)-\psi(x_0)}
or \eqn{\psi(x) / \psi(x_0)} are constructed, where \eqn{x_0} is a reference 
level specified by the \code{reference} argument. 
}

\item{reference}{
%%     ~~Describe \code{CI.type} here~~
must be specified if \code{contrast} is specified. 
}

 \item{\dots}{
%%     ~~Describe \code{\dots} here~~
not used.
}

}
\author{
%%  ~~who you are~~
Arvid Sjolander
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{stdGlm}}
}
\examples{

##See documentation for stdGlm

}
