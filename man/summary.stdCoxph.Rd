\name{summary.stdCoxph}
\alias{summary.stdCoxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summarizes Cox regression standardization fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a \code{summary} method for class \code{"stdCoxph"}.
}
\usage{
\method{summary}{stdCoxph}(object, t, CI.type = "plain", CI.level = 0.95,
  transform = NULL, contrast = NULL, reference = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
an object of class \code{"stdCoxph"}.
}
  
  \item{t}{
%%     ~~Describe \code{t} here~~
numeric, indicating the times at which to summarize. It defaults to the specified 
value(s) of the argument \code{t} in the \code{stdCox} function.   
}
  \item{CI.type}{
%%     ~~Describe \code{CI.type} here~~
string, indicating the type of confidence intervals. Either "plain", which
gives untransformed intervals, or "log", which gives log-transformed intervals.
}
  \item{CI.level}{
%%     ~~Describe \code{CI.type} here~~
desired coverage probability of confidence intervals, on decimal form. 
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
not used.
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
