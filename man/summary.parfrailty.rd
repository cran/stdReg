\name{summary.parfrailty}
\alias{summary.parfrailty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Summarizes parfrailty fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a \code{summary} method for class \code{"parfrailty"}.
}
\usage{
\method{summary}{parfrailty}(object, CI.type = "plain", CI.level = 0.95, 
  digits=max(3L, getOption("digits") - 3L), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
an object of class \code{"parfrailty"}.
}
  
  \item{CI.type}{
%%     ~~Describe \code{CI.type} here~~
string, indicating the type of confidence intervals. Either "plain", which
gives untransformed intervals, or "log", which gives log-transformed intervals.
}
  \item{CI.level}{
%%     ~~Describe \code{CI.type} here~~
desired coverage probability of confidence intervals, in decimal form. 
}
\item{digits}{
%%     ~~Describe \code{CI.type} here~~
the number of significant digits to use when printing..  
}
\item{\dots}{
%%     ~~Describe \code{\dots} here~~
not used.
}

}
\author{
%%  ~~who you are~~
Arvid Sjolander and Elisabeth Dahlqwist.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{parfrailty}}
}
\examples{

##See documentation for frailty

}
