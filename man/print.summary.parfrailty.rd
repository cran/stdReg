\name{print.summary.parfrailty}
\alias{print.summary.parfrailty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Prints summary of parfrailty fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a \code{print} method for class \code{"summary.parfrailty"}.
}
\usage{
%%print.summary.stdFrailty(x, ...)
\method{print}{summary.parfrailty}(x, digits = max(3L, getOption("digits") - 3L),
                             ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
an object of class \code{"summary.parfrailty"}.
}
\item{digits}{
%%     ~~Describe \code{x} here~~
the number of significant digits to use when printing.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
not used.
}
}
\author{
%%  ~~who you are~~
Arvid Sjolander and Elisabeth Dahlqwist
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{parfrailty}}
}
\examples{

##See documentation for frailty

}
