\name{stdGlm}
\alias{stdGlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Regression standardization in generalized linear models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{stdGlm} performs regression standardization in generalized linear models,
at specified values of the exposure, over the sample covariate distribution.
Let \eqn{Y}, \eqn{X}, and \eqn{Z} be the outcome, the exposure, and a
vector of covariates, respectively. \code{stdGlm} uses a fitted generalized linear
model to estimate the standardized 
mean \eqn{\theta(x)=E\{E(Y|X=x,Z)\}}, where \eqn{x} is a specific value of \eqn{X}, 
and the outer expectation is over the marginal distribution of \eqn{Z}. 
}
\usage{
stdGlm(fit, data, X, x, clusters, case.control = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{formula} here~~
an object of class \code{"glm"}, as returned by the \code{glm} function 
  in the \code{stats} package. 
}
  
  \item{data}{
%%     ~~Describe \code{data} here~~
a data frame containing the variables in the model. This should be the same 
data frame as was used to fit the model in \code{fit}.
}
  \item{X}{
%%     ~~Describe \code{X} here~~
a string containing the name of the exposure variable \eqn{X} in \code{data}.
}
  \item{x}{
%%     ~~Describe \code{x} here~~
an optional vector containing the specific values of \eqn{X} at which to estimate 
the standardized mean. If \eqn{X} is binary (0/1) or
a factor, then \code{x} defaults to all values of \eqn{X}. If \eqn{X} is numeric,
then \code{x} defaults to the mean of \eqn{X}. 
}
 
  \item{clusters}{
%%     ~~Describe \code{clusters} here~~
an optional string containing the name of a cluster identification variable when data are clustered. 
}
  \item{case.control}{
%%     ~~Describe \code{case.control} here~~
logical. Do data come from a case-control study? Defaults to FALSE.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{stdGlm} assums that a generalized linear model 
\deqn{\eta\{E(Y|X,Z)\}=h(X,Z;\beta)}
has been fitted. The maximum likelihood estimate of \eqn{\beta} is used to obtain 
estimates of the mean \eqn{E(Y|X=x,Z)}:
\deqn{\hat{E}(Y|X=x,Z)=\eta^{-1}\{h(X=x,Z;\hat{\beta})\}.} 
For each \eqn{x} in the \code{x} argument, these estimates are averaged across 
all subjects (i.e. all observed values of \eqn{Z}) to produce estimates  
\deqn{\hat{\theta}(x)=\sum_{i=1}^n \hat{E}(t|X=x,Z_i)/n.}
The variance for \eqn{\hat{\theta}(x)} is obtained by the sandwich formula. 

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{"stdGlm"} is a list containing all input arguments to the \code{stdGlm}
function, except \code{data}. In addition the list contains:
\item{fit}{
 the fitted \code{glm} object.
}
\item{means.est}{
  a vector with length equal to \code{length(x)}, where element \code{j} is equal to 
  \eqn{\hat{\theta}}(\code{x[j]}).
  }
\item{means.vcov}{
  a square matrix with \code{length(x)} rows, where the element 
  on row \code{i} and column \code{j} is the (estimated) covariance of 
  \eqn{\hat{\theta}}(\code{x[i]}) and 
  \eqn{\hat{\theta}}(\code{x[j]}).
  }
}
\references{
%% ~put references to the literature/web site here ~
Rothman K.J., Greenland S., Lash T.L. (2008). \emph{Modern Epidemiology}, 3rd edition. 
Lippincott, Williams \& Wilkins.
}
\author{
%%  ~~who you are~~
Arvid Sjolander. 
}
\note{
%%  ~~further notes~~
The variance calculation performed by \code{stdGlm} does not condition on 
the observed covariates \eqn{\bar{Z}=(Z_1,...,Z_n)}. To see how this matters, note that 
\deqn{var\{\hat{\theta}(x)\}=E[var\{\hat{\theta}(x)|\bar{Z}\}]+var[E\{\hat{\theta}(x)|\bar{Z}\}].} 
The usual parameter \eqn{\beta} in a generalized linear model does not depend 
on \eqn{\bar{Z}}. Thus, \eqn{E(\hat{\beta}|\bar{Z})} is 
independent of \eqn{\bar{Z}} as well (since \eqn{E(\hat{\beta}|\bar{Z})=\beta}), so that the 
term \eqn{var[E\{\hat{\beta}|\bar{Z}\}]} in the corresponding variance decomposition 
for \eqn{var(\hat{\beta})} becomes equal to 0. However, \eqn{\theta(x)} depends 
on \eqn{\bar{Z}} through the average over the sample distribution for \eqn{Z}, 
and thus the term \eqn{var[E\{\hat{\theta}(x)|\bar{Z}\}]} is not 0, unless one 
conditions on \eqn{\bar{Z}}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

##Example 1: continuous outcome
n <- 1000
Z <- rnorm(n)
X <- rnorm(n, mean = Z)
Y <- rnorm(n, mean = X + Z + 0.1 * X^2)
dd <- data.frame(Z,X,Y)
fit <- glm(formula = Y ~ X + Z + I(X^2), data = dd)
fit.std <- stdGlm(fit = fit, data = dd, X = "X", x = seq(-3,3,0.5))
print(summary(fit.std))
plot(fit.std)

##Example 2: binary outcome
n <- 1000
Z <- rnorm(n)
X <- rnorm(n, mean = Z)
Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
dd <- data.frame(Z,X,Y)
fit <- glm(formula = Y ~ X + Z + X*Z, family = "binomial", data = dd)
fit.std <- stdGlm(fit = fit, data = dd, X = "X", x = seq(-3,3,0.5))
print(summary(fit.std))
plot(fit.std)
}
