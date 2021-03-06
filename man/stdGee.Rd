\name{stdGee}

\alias{stdGee}

\title{
Regression standardization in conditional generalized estimating equations
}

\description{
\code{stdGee} performs regression standardization in linear and log-linear 
fixed effects models, at specified values of the exposure, over the sample 
covariate distribution. Let \eqn{Y}, \eqn{X}, and \eqn{Z} be the outcome, 
the exposure, and a vector of covariates, respectively. It is assumed that data
are clustered with a cluster indicator \eqn{i}. \code{stdGee} uses 
fitted fixed effects model, with cluster-specific intercept \eqn{a_i} 
(see \code{details}), to estimate the standardized mean 
\eqn{\theta(x)=E\{E(Y|i,X=x,Z)\}}, where \eqn{x} is a specific value of \eqn{X}, 
and the outer expectation is over the marginal distribution of \eqn{(a_i,Z)}. 
}

\usage{
stdGee(fit, data, X, x, clusterid, subsetnew)
}

\arguments{
  \item{fit}{
an object of class \code{"gee"}, with argument \code{cond = TRUE}, as returned 
by the \code{gee} function in the \pkg{drgee} package. If arguments 
\code{weights} and/or \code{subset} are used when fitting the model, then the 
same weights and subset are used in \code{stdGee}.  
} 
  \item{data}{
a data frame containing the variables in the model. This should be the same 
data frame as was used to fit the model in \code{fit}.
}
  \item{X}{
a string containing the name of the exposure variable \eqn{X} in \code{data}.
}
  \item{x}{
an optional vector containing the specific values of \eqn{X} at which to estimate 
the standardized mean. If \eqn{X} is binary (0/1) or
a factor, then \code{x} defaults to all values of \eqn{X}. If \eqn{X} is numeric,
then \code{x} defaults to the mean of \eqn{X}. If \code{x} is set to \code{NA},
then \eqn{X} is not altered. This produces an estimate of the marginal mean
 \eqn{E(Y)=E\{E(Y|X,Z)\}}. 
}
  \item{clusterid}{
an mandatory string containing the name of a cluster identification variable.
Must be identical to the clusterid variable used in the gee call. 
}

  \item{subsetnew}{
an optional logical statement specifying a subset of observations to be used in 
the standardization. This set is assumed to be a subset of the subset (if any)
that was used to fit the regression model. 
}
}

\details{
\code{stdGee} assumes that a fixed effects model 
\deqn{\eta\{E(Y|i,X,Z)\}=a_i+h(X,Z;\beta)} 
has been fitted. The link function \eqn{\eta} is assumed to be the identity link
or the log link. The conditional generalized estimating equation (CGGE)
estimate of \eqn{\beta} is used to obtain estimates of the cluster-specific 
means:
\deqn{\hat{a}_i=\sum_{j=1}^{n_i}r_{ij}/n_i,}
where 
\deqn{r_{ij}=Y_{ij}-h(X_{ij},Z_{ij};\hat{\beta})}
if \eqn{\eta} is the identity link, and
\deqn{r_{ij}=Y_{ij}exp\{-h(X_{ij},Z_{ij};\hat{\beta})\}}
if \eqn{\eta} is the log link, and \eqn{(X_{ij},Z_{ij})} is the value of 
\eqn{(X,Z)} for subject \eqn{j} in cluster \eqn{i}, \eqn{j=1,...,n_i}, 
\eqn{i=1,...,n}. The CGEE estimate of \eqn{\beta} and the estimate of 
\eqn{a_i} are used to estimate the mean \eqn{E(Y|i,X=x,Z)}:
\deqn{\hat{E}(Y|i,X=x,Z)=\eta^{-1}\{\hat{a}_i+h(X=x,Z;\hat{\beta})\}.} 
For each \eqn{x} in the \code{x} argument, these estimates are averaged across 
all subjects (i.e. all observed values of \eqn{Z} and all estimated values of 
\eqn{a_i}) to produce estimates  
\deqn{\hat{\theta}(x)=\sum_{i=1}^n \sum_{j=1}^{n_i} \hat{E}(Y|i,X=x,Z_i)/N,}
where \eqn{N=\sum_{i=1}^n n_i}. The variance for \eqn{\hat{\theta}(x)} is 
obtained by the sandwich formula. 
}

\value{
An object of class \code{"stdGee"} is a list containing 
\item{call}{
  the matched call. 
  }
\item{input}{
  \code{input} is a list containing all input arguments.     
  }
\item{est}{
  a vector with length equal to \code{length(x)}, where element \code{j} is equal to 
  \eqn{\hat{\theta}}(\code{x[j]}).
  }
\item{vcov}{
  a square matrix with \code{length(x)} rows, where the element 
  on row \code{i} and column \code{j} is the (estimated) covariance of 
  \eqn{\hat{\theta}}(\code{x[i]}) and 
  \eqn{\hat{\theta}}(\code{x[j]}).
  }
}

\references{
Goetgeluk S. and Vansteelandt S. (2008). Conditional generalized estimating 
equations for the analysis of clustered and longitudinal data. 
\emph{Biometrics} \bold{64}(3), 772-780.

Martin R.S. (2017). Estimation of average marginal effects in multiplicative 
unobserved effects panel models. \emph{Economics Letters} \bold{160}, 16-19.

Sjolander A. (2019). Estimation of marginal causal effects in the presence of 
confounding by cluster. \emph{Biostatistics} doi: 10.1093/biostatistics/kxz054
}

\author{
Arvid Sjolander. 
}

\note{
The variance calculation performed by \code{stdGee} does not condition on 
the observed covariates \eqn{\bar{Z}=(Z_{11},...,Z_{nn_i})}. To see how this 
matters, note that 
\deqn{var\{\hat{\theta}(x)\}=E[var\{\hat{\theta}(x)|\bar{Z}\}]+var[E\{\hat{\theta}(x)|\bar{Z}\}].} 
The usual parameter \eqn{\beta} in a generalized linear model does not depend 
on \eqn{\bar{Z}}. Thus, \eqn{E(\hat{\beta}|\bar{Z})} is 
independent of \eqn{\bar{Z}} as well (since \eqn{E(\hat{\beta}|\bar{Z})=\beta}), 
so that the term \eqn{var[E\{\hat{\beta}|\bar{Z}\}]} in the corresponding 
variance decomposition for \eqn{var(\hat{\beta})} becomes equal to 0. However, 
\eqn{\theta(x)} depends on \eqn{\bar{Z}} through the average over the sample 
distribution for \eqn{Z}, and thus the term \eqn{var[E\{\hat{\theta}(x)|\bar{Z}\}]} 
is not 0, unless one conditions on \eqn{\bar{Z}}. 
}

\examples{

require(drgee)

n <- 1000
ni <- 2
id <- rep(1:n, each=ni)
ai <- rep(rnorm(n), each=ni)
Z <- rnorm(n*ni)
X <- rnorm(n*ni, mean=ai+Z)
Y <- rnorm(n*ni, mean=ai+X+Z+0.1*X^2)
dd <- data.frame(id, Z, X, Y)
fit <- gee(formula=Y~X+Z+I(X^2), data=dd, clusterid="id", link="identity",
  cond=TRUE)
fit.std <- stdGee(fit=fit, data=dd, X="X", x=seq(-3,3,0.5), clusterid="id")
print(summary(fit.std, contrast="difference", reference=2))
plot(fit.std)

}
