\name{stdCoxph}

\alias{stdCoxph}

\title{
Regression standardization in Cox proportional hazards models
}

\description{
\code{stdCoxph} performs regression standardization in Cox proportional hazards models,
at specified values of the exposure, over the sample covariate distribution.
Let \eqn{T}, \eqn{X}, and \eqn{Z} be the survival outcome, the exposure, and a
vector of covariates, respectively. \code{stdCoxph} uses a fitted Cox 
proportional hazards model to estimate the standardized 
survival function \eqn{\theta(t,x)=E\{S(t|X=x,Z)\}}, where \eqn{t} is a specific value of \eqn{T}, 
\eqn{x} is a specific value of \eqn{X}, and the expectation is over the marginal 
distribution of \eqn{Z}.     
}

\usage{
stdCoxph(fit, data, X, x, t, clusterid, subsetnew)
}

\arguments{
  \item{fit}{
an object of class \code{"coxph"}, as returned by the \code{coxph} function 
  in the \pkg{survival} package, but without special terms \code{strata}, \code{cluster} or \code{tt}. 
  Only \code{breslow} method for handling ties is allowed. If arguments 
  \code{weights} and/or \code{subset} are used when fitting the model, 
  then the same weights and subset are used in \code{stdGlm}.  
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
the standardized survival function. If \eqn{X} is binary (0/1) or
a factor, then \code{x} defaults to all values of \eqn{X}. If \eqn{X} is numeric,
then \code{x} defaults to the mean of \eqn{X}. If \code{x} is set to \code{NA},
then \eqn{X} is not altered. This produces an estimate of the marginal survival 
function \eqn{S(t)=E\{S(t|X,Z)\}}.   
}
  \item{t}{
an optional vector containing the specific values of \eqn{T} at which to estimate 
the standardized survival function. It defaults to all the observed event times
in \code{data}. 
}
  \item{clusterid}{
an optional string containing the name of a cluster identification variable when data are clustered.
}
  \item{subsetnew}{
an optional logical statement specifying a subset of observations to be used in 
the standardization. This set is assumed to be a subset of the subset (if any)
that was used to fit the regression model. 
}
}

\details{
\code{stdCoxph} assumes that a Cox proportional hazards model 
\deqn{\lambda(t|X,Z)=\lambda_0(t)exp\{h(X,Z;\beta)\}}
has been fitted. Breslow's 
estimator of the cumulative baseline hazard \eqn{\Lambda_0(t)=\int_0^t\lambda_0(u)du}
is used together with the partial likelihood estimate of \eqn{\beta} to obtain 
estimates of the survival function \eqn{S(t|X=x,Z)}:
\deqn{\hat{S}(t|X=x,Z)=exp[-\hat{\Lambda}_0(t)exp\{h(X=x,Z;\hat{\beta})\}].} 
For each \eqn{t} in the \code{t} argument and for each \eqn{x} in the \code{x} argument, 
these estimates are averaged across all subjects (i.e. all observed values of \eqn{Z})
to produce estimates 
\deqn{\hat{\theta}(t,x)=\sum_{i=1}^n \hat{S}(t|X=x,Z_i)/n,} 
where \eqn{Z_i} is the value of \eqn{Z} for subject \eqn{i}, \eqn{i=1,...,n}. 
The variance for \eqn{\hat{\theta}(t,x)} is obtained by the sandwich formula. 
}

\value{
An object of class \code{"stdCoxph"} is a list containing 
\item{call}{
  the matched call. 
  }
\item{input}{
  \code{input} is a list containing all input arguments.     
  }
\item{est}{
  a matrix with \code{length(t)} rows and \code{length(x)} columns, where the element 
  on row \code{i} and column \code{j} is equal to \eqn{\hat{\theta}}(\code{t[i],x[j]}). 
  }
\item{vcov}{
  a list with \code{length(t)} elements. Each element is a square matrix with 
  \code{length(x)} rows. In the \code{k:th} matrix, the element on row \code{i} 
  and column \code{j} is the (estimated) covariance of \eqn{\hat{\theta}}(\code{t[k]},\code{x[i]}) 
  and \eqn{\hat{\theta}}(\code{t[k]},\code{x[j]}).
  }
}

\references{

Chang I.M., Gelman G., Pagano M. (1982). Corrected group prognostic curves
and summary statistics. \emph{Journal of Chronic Diseases} \bold{35}, 669-674.

Gail M.H. and Byar D.P. (1986). Variance calculations for direct adjusted survival
curves, with applications to testing for no treatment effect. \emph{Biometrical Journal}
\bold{28}(5), 587-599. 

Makuch R.W. (1982). Adjusted survival curve estimation using covariates.
\emph{Journal of Chronic Diseases} \bold{35}, 437-443.

Sjolander A. (2016). Regression standardization with the R-package stdReg.
\emph{European Journal of Epidemiology} \bold{31}(6), 563-574.

Sjolander A. (2016). Estimation of causal effect measures with the R-package stdReg.
\emph{European Journal of Epidemiology} \bold{33}(9), 847-858.

}

\author{
Arvid Sjolander
}

\note{
Standardized survival functions are sometimes referred to as (direct) adjusted
survival functions in the literature.

\code{stdCoxph} does not currently handle time-varying exposures or covariates. 

\code{stdCoxph} internally loops over all values in the \code{t} argument. Therefore,
the function will usually be considerably faster if \code{length(t)} is small.

The variance calculation performed by \code{stdCoxph} does not condition on 
the observed covariates \eqn{\bar{Z}=(Z_1,...,Z_n)}. To see how this matters, 
note that 
\deqn{var\{\hat{\theta}(t,x)\}=E[var\{\hat{\theta}(t,x)|\bar{Z}\}]+var[E\{\hat{\theta}(t,x)|\bar{Z}\}].} 
The usual parameter \eqn{\beta} in a Cox proportional hazards model does not 
depend on \eqn{\bar{Z}}. Thus, \eqn{E(\hat{\beta}|\bar{Z})} is independent of 
\eqn{\bar{Z}} as well (since \eqn{E(\hat{\beta}|\bar{Z})=\beta}), so that the 
term \eqn{var[E\{\hat{\beta}|\bar{Z}\}]} in the corresponding variance 
decomposition for \eqn{var(\hat{\beta})} becomes equal to 0. However, 
\eqn{\theta(t,x)} depends on \eqn{\bar{Z}} through the average over the sample 
distribution for \eqn{Z}, and thus the term \eqn{var[E\{\hat{\theta}(t,x)|\bar{Z}\}]} 
is not 0, unless one conditions on \eqn{\bar{Z}}. The variance calculation by 
Gail and Byar (1986) ignores this term, and thus effectively conditions on 
\eqn{\bar{Z}}.      
}

\examples{

require(survival)

n <- 1000
Z <- rnorm(n)
X <- rnorm(n, mean=Z)
T <- rexp(n, rate=exp(X+Z+X*Z)) #survival time
C <- rexp(n, rate=exp(X+Z+X*Z)) #censoring time
U <- pmin(T, C) #time at risk
D <- as.numeric(T < C) #event indicator
dd <- data.frame(Z, X, U, D)
fit <- coxph(formula=Surv(U, D)~X+Z+X*Z, data=dd, method="breslow")
fit.std <- stdCoxph(fit=fit, data=dd, X="X", x=seq(-1,1,0.5), t=1:5)
print(summary(fit.std, t=3))
plot(fit.std)

}
