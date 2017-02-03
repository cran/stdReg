\name{stdParfrailty}
\alias{stdParfrailty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Regression standardization in shared gamma-Weibull frailty models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{stdParfrailty} performs regression standardization in shared gamma-Weibull frailty models,
at specified values of the exposure, over the sample covariate distribution.
Let \eqn{T}, \eqn{X}, and \eqn{Z} be the survival outcome, the exposure, and a
vector of covariates, respectively. \code{stdParfrailty} uses a fitted Cox 
proportional hazards model to estimate the standardized 
survival function \eqn{\theta(t,x)=E\{S(t|X=x,Z)\}}, where \eqn{t} is a specific value of \eqn{T}, 
\eqn{x} is a specific value of \eqn{X}, and the expectation is over the marginal 
distribution of \eqn{Z}.     
}
\usage{
stdParfrailty(fit, data, X, x, t, clusterid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
%%     ~~Describe \code{formula} here~~
an object of class \code{"parfrailty"}, as returned by the \code{parfrailty} function 
  in the \code{stdReg} package..  
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
the standardized survival function. If \eqn{X} is binary (0/1) or
a factor, then \code{x} defaults to all values of \eqn{X}. If \eqn{X} is numeric,
then \code{x} defaults to the mean of \eqn{X}. If \code{x} is set to \code{NA},
then \eqn{X} is not altered. This produces an estimate of the marginal survival 
function \eqn{S(t)=E\{S(t|X,Z)\}}. 
}
  \item{t}{
%%     ~~Describe \code{t} here~~
an optional vector containing the specific values of \eqn{T} at which to estimate 
the standardized survival function. It defaults to all the observed event times
in \code{data}. 
}
  \item{clusterid}{
%%     ~~Describe \code{clusters} here~~
an string containing the name of the cluster identification variable.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{stdParfrailty} assumes that a shared gamma-Weibull frailty model 
\deqn{\lambda(t_{ij}|X_{ij},Z_{ij})=\lambda(t_{ij};\alpha,\eta)U_iexp\{h(X_{ij},Z_{ij};\beta)\}}
has been fitted, with parametrization as descibed in the help section for \code{parfrailty}.
Integrating out the gamma frailty gives the survival function
\deqn{S(t|X,Z)=[1+\phi\Lambda_0(t;\alpha,\eta)exp\{h(X,Z;\beta)\}]^{-1/\phi},}
where \eqn{\Lambda_0(t;\alpha,\eta)} is the cumulative baseline hazard
\deqn{(t/\alpha)^{\eta}.}
The ML estimates of \eqn{(\alpha,\eta,\phi,\beta)} are used to obtain 
estimates of the survival function \eqn{S(t|X=x,Z)}:
\deqn{\hat{S}(t|X=x,Z)=[1+\hat{\phi}\Lambda_0(t;\hat{\alpha},\hat{\eta})exp\{h(X,Z;\hat{\beta})\}]^{-1/\hat{\phi}}.} 
For each \eqn{t} in the \code{t} argument and for each \eqn{x} in the \code{x} argument, 
these estimates are averaged across all subjects (i.e. all observed values of \eqn{Z})
to produce estimates 
\deqn{\hat{\theta}(t,x)=\sum_{i=1}^n \hat{S}(t|X=x,Z_i)/n.} 
The variance for \eqn{\hat{\theta}(t,x)} is obtained by the sandwich formula. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{"stdParfrailty"} is a list containing 
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
%% ~put references to the literature/web site here ~
Makuch R.W. (1982). Adjusted survival curve estimation using covariates.
\emph{Journal of Chronic Diseases} \bold{35}, 437-443.

Chang I.M., Gelman G., Pagano M. (1982). Corrected group prognostic curves
and summary statistics. \emph{Journal of Chronic Diseases} \bold{35}, 669-674.

Gail M.H. and Byar D.P. (1986). Variance calculations for direct adjusted survival
curves, with applications to testing for no treatement effect. \emph{Biometrical Journal}
\bold{28}(5), 587-599. 
}
\author{
%%  ~~who you are~~
Arvid Sjolander
}
\note{
%%  ~~further notes~~
Standardized survival functions are sometimes referred to as (direct) adjusted
survival functions in the literature.

\code{stdParfrailty} does not currently handle time-varying exposures or covariates. 

\code{stdParfrailty} internally loops over all values in the \code{t} argument. Therefore,
the function will usually be considerably faster if \code{length(t)} is small.

The variance calculation performed by \code{stdParfrailty} does not condition on 
the observed covariates \eqn{\bar{Z}=(Z_1,...,Z_n)}. To see how this matters, note that 
\deqn{var\{\hat{\theta}(t,x)\}=E[var\{\hat{\theta}(t,x)|\bar{Z}\}]+var[E\{\hat{\theta}(t,x)|\bar{Z}\}].} 
The usual parameter \eqn{\beta} in a Cox proportional hazards model does not depend 
on \eqn{\bar{Z}}. Thus, \eqn{E(\hat{\beta}|\bar{Z})} is 
independent of \eqn{\bar{Z}} as well (since \eqn{E(\hat{\beta}|\bar{Z})=\beta}), so that the 
term \eqn{var[E\{\hat{\beta}|\bar{Z}\}]} in the corresponding variance decomposition 
for \eqn{var(\hat{\beta})} becomes equal to 0. However, \eqn{\theta(t,x)} depends 
on \eqn{\bar{Z}} through the average over the sample distribution for \eqn{Z}, 
and thus the term \eqn{var[E\{\hat{\theta}(t,x)|\bar{Z}\}]} is not 0, unless one 
conditions on \eqn{\bar{Z}}. The variance calculation by Gail and Byar (1986) ignores this term,
and thus effectively conditions on \eqn{\bar{Z}}.      
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

require(survival)

#simulate data
n <- 1000
m <- 3
alpha <- 1.5
eta <- 1
phi <- 0.5
beta <- 1
id <- rep(1:n,each=m)
U <- rep(rgamma(n,shape=1/phi,scale=phi), each=m)
X <- rnorm(n*m)
#reparametrize scale as in rweibull function
weibull.scale <- alpha/(U*exp(beta*X))^(1/eta)
T <- rweibull(n*m, shape=eta, scale=weibull.scale)

#right censoring
C <- runif(n*m,0,10)
D <- as.numeric(T<C)
T <- pmin(T,C)
  
#strong left-truncation
L <- runif(n*m,0,2)
incl <- T>L
incl <- ave(x=incl, id, FUN=sum)==m
dd <- data.frame(L, T, D, X, id)
dd <- dd[incl, ]  
 
fit <- parfrailty(formula=Surv(L, T, D) ~ X, data=dd, clusterid="id")
fit.std <- stdParfrailty(fit=fit, data=dd, X="X", x=seq(-1,1,0.5), t=1:5, clusterid="id")
print(summary(fit.std, t=3))
plot(fit.std)

}
