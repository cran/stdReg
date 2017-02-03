\name{parfrailty}
\alias{parfrailty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fits shared gamma-Weibull frailty models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{parfrailty} fits shared gamma-Weibull frailty models. It is specifically 
designed to work with the function \code{stdParfrailty}, which performs regression 
standardization in shared gamma-Weibull frailty models.   
}
\usage{
parfrailty(formula, data, clusterid, init)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
an object of class "\code{formula}", on the same format as accepted by the 
\code{coxph} function in the \code{survival} package.  
}
  
  \item{data}{
%%     ~~Describe \code{data} here~~
a data frame containing the variables in the model. 
}
  
  \item{clusterid}{
%%     ~~Describe \code{clusters} here~~
an string containing the name of a cluster identification variable. 
}
  \item{init}{
%%     ~~Describe \code{case.control} here~~
an optional vector of initial values for the model parameters.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{parfrailty} fits the shared gamma-Weibull frailty model 
\deqn{\lambda(t_{ij}|C_{ij})=\lambda(t_{ij};\alpha,\eta)U_iexp\{h(C_{ij};\beta)\},}
where \eqn{t_{ij}} and \eqn{C_{ij}} are the survival time and covariate vector
for subject \eqn{j} in cluster \eqn{i}, respectively. \eqn{\lambda(t;\alpha,\eta)} is the 
Weibull baseline hazard function
\deqn{\eta t^{\eta-1}\alpha^{-\eta},}  
where \eqn{\eta} is the shape parameter and \eqn{\alpha} is the scale parameter.
Note that this is a different parametrization than in the \code{rweibull} function
in the \code{stats} package. \eqn{U_i} is the unobserved frailty term for 
cluster \eqn{i}, which is assumed to have a gamma distribution with scale = 
1/shape = \eqn{\phi}. \eqn{h(X;\beta)} is the regression function as specified
by the \code{formula} argument, parametrized by a vector \eqn{\beta}. The 
ML estimates \eqn{\{log(\hat{\alpha}),log(\hat{\eta}),log(\hat{\phi}),
\hat{\beta}\}} are obtained by maximizing the marginal (over \eqn{U}) likelihood.  
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of class \code{"parfrailty"} is a list containing:
\item{est}{
 the ML estimates \eqn{\{log(\hat{\alpha}),log(\hat{\eta}),
 log(\hat{\phi}),\hat{\beta}\}}.
}

\item{vcov}{
  the variance-covariance vector of the ML estimates.
  }
  
\item{score}{
a matrix containing the cluster-specific contributions to the ML score equations. 
}
}
\references{
%% ~put references to the literature/web site here ~
Van den Berg G.J., Drepper B. (2016). Inference for shared-frailty survival 
models with left-truncated data. \emph{Econometric Reviews}, 35(6), 1075-1098.
}
\author{
%%  ~~who you are~~
Arvid Sjolander and Elisabeth Dahlqwist. 
}
\note{
%%  ~~further notes~~
If left truncation is present, it is assumed that it is strong left truncation. 
This means that, even if the truncation time may be subject-specific, the whole 
cluster is unobserved if at least one subject in the cluster dies before 
his/her truncation time. If all subjects in the cluster survive beyond their 
subject-specific truncation times, then the whole cluster is observed 
(Van den Berg and Drepper, 2016). 

  
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
print(summary(fit))

}
