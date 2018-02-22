#This file contains code for the paper
#"Estimation of causal effect measures with the R-package stdReg"

#-----AF; logistic regression-----
  
fit <- glm(formula=lbw~smoker+race+age, family="binomial", data=clslowbwt)
options(show.signif.stars=FALSE)
summary(fit)
fit.std <- stdGlm(fit=fit, data=clslowbwt, X="smoker",  x=c(NA,0), clusterid="id")
summary(fit.std)
AF <- function(est){
  p <- est[1]
  p0 <- est[2]
  af <- c(1-p0/p)
  return(af)
}
AFest <- AF(fit.std$est)
AFest
confint(object=fit.std, fun=AF)

#-----AF; Cox PH regression-----

rott2$nochemo <- as.numeric(rott2$chemo=="no")
fit <- coxph(Surv(rf,rfi)~nochemo+age+meno+size+factor(grade)+
	I(exp(-0.12*nodes))+pr+er, data=rott2, method="breslow")
summary(fit)
t <- 10:60
fit.std <- stdCoxph(fit=fit, data=rott2, X="nochemo", t=10:60, x=c(NA,0))
AF <- function(est){
  p <- 1-est[, 1]
  p0 <- 1-est[, 2]
  af <- 1-p0/p
  return(af)
} 
AFest <- AF(fit.std$est)
AFci <- confint(object=fit.std, fun=AF)
plot(t, AFest, type="l", ylab="AF", xlab="time since diagnosis (months)",
 ylim=c(0.05,0.3))
matlines(t, AFci, lty="dashed", col="black")

#-----NNT; logistic regression-----

clslowbwt$nosmoke <- 1-clslowbwt$smoker
fit <- glm(formula=lbw~nosmoke+race+age, family="binomial", data=clslowbwt)
fit.std <- stdGlm(fit=fit,	X="nosmoke", data=clslowbwt, x=c(NA,1), 
  clusterid="id", subsetnew=nosmoke==0)
summary(fit.std)
NNT <- function(est){
  p <- est[1]
  p1 <- est[2]
  nnt <- 1/(p-p1)
  return(nnt)
}
NNTest <- NNT(fit.std$est)
NNTest
confint(object=fit.std, fun=NNT, type="log")

#-----NNT; Coxph regression-----

fit <- coxph(Surv(rf,rfi)~chemo+age+meno+size+factor(grade)+
	I(exp(-0.12*nodes))+pr+er, data=rott2, method="breslow")
fit.std <- stdCoxph(fit=fit, data=rott2, X="chemo", t=10:60, x=c(NA,"yes"), 
  subsetnew=chemo=="no")
NNT <- function(est){
  p <- 1-est[, 1]
  p1 <- 1-est[, 2]
  nnt <- 1/(p-p1)
  return(nnt)
}
NNTest <- NNT(fit.std$est)
NNTci <- confint(object=fit.std, fun=NNT)
plot(t, NNTest, type="l", ylab="NNT (patients)", xlab="time since diagnosis (months)",
 ylim=c(5,100))
matlines(t, NNTci, lty="dashed", col="black")

#-----RERI; logistic regression-----

clslowbwt$smokerrace <- factor(paste0(clslowbwt$smoker, as.numeric(clslowbwt$race=="2. Black")))
fit <- glm(formula=lbw~smokerrace+age, family="binomial", data=clslowbwt, subset=race!="3. Other")
fit.std <- stdGlm(fit=fit, X="smokerrace", data=clslowbwt, clusterid="id")
RERI <- function(est){
  p00 <- est[1]
  p01 <- est[2]
  p10 <- est[3]
  p11 <- est[4]
  reri <- (p11-p10-p01+p00)/p00
  return(reri)
}
RERIest <- RERI(fit.std$est)
RERIest
confint(object=fit.std, fun=RERI)

#-----RERI; Cox PH regression-----

rott2$nochemo <- as.numeric(rott2$chemo=="no")
rott2$chemograde <- factor(paste0(rott2$nochemo, rott2$grade-2))
fit <- coxph(Surv(rf,rfi)~chemograde+age+meno+size+
	I(exp(-0.12*nodes))+pr+er, data=rott2, method="breslow")
fit.std <- stdCoxph(fit=fit, data=rott2, X="chemograde", t=10:60)
RERI <- function(est){
  p00 <- 1-est[, 1]
  p01 <- 1-est[, 2]
  p10 <- 1-est[, 3]
  p11 <- 1-est[, 4]
  reri <- (p11-p10-p01+p00)/p00
  return(reri)
}
RERIest <- RERI(fit.std$est)
RERIci <- confint(object=fit.std, fun=RERI)
plot(t, RERIest, type="l", ylab="RERI", xlab="time since diagnosis (months)",
 ylim=c(-0.1,0.6))
matlines(t, RERIci, lty="dashed", col="black")





