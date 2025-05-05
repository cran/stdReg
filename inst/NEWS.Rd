\name{NEWS}
\encoding{UTF-8}
\title{News for Package \pkg{stdReg}}
\section{Changes in version 3.4.2 (2025-05-05)}{
  \subsection{Bug fixes}{
    \itemize{
      \item The functionality for computing standard errors in stdCoxph has been updated,
      to reflect an update that was made in the survival package version 3.2-9
      for handling centering of binary covariates.  
    }
  } 
}
\section{Changes in version 3.4.1 (2021-05-11)}{
  \subsection{Bug fixes}{
    \itemize{
      \item A minor bug has been corrected in the handling of the argument subsetnew.  
    }
  } 
}
\section{Changes in version 3.4.0 (2019-03-21)}{
  \subsection{Bug fixes}{
    \itemize{
      \item An error has been corrected in the summary.parfrailty function.  
    }
  } 
  \subsection{User-visible changes}{
    \itemize{
      \item \code{parfrailty} now returns the \code{data} argument.
      \item Added output elements \code{call} and \code{input} to 
        \code{stdCoxph}, \code{stdGee}, \code{stdGlm} and \code{stdParfrailty}. 
      \item The \code{est}, \code{score} and \code{hessian} elements of 
      \code{class("parfrailty")} have been named.  
      \item A help entry for \code{summary.stdGlm} has been added.
      \item Informative error messages have been added in places.
    }
  } 
  \subsection{Non-visible changes}{
    \itemize{
      \item Less internal calls to \code{coxph.detail} in \code{stdCoxph},
      which makes the function slightly faster.
      \item Replaced \code{class(fit)=="name"} with 
      \code{inherits(x=fit, what="name")} in places, which should be safer.
      \item Minor changes in the internal handling of formula arguments,
      to avoid errors if formulas are given implicitly to models,
      e.g. \code{f <- y~x}, \code{fit <- glm(formula=f)}  
    }
  } 
}
\section{Changes in version 3.3.0 (2019-03-15)}{
  \subsection{Bug fixes}{
    \itemize{
      \item A minor bug has been corrected in the calculation of \code{vcov}
      for \code{stdCox}. This may have a minor effect on the estimated variance-covariance
      matrix.  
    }
  } 
}

