\name{pFstopEcall}
\alias{pFstopEcall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
operating characteristics of a single-arm trial with a binary endpoint with futility stopping and calling efficacy at interim
}
\description{
Determines the operating characteristics of a single-arm trial with a binary endpoint (success - failure) on the basis of analytical derivations.  The design allows for futility stopping and evaluates efficacy at interim.
Outputs are 
\itemize{
\item         probability of calling efficacy (P.effic)  at each interim analysis 
\item         cumulative probability of calling efficacy up to the interim (P.effic.cum).
\item          probability of futility stopping (P.futil)  at each interim analysis 
\item          cumulative stopping probability up to the interim (P.futil.cum).
}
}
\usage{
pFstopEcall(p, vn.int, v.critE, v.critF)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{true response rate
%%     ~~Describe \code{p} here~~
}
  \item{vn.int}{vector of sample sizes at the interim analyses 
%%     ~~Describe \code{vn.int} here~~
}
  \item{v.critE}{ vector of critical boundaries for calling efficacy at the interim analyses (call efficacy if number of successes >= boundary)
%%     ~~Describe \code{v.critE} here~~
}
  \item{v.critF}{vector of critical boundaries for futility stopping at the interim analyses (stop for futility if number of successes <= boundary). Choose boundary=-1 if no stopping for futility is allowed.
%%     ~~Describe \code{v.critF} here~~
}
}

\examples{
pFstopEcall(p=0.3,vn.int=c(3,9,13,18,20),v.critE=4:8, v.critF=c(0,1,2,3,3))
}
