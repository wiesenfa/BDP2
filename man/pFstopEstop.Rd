\name{pFstopEstop}
\alias{pFstopEstop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Operating characteristics of a single-arm trial with a binary endpoint with stopping for futility and stopping for efficacy}
\description{
Determines the operating characteristics of a single-arm trial with a binary endpoint (success - failure) on the basis of analytical derivations.  The design allows for stopping for futility and stopping for efficacy.
Outputs are 
\itemize{
\item          probability of stopping for efficacy (P.effic)  at each interim analysis 
\item          cumulative probability of stopping for efficacy up to the interim (P.effic.cum).
\item          probability of futility stopping (P.futil)  at each interim analysis 
\item         cumulative stopping probability up to the interim (P.futil.cum).
}
}
\usage{
pFstopEstop(p, vn.int, v.critE, v.critF)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{true response rate
%%     ~~Describe \code{p} here~~
}
  \item{vn.int}{vector of sample sizes at the interim analyses 
%%     ~~Describe \code{vn.int} here~~
}
  \item{v.critE}{vector of critical boundaries for calling efficacy at the interim analyses (stop for efficacy if number of successes >= boundary)
%%     ~~Describe \code{v.critE} here~~
}
  \item{v.critF}{ vector of critical boundaries for futility stopping at the interim analyses (stop for futility if number of successes <= boundary). Choose boundary=-1 if no stopping for futility is allowed.
%%     ~~Describe \code{v.critF} here~~
}
}
\examples{
pFstopEstop(p=0.3,vn.int=c(3,9,13,18,20),v.critE=4:8, v.critF=c(0,1,2,3,3))
}
