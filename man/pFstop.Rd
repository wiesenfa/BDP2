\name{pFstop}
\alias{pFstop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Operating characteristics of a single-arm trial with a binary endpoint with futility stopping
}
\description{
Determines the operating characteristics of a single-arm trial with a binary endpoint (success - failure) on the basis of analytical derivations. The design allows for futility stopping.
Outputs are
\itemize{
\item probability of futility stopping (P.futil) at each interim analysis 
\item cumulative stopping probability up to the interim (P.futil.cum).
}
}
\usage{
pFstop(p, vn.int, v.crit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{p}{true response rate
%%     ~~Describe \code{p} here~~
}
  \item{vn.int}{vector of sample sizes at the interim analyses 
%%     ~~Describe \code{vn.int} here~~
}
  \item{v.crit}{vector of critical boundaries for futility stopping at the interim analyses (stop for futility if number of successes <= boundary). Choose boundary=-1 if no stopping is allowed.
%%     ~~Describe \code{v.crit} here~~
}
}
\examples{
pFstop(p=0.3,vn.int=c(3,9,13,18,20),v.crit=c(0,1,2,3,3))
}
