\name{plot.cE_vs_pEcall}
\alias{plot.cE_vs_pEcall}
\alias{plot.cE_vs_pEstop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot objects returned by plotBDP2()
}
\description{
Plot objects returned by \code{\link{plotBDP2}} with arguments \code{x="cE"} and either \code{y="PEcall"} or \code{y="PEstop"} which can computationally be relatively expensive.
}
\usage{
\method{plot}{cE_vs_pEcall}(x, ...)
\method{plot}{cE_vs_pEstop}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object returned by \code{\link{plotBDP2}} with arguments \code{x="cE"} and either \code{y="PEcall"} or \code{y="PEstop"}
}
  \item{\dots}{
arguments passed to \code{\link{plot.default}}
}
}
