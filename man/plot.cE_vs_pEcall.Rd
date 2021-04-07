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
\examples{
shape1F=0.3
shape2F=0.7
shape1E=0.12
shape2E=0.88
res=plotBDP2(x="cE",y="PEcall",n=30,interim.at=15,
             pF=0.3,cF=0.01,pE=0.12,cE=seq(.5,1,by=.01),p0=0.3,p1=0.12,
             shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
             col=c("green","red"),cex.sub=.8)
plot(res)
}
