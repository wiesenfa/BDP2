\name{BDP2}
\alias{BDP2}
\alias{BDP2_Fstop}
\alias{BDP2_FstopEcall}
\alias{BDP2_FstopEstop}
\alias{BDP2_simulateEF_PostProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Operating characteristics of a single-arm trial with a binary
endpoint 
}
\description{
Determines the operating characteristics of a single-arm trial with a binary
endpoint (response, success) and interim efficacy and futility analyses.
Declaration of efficacy and futility (including possibly early stopping)
is based on the posterior probability that the true response rate
is at least pE , pF respectively.
}
\usage{
BDP2(n, interim.at, ptrue, 
      eff.stop = FALSE, 
      pF, cF, pE = NULL, cE = NULL, 
      shape1F, shape2F, shape1E = NULL, shape2E = NULL, 
      simulate = FALSE, nsim = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size at the final analysis}
  \item{interim.at}{vector of sample sizes at the interim analyses}
  \item{ptrue}{true (assumed) response rate used for analytical evaluations or simulating the trial}
  \item{eff.stop}{
    \code{FALSE}: No evaluation of efficacy. 
    \code{"call"}:  no stop for
    efficacy; in this case the program merely calculates the probability
    that the efficicacy criterion is satisfied (possibly triggering a
    notification of the DMC and the start of the planning of a subsequent
    trial).
    \code{stop}: the study ends if the efficacy criterion is reached at an interim analysis.
  }
  \item{pF}{response rate used for the futility criterion (may be identical to pE)}
  \item{cF}{critical level of posterior probabilities used for declaring futility }
  \item{pE}{response rate used for the efficacy criterion}
  \item{cE}{critical level of posterior probabilities used for declaring efficacy}
  \item{shape1F}{first parameter of the Beta prior for futility analysis}
  \item{shape2F}{second parameter of the Beta prior for futility analysis}
  \item{shape1E}{first parameter of the Beta prior for efficacy analysis}
  \item{shape2E}{second parameter of the Beta prior for efficacy analysis }
  \item{simulate}{\code{FALSE} for analytical evaluation and \code{TRUE} for simulation}
  \item{nsim}{number of simulation runs (only used if \code{simulate==TRUE})}
}
\details{

Assumptions:
Endpoint (response/no response) data available for all study patients.
Beta-binomial model. Prior distribution = Beta(shape1, shape2).

The posterior distribution at interim analysis with n.int  patients and k.int
successes is Beta(k.int + shape1F, n.int + shape2F - k.int) and Beta(k.int + shape1E, n.int + shape2E - k.int), respectively.
Efficacy is declared if the posterior probability P(true response rate > pE)
is >= cE.
Futility is declared if the posterior probability P(true success rate > pF)
is < cF.
cF, cE translate into futility/efficacy boundaries (maximum number of
responses leading to early termination for futility/ minimum number of responses
leading to declaring of, or early termination for, efficacy).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
