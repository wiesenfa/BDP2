\name{plotBDP2}
\alias{plotBDP2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plots
}
\description{
Output of desired plots
}
\usage{
plotBDP2(x = c("n", "k", "ptrue", "cE", "cF"), 
          y = c("Prob0Successes", "PostProb0or1Successes", "bFbE", 
                "PEcall_p0_p1", "PEstop_p0_p1", 
                "PFstopEcall", "PFstopEstop",
                "PEcall", "PEstop", "PFstop", "PFstopEstop", 
                "ExpectedNumber", 
                "PredictivePower"), 
          n, interim.at, ptrue, 
          pF, cF, pE, cE, p0, p1, Estop=FALSE,
          shape1F, shape2F, shape1E = NULL, shape2E = NULL, 
          col = c("green", "red"), cex.legend=1, add = FALSE, 
          show=TRUE, progress = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
character string specifying what is given by the x axis

\itemize{
\item "n" for number of patients at final or at first interim, depending on plot
\item "ptrue" for the true response rate
\item "k" number of successes at interim (only for \code{y=="PredictivePower"})
\item "cE" for the critical level of posterior probabilities used for declaring efficacy
\item "cF" for the critical level of posterior probabilities used for declaring futility
}
}
  \item{y}{
character string specifying what is given by the y axis
\itemize{
\item "Prob0Successes" for plot of probability of 0 successes out of n at first interim, i.e. (1-p)^n. Two curves are generated: one for p0 and one for p1. Choose x="n" and set n to a vector of number patients at first interim.
\item "PostProb0or1Successes" generates two curves of posterior probability that reponse rate exceeds cF. One given 0 successes and one given 1 success observed in n, i.e. P(p>pF| 0 successes out of n) and P(p>pF| 1 success out of n). Choose x="n" and set n to a vector of number patients at first interim.
\item "bFbE" for plot of boundaries for futility bF and for efficacy bE in terms of number of successes per n number of patients. Futility stop if number of successes <= bF, call efficacy if number of successes >= bE. Choose x="n" and set n to the maximal number of patients (not a vector).
\item "PEcall_p0_p1" for plot of probability to call efficacy at final (i.e. P(p>pE|Data)>=cE ), evaluated for data generated with control response rate p0 (corresponds to type I error) and for data generated with target response rate p1 (corresponds to power) for varying number of patients at final. Choose x="n" and set n to a vector of number patients at final.
\item "PEstop_p0_p1" as above but with stopping for efficacy instead of calling efficacy 
\item "PFstopEcall" for plotting operating characteristics for a given response rate as a function of number of patients at final, n. Shows one curve (default: in red) for the probability of stopping for futility up to final analysis and one curve (default: in green) for the probability of calling efficacy at final analysis. For the setting ptrue=p0, the curves show the probability of true stopping (default: in red) and type I error (default: in green). For the setting ptrue=p1, the curves show the probability of false stopping (default: in red) and power (default: in green). Choose x="n" and set n to a vector of number patients at final.
\item "PFstopEstop" as above but with stopping for efficacy instead of calling efficacy 
\item "PEcall" for plot of probability to call efficacy at final (i.e. P(p>pE|Data)>=cE ). 

\itemize{
\item For x="ptrue" then this is the power function. 
\item For x="cE" this gives plots of probability to call efficacy at final, evaluated for data generated with control response rate p0 (corresponds to type I error) and for data generated with target response rate p1 (corresponds to power).
}
\item "PEstop" 
\itemize{
\item For x="ptrue" then this is the power function. 
\item For x="cE" this gives plots of probability to call efficacy at final, evaluated for data generated with control response rate p0 (corresponds to type I error) and for data generated with target response rate p1 (corresponds to power).
}


\item "PFstop" for plot of cumulative probability to stop for futility up to final.
\itemize{
\item For x="ptrue" this gives the futility stopping probability as function of ptrue. 
\item For x="cF" this gives 2 curves, evaluated for data generated with control response rate p0 (corresponds to true stopping probability, default: in green) and for data generated with target response rate p1 (corresponds to false stopping probability, default: in red).
}
\item "ExpectedNumber" for expected number of patients in the trial. Choose x= "ptrue". 
      Takes stopping for efficacy into account if \code{Estop==TRUE}
\item "PredictivePower" for predictive power (only for \code{x=="k"})
}
}
  \item{n}{sample size at the final analysis, vector if \code{x=="n"}}
  \item{interim.at}{vector of sample sizes at the interim analyses}
  \item{ptrue}{true (assumed) response rate used for simulating the trial, vector if \code{x=="ptrue"}}
  \item{pF}{response rate used for the futility criterion P(p>pF|Data) < cF (may be identical to pE)}
  \item{cF}{critical level of posterior probabilities used for declaring futility, vector if \code{x=="cF"} }
  \item{pE}{response rate used for the efficacy criterion P(p>pE|Data) >= cE}
  \item{cE}{critical level of posterior probabilities used for declaring efficacy, vector if \code{x=="cE"}}
  \item{p0}{response rate corresponding to H0}
  \item{p1}{response rate corresponding to H1 (p1 > p0)}
  \item{shape1F}{first parameter of the Beta prior for futility analysis}
  \item{shape2F}{second parameter of the Beta prior for futility analysis}
  \item{shape1E}{first parameter of the Beta prior for efficacy analysis}
  \item{shape2E}{second parameter of the Beta prior for efficacy analysis }
  \item{Estop}{Stop for efficacy? Defaults to \code{FALSE}. Only relevant if \code{y=="ExpectedNumber"}. }
  \item{col}{line color, for some plots vector of length 2.}
  \item{add}{add line to existing plot. Only supported if \code{x=="ptrue"}}
  \item{show}{show plot (otherwise computed objects are invisibly returned)}
  \item{progress}{only used by shiny app}
  \item{cex.legend}{size of legend text relative to cex}
  \item{\dots}{additional arguments passed to \code{plot.default()}}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Kopp-Schneider, A., Wiesenfarth, M., Witt, R., Edelmann, D., Witt, O. and Abel, U. (2018).\cr Monitoring futility and efficacy in phase II trials with Bayesian
posterior distributions - a calibration approach. 
\emph{Biometrical Journal}, to appear.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# See vignette for more details and examples
pF=0.3
pE=0.12
shape1F=0.3
shape2F=0.7
shape1E=0.12
shape2E=0.88
cF=0.01
cE=0.9
nvec=c(18:40)
interim.at=c(10,20,30)

# Type I error and probability of true stopping for the uninteresting response rate. 
ptrue=0.12
plotBDP2(x="n", y="PFstopEcall",
         n =nvec, interim.at = interim.at,
         pF=pF,cF=cF,pE=pE,cE=cE,ptrue=ptrue,
         shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E)


# Power and probability of false stopping for the target response rate. 
ptrue=0.3
plotBDP2(x="n", y="PFstopEcall",
         n =nvec, interim.at = interim.at,
         pF=pF,cF=cF,pE=pE,cE=cE,ptrue=ptrue,
         shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
