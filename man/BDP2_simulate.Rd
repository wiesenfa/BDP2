\name{BDP2_simulate}
\alias{BDP2_simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulated operating characteristics of a single-arm trial with a binary
endpoint for different types of futility criteria}
\description{
NOTE: Usually function \code{\link{BDP2}} will be preferred.


Determines the operating characteristics of a single-arm trial with a binary
endpoint (success - failure) and interim futility analyses.
The user can choose among 10 futility criteria, which are based on predictive
or conditional power (the latter either assuming H1 or the MLE), posterior
or predictive probabilities, tail probabilites (under H0 or H1), constant success
rates, or arbitrary user-defined futility bounds.

Assumptions:
Endpoint (success/no success) data available for all study patients.
In case of Bayesian analysis: Beta-binomial model. Prior distribution = Beta(shape1,shape2) (uniform if prior.mean=0.5,prior.sampleSize=2)
In case of predictive or conditional power, tail probabilities, or rates: One-sided
testing in the final analysis.
}
\usage{
BDP2_simulate(n, vn.int, p, p0, p1, 
              alpha = 0.05, crit, type = 5, nsim, 
              shape1 = 1, shape2 = 1)
}
\arguments{
  \item{n}{sample size at the final analysis
 }
  \item{vn.int}{
vector of sample sizes at the interim analyses (the vector may be equal
to 1:(n-1) = continuous monitoring of futility)
  }
  \item{p}{true (assumed) success rate used for simulating the trial
}
  \item{p0}{success rate corresponding to H0
}
  \item{p1}{success rate corresponding to H1 (p1 > p0)
}
  \item{alpha}{nominal probability of type 1 error used for the final test
}
  \item{crit}{
critical level(s) of predictive/conditional power, (posterior)
probabilities, rates, or patient numbers used for early termination.
crit translates into futility boundaries (maximum number of successes
leading to early termination).
If simple rates are used for monitoring futility (type=7), crit is the
critical success rate.
Rates or probabilites must be input as percentages (???). In case of type 10
analyses, crit must be  a vector of numbers of successes indicating futility
bounds at each (interim or final)analysis.
}
  \item{type}{
type of futility analysis: 1=predictive power; 2=conditional power under H1;
3=conditional power under the MLE; 4=posterior probability of a success
rate > p0;  5=posterior probability of a success rate >= p1;
6=predictive probability of reaching a "positive" final study result,
where "positive" is defined in terms of the posterior probability of a
success rate > p0. (Here, the ciritcal level must be >= 1- alpha).
7=estimated success rate; 8=p-values; 9="p-values under H1"; 10=user-
defined bounds (a vector of critical numbers of successes for each analysis)}
  \item{nsim}{number of simulation runs
}
  \item{shape1}{shape parameter for prior distribution
}
  \item{shape2}{shape parameter for prior distribution
}
}
\details{
Some methodological details on the 10 types of futility criteria:

Type 1: Predictive power.
Given the results of the interim analysis the predictive power at the final analysis
(n patients, critical number of successes k.crit) is P(X >= k.crit - k.int), where
X  follows a beta-binomial distribution with parameters n'= n - n.int, a = k.int + shape1,
and b= n.int - k.int + shape2.

Type 2,3: Conditional power.
The conditional power at the interim analysis is P(X >= k.crit - k.int), where X
follows a binomial distribution with parameters n'= n - n.int, and success probability
either equal to p1 (futility analysis type 2) or to the estimated success rate (MLE)
at the interim analysis (type 3.)

Type 4,5:Posterior probabilities.
The posterior distribution at interim analysis with n.int  patients and k.int
successes is Beta(k.int + shape1, n + shape2 - k.int)
Type  4: Futility is declared if the posterior probability P(true success rate > p0)
is < crit. (Here, crit must be large,e.g. 70\%).
Type  5: Futility is declared if the posterior probability P(true success rate >= p1)
is < crit. (Here, crit must be small, e.g. 10\%).

Type 6: Predictive probability combined with posterior probability.
Futility is declared if the posterior predictive probability that the study
will be a success is < crit (e.g. 10\%). Here, the success is defined by the total number
of successes in the trial yielding a posterior probability of at least 1 - alpha
(when evaluated in the final analysis) that the true success rate is > p0.

Type 7: Estimated success rates.
Futility is declared if the success rate is smaller than a fixed benchmark crit. The
final analysis is test-based.

Type 8,9: Tail probabilites under H0,H1.
Type 8: The futility criterion uses an alpha level crit that is constant across all
interim analysis. The final analysis is test-based. Futility is declared if the p-value
(upper tail) is >= crit.
Type 8 futility analyses should only be used if the number of patients at the first
interim analysis is not too low (say, at least 5 to 10). The value of crit is not
identicalt o the alpha level used in the final test. Generally, a fairly high value of
crit will be appropriate (e.g. 70\%).
Type 9: Similar to type 8, but with lower-tail probabilites calculated under H1
("p-values under H1").
I.e., futility is declared if, under H1, the probability of obtaining at most as many
successes as the observed number is < crit ("observed number of successes "too low" to
be compatible with H1 at one-sided significance level = crit).
Generally, a small value of crit (e.g. 5\% or 10\%) should be chosen.

Type 10: User-defined boundaries.
Here, the futility boundaries (maximum numbers of successes leading to early termination)
are directly input by the user. crit ist the vector of these boundaries at each
(interim or final) analysis. The study is terminated if the number of successes is
at analysis no. m is <= the crit[m].
}

\examples{
BDP2_simulate(n=30, vn.int=c(10,20), 
              p=0.2, p0=0.1, p1=0.3, crit=0.1, type=5, nsim=1000,shape1=1,shape2=1)
}

