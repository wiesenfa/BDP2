\name{crit_general}
\alias{crit_general}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates CritBoundaries for different types of futility criteria
}
\description{
Calculates CritBoundaries for different types of futility criteria
}
\usage{
crit_general(n, p0, p1, vn.int, alpha, crit, type=5, shape1 = 1, shape2 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
       
  \item{n}{sample size at the final analysis
%%     ~~Describe \code{n} here~~
}
  \item{p0}{success rate corresponding to H0
%%     ~~Describe \code{p0} here~~
}
  \item{p1}{success rate corresponding to H1 (p1 > p0)
%%     ~~Describe \code{p1} here~~
}
  \item{vn.int}{vector of sample sizes at the interim analyses (the vector may be equal to 1:(n-1) = continuous monitoring of futility)
}
  \item{alpha}{nominal probability of type 1 error used for the final test
}
  \item{crit}{critical level(s) of predictive/conditional power, (posterior)
probabilities, rates, or patient numbers used for early termination.
crit translates into futility boundaries (maximum number of successes
leading to early termination).
If simple rates are used for monitoring futility (type=7), crit is the
critical success rate.
Rates or probabilites must be input as percentages. In case of type 10
analyses, crit must be  a vector of numbers of successes indicating futility
bounds at each (interim or final)analysis.

}
  \item{type}{see \code{\link{BDP2_simulate}} for details. Type of futility analysis: 1=predictive power; 2=conditional power under H1;
 3=conditional power under the MLE; 4=posterior probability of a success
 rate > p0;  5=posterior probability of a success rate >= p1;
 6=predictive probability of reaching a "positive" final study result,
 where "positive" is defined in terms of the posterior probability of a
 success rate > p0. (Here, the ciritcal level must be >= 1- alpha).
 7=estimated success rate; 8=p-values; 9="p-values under H1"; 10=user-
 defined bounds (a vector of critical numbers of successes for each analysis)

%%     ~~Describe \code{type} here~~
}
  \item{shape1}{for prior distribution
%%     ~~Describe \code{shape1} here~~
}
  \item{shape2}{for prior distribution
%%     ~~Describe \code{shape2} here~~
}
}
