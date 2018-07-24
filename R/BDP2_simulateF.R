BDP2_simulate=function (n, vn.int, p, p0, p1, alpha=0.05, crit, type=5, nsim,shape1=1,shape2=1) {

  # Determines the operating characteristics of a single-arm trial with a binary
  # endpoint (success - failure) and interim futility analyses.
  # The user can choose among 10 futility criteria, which are based on predictive
  # or conditional power (the latter either assuming H1 or the MLE), posterior
  # or predictive probabilities, tail probabilites (under H0 or H1), constant success
  # rates, or arbitrary user-defined futility bounds.

  # Assumptions:
  # Endpoint (success/no success) data available for all study patients.
  # In case of Bayesian analysis: Beta-binomial model. Prior distribution = Beta(prior.mean*prior.sampleSize,(1-prior.mean)*prior.sampleSize) (uniform if prior.mean=0.5,prior.sampleSize=2)
  # In the following, for all Bayesian analyses, let shape parameters of the prior distribution be shape1=prior.mean*prior.sampleSize, shape2=(1-prior.mean)*prior.sampleSize.
  # In case of predictive or conditional power, tail probabilities, or rates: One-sided
  # testing in the final analysis.
   # Some methodological details on the 10 types of futility criteria:

  # Type 1: Predictive power.
  # Given the results of the interim analysis the predictive power at the final analysis
  # (n patients, critical number of successes k.crit) is P(X >= k.crit - k.int), where
  # X  follows a beta-binomial distribution with parameters n'= n - n.int, a = k.int + shape1,
  # and b= n.int - k.int + shape2.

  # Type 2,3: Conditional power.
  # The conditional power at the interim analysis is P(X >= k.crit - k.int), where X
  # follows a binomial distribution with parameters n'= n - n.int, and success probability
  # either equal to p1 (futility analysis type 2) or to the estimated success rate (MLE)
  # at the interim analysis (type 3.)

  # Type 4,5:Posterior probabilities.
  # The posterior distribution at interim analysis with n.int  patients and k.int
  # successes is Beta(k.int + shape1, n + shape2 - k.int)
  # Type  4: Futility is declared if the posterior probability P(true success rate > p0)
  # is < crit. (Here, crit must be large,e.g. 70%).
  # Type  5: Futility is declared if the posterior probability P(true success rate >= p1)
  # is < crit. (Here, crit must be small, e.g. 10%).

  # Type 6: Predictive probability combined with posterior probability.
  # Futility is declared if the posterior predictive probability that the study
  # will be a success is < crit (e.g. 10%). Here, the success is defined by the total number
  # of successes in the trial yielding a posterior probability of at least 1 - alpha
  # (when evaluated in the final analysis) that the true success rate is > p0.

  # Type 7: Estimated success rates.
  # Futility is declared if the success rate is smaller than a fixed benchmark crit. The
  # final analysis is test-based.

  # Type 8,9: Tail probabilites under H0,H1.
  # Type 8: The futility criterion uses an alpha level crit that is constant across all
  # interim analysis. The final analysis is test-based. Futility is declared if the p-value
  # (upper tail) is >= crit.
  # Type 8 futility analyses should only be used if the number of patients at the first
  # interim analysis is not too low (say, at least 5 to 10). The value of crit is not
  # identicalt o the alpha level used in the final test. Generally, a fairly high value of
  # crit will be appropriate (e.g. 70%).
  # Type 9: Similar to type 8, but with lower-tail probabilites calculated under H1
  # ("p-values under H1").
  # I.e., futility is declared if, under H1, the probability of obtaining at most as many
  # successes as the observed number is < crit ("observed number of successes "too low" to
  # be compatible with H1 at one-sided significance level = crit).
  # Generally, a small value of crit (e.g. 5% or 10%) should be chosen.

  # Type 10: User-defined boundaries.
  # Here, the futility boundaries (maximum numbers of successes leading to early termination)
  # are directly input by the user. crit ist the vector of these boundaries at each
  # (interim or final) analysis. The study is terminated if the number of successes is
  # at analysis no. m is <= the crit[m].

  # Args:
  # n:           sample size at the final analysis
  # p0:          success rate (not in %!) corresponding to H0
  # p1:          success rate (not in %!) corresponding to H1 (p1 > p0)
  # p:           true (assumed) success rate (%) used for simulating the trial
  # vn.int:      vector of sample sizes at the interim analyses (the vector may be equal
  #                to 1:(n-1) = continuous monitoring of futility)
  # alpha:       nominal probability of type 1 error (not in %!) used for the final test
  # crit:        critical level(s) of predictive/conditional power, (posterior)
  #                probabilities, rates, or patient numbers used for early termination.
  #                crit translates into futility boundaries (maximum number of successes
  #                leading to early termination).
  #                If simple rates are used for monitoring futility (type=7), crit is the
  #                critical success rate.
  #                Rates or probabilites must be input as percentages. In case of type 10
  #                analyses, crit must be  a vector of numbers of successes indicating futility
  #                bounds at each (interim or final)analysis.
  # type:        type of futility analysis: 1=predictive power; 2=conditional power under H1;
  #                3=conditional power under the MLE; 4=posterior probability of a success
  #                rate > p0;  5=posterior probability of a success rate >= p1;
  #                6=predictive probability of reaching a "positive" final study result,
  #                where "positive" is defined in terms of the posterior probability of a
  #                success rate > p0. (Here, the ciritcal level must be >= 1- alpha).
  #                7=estimated success rate; 8=p-values; 9="p-values under H1"; 10=user-
  #                defined bounds (a vector of critical numbers of successes for each analysis)
  # nsim:        number of simulation runs
  # shape1         for prior distribution
  # shape2         for prior distribution

  # Output:
  # Prints the results of the analysis.

  # Functions invoked:
  # CritBoundary, CritBoundaryPost, CritBoundaryTest, SimOC

  #--------------------------------------------------------------------------------------------
  # Input check:
  if (p0 >= p1) stop('p1 must be > p0.')
  if (!(type %in% 1:10)) stop('type must be an integer between 1 and 10')
  if (max(vn.int) >= n) stop('Interim analyses must be performed at patient numbers < n.')
  if (is.unsorted(vn.int) == TRUE) stop('vn.int must be sorted in ascending order')
  if ((type == 10) && (length(crit) != 1 + length(vn.int))) stop('incorrect length of vector crit')

    # Preliminary commands:

  if (type != 10) crit <- crit
  n.int <- length(vn.int)
  n.an <- n.int + 1
  vn.an <- c(vn.int, n)


  #--------------------------------------------------------------------------------------------
  # 1. Basic structure of the decision process (final analysis):

  # 1a. Critical values for testing or type 6 futility criterion used in the final analysis:

  # Critical value for the one-sided binomial test:

  k.crit.fin <- qbinom(size = n, prob = p0, p = alpha, lower.tail = FALSE) + 1

  # Critical value for the type 6 futility criterion:
  if (type == 6) {

    for (k in 1:n) {
      # k observed successes among n patients

      z <- pbeta(q = p0, shape1 = shape1 + k, shape2 = shape2 + n - k, lower.tail = FALSE)
      # = posterior probability (calculated from k successes observed among
      # n patients) that psuccess > p0.

      k.crit.fin <- k
      if (z > 1- alpha) break
      # k is high enough to declare the treatment successful"

    } # End loop over number of observed successes
  }

  # 1b, Power of the final test (type 1-3, 7-9 futility criteria):

  power.fin <- pbinom(size = n, q = k.crit.fin - 1, prob = p1, lower.tail = FALSE)

  if ((type <= 3) || (type %in% 7:9))  {
    power.fin <- round(100*power.fin, digits=1)
    cat('\n')
    cat('Exact power of the test used in the final analysis (%):  ', power.fin, '\n\n')

    if ((type <= 3) && (power.fin <= crit)) {
      stop('The critical predictive or conditional power must be lower than the nominal one')
    }
  }

  #--------------------------------------------------------------------------------------------
  # 2. Vectors of critical boundaries corresponding to interim analyses at each possible
  #    patient number and final analysis:

  # Critical boundaries for type 1,2,3,6. (Here, critical boundaries make use of the final
  # analysis):
  if ((type <= 3) || (type == 6)) {
    v.boundary <- CritBoundary(type, n, crit, k.crit.fin, p1,shape1=shape1,shape2=shape2) #shape parameters only used in type 6
  }

  # Critical boundaries for type 4,5 (Here, critical boundaries are based on the posterior
  # probabilites of the true success rates):
  if (type == 4) {
    v.boundary <- CritBoundaryPost(n, p0, crit,shape1=shape1,shape2=shape2)
  }
  if (type == 5) {
    v.boundary <- CritBoundaryPost (n, p1, crit,shape1=shape1,shape2=shape2)
  }

  # Critical boundaries for type 7 (futility declared if the success rate is < crit):
  if (type == 7) {
    v.boundary <- c(ceiling(crit*1:(n - 1)) - 1, k.crit.fin -1)
    # critical boundaries based on constant benchmark rate crit, except for the
    # final analysis
  }

  # Critical boundaries for type 8 (futility declared p < crit):
  if (type == 8) {
    v.boundary <- c(CritBoundaryTest (n, p0, crit, type), k.crit.fin - 1)
  }

  # Critical boundaries for type 9 (futility declared if, under H1, the probability
  # of obtaining at most the observed number of successes is < crit):
  if (type == 9) {
    v.boundary <- c(CritBoundaryTest (n, p1, crit, type), k.crit.fin - 1)
  }

  # Critical boundaries for type 10 (user-defined boundaries):
  if (type == 10) {
    v.boundary <- rep(-1, n)
    v.boundary[c(vn.an)] <- crit
  }
  v.crit <- v.boundary[vn.an]

  #--------------------------------------------------------------------------------------------

  # 3. Simulations determining the OC:

  # Type 1 error:
  alpha.est <- round(1 - (SimOC (n, p0, vn.int, v.boundary, nsim)[[2]])[n.an], digits = 4)

  # Power:
  power.est <- round(1 - (SimOC (n, p1, vn.int, v.boundary, nsim)[[2]])[n.an], digits = 4)

  # Probabilities of early stopping/null result for the assumed success rate p:
  SimRes <- SimOC (n, p, vn.int, v.boundary, nsim)

  # Data frame for output (early stopping):
  SimResDF <- data.frame(Int.Analysis = 1:n.an, Pat.number = vn.an, Crit.boundary = v.crit,
                         P.stop = SimRes[[1]], P.stop.cum = SimRes[[2]])

  # Expected number of patients:
  mean.npat <- SimRes[[3]]

  #--------------------------------------------------------------------------------------------
  # 4. Output:
  if ((type <= 3) || (type %in% 7:9))  {
    cat('\n')
    cat('Critical value of the test used in the final analysis:   ', k.crit.fin, '\n\n')
  }
  if (type == 6) {
    cat('\n')
    cat('Critical value of the futility criterion used in the final analysis:   ', k.crit.fin, '\n\n')
  }

  cat('\n\n')
  cat('#######################################################################\n\n')
  if (type == 1 ) cat('Futility analysis based on predictive power\n\n')
  if (type == 2 ) cat('Futility analysis based on conditional power calculated under H1\n\n')
  if (type == 3 ) cat('Futility analysis based on conditional power calculated under the MLE\n\n')
  if (type == 4 ) {
    cat('Futility analysis based on the posterior probability of\n')
    cat('success rate > p0\n\n')
  }
  if (type == 5 ) {
    cat('Futility analysis based on the posterior probability of\n')
    cat('success rate >= p1\n\n')
  }
  if (type == 6 ) {
    cat('Futility analysis based on the predictive probability of attaining a\n')
    cat('posterior probability of at least 1-alpha (when evaluated in the final\n')
    cat('analysis) that the success rate is > p0\n\n')
  }
  if (type == 7 ) cat('Futility analysis based on success rates\n\n')
  if (type == 8 ) {
    cat('Futility analysis based on p-values (upper tail)\n\n')
  }
  if (type == 9 ) {
    cat('Futility analysis based on lower tail probabilites under H1\n\n')
  }
  if (type == 10 ) {
    cat('Futility analysis based on user-defined bounds\n\n')
  }

  cat('#######################################################################\n\n')

  cat('Simulation results (Operating Characteristics of the design):\n\n')

  if ((type <= 3) || (type == 7)) {
    cat('  Estimated alpha level, accounting for futility analyses: ', alpha.est, '\n\n')
    cat('  Estimated power, accounting for futility analyses:       ', power.est, '\n\n')
  }

  cat('\n')
  cat('Critical boundaries for early stopping and probabilities of early stopping.\n')
  cat('The recruitment stops if the number of successes at an interim analysis\n')
  cat('is <= the critical boundary.\n')
  cat('The last row corresponds to the final analysis\n')
  cat('(in case of statistical testing: non-rejection of H0).\n\n')
  print(SimResDF)
  cat('\n')
  cat('Expected number of patients enrolled in the trial   :', mean.npat, '\n\n')

}
