crit=function (n, pThreshold,  vn.int, crit, shape1=1, shape2=1) {

  # CritBoundaries for P(p > pThreshold) < crit

  # Args:
  # n:           sample size at the final analysis
  # pThreshold:  threshold that the response rate p should exceed
  # vn.int:      vector of sample sizes at the interim analyses
  # crit:        threshold for posterior probability
  # shape1:      shape1 parameter for prior distribution
  # shape2:      shape2 parameter for prior distribution

  #--------------------------------------------------------------------------------------------
  # Input check:
  if (max(vn.int) >= n) stop('Interim analyses must be performed at patient numbers < n.')
  if (is.unsorted(vn.int) == TRUE) stop('vn.int must be sorted in ascending order')

  # Preliminary commands:
  vn.an <- c(vn.int, n)


  v.boundary <- CritBoundaryPost (n, pThreshold, crit,shape1=shape1,shape2=shape2)

  v.crit <- v.boundary[vn.an]
  v.crit
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

crit_general=function (n, p0, p1,  vn.int, alpha, crit, type=5, shape1=1 , shape2=1) {

  ##Calculates CritBoundaries for all 10 types of futility criteria:
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
  # shape1         for prior distribution
  # shape2         for prior distribution

  #--------------------------------------------------------------------------------------------
  # Input check:
  if (p0 >= p1) stop('p1 must be > p0.')
  if (!(type %in% 1:10)) stop('type must be an integer between 1 and 10')
  if (max(vn.int) >= n) stop('Interim analyses must be performed at patient numbers < n.')
  if (is.unsorted(vn.int) == TRUE) stop('vn.int must be sorted in ascending order')
  if ((type == 10) && (length(crit) != 1 + length(vn.int))) stop('incorrect length of vector crit')

  # Preliminary commands:

  if (type != 10) crit <- crit
  # p <- p/100
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
  v.crit
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

CritBoundary=function (type, n, crit, k.crit.fin, p1,shape1=1,shape2=1) {

  # CritBoundary

  # Vector of critical boundaries corresponding to interim analyses of type 1-3, 6,
  # and final analysis:
  # With an analysis performed at i patients, an observed number k of successes
  # leads to a negative decision (suspension/stopping of the trial or a "null result"
  # at the final analysis) iff k <= v.boundary[i]
  # The critical boundaries calculated for early stopping at interim analyses are based
  # on the critical boundaries (k.crit.fin) for success of the trial at the final analysis.

  # Args:
  # type:         type of futility analysis: 1=predictve power; 2= conditional power
  #                 under H1; 3=conditional power under the MLE;  6= predictive probabi-
  #                 lity of reaching a "positive" final study result in terms, where
  #                 "positive" is defined in terms of the posterior probability of a success
  #                 rate > p0 (must be >= 1- alpha).
  # n:            sample size at the final analysis
  # crit:         critical level (%) for probability of a positive judgement on the treatment
  #                 based on predicitve/conditional power or the predictive probability
  # k.crit.fin    critical value for trial success used in the final analysis
  # p1:           probability of success corresponding to H1 (needed for futility analyses
  #                 of type 2)

  # Returns:
  # Vector v.boundary of of critical boundaries for an analysis at i=1,...,n patients.

  # Functions invoked: CondPower, PredPower, CritBoundaryPost
  #--------------------------------------------------------------------------------------------
  v.boundary <- rep(-1, n)
  # Initially: no stopping for futility

  # Loop over patient numbers 1,...,n-1 (corresponding to number of patients avalaible for
  # interim analyses):

  for (i in 1:(n - 1)) {

    # Loop over the number of successes observed among the i patients
    for (k in 0:i) {

      if ((type == 1) || (type == 6)) {
        z <- PredPower(n, i, k, k.crit.fin,shape1=shape1,shape2=shape2)
      }
      if (type == 2) {
        z <- CondPower(n, i, k, k.crit.fin, p1)
      }
      if (type == 3) {
        pmle <- k/i
        z <-   z <- CondPower(n, i, k, k.crit.fin, pmle)
      }

      if (z >= crit) break
      # k is high enough to warrant continuation of the trial
      # the last value of k is then critical for early stopping
      v.boundary[i] <- k

    }  # end loop over observed successes k
  }  # end loop over patient numbers i

  # Final analysis:
  v.boundary[n] <- k.crit.fin - 1

  return(v.boundary)

}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

CritBoundary1=function (n, crit, k.crit.fin, EF, shape1F=1, shape2F=1,shape1E=NULL, shape2E=NULL) {

  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F

  # CritBoundary1

  # Vector of critical boundaries corresponding to interim analyses based on predictive
  # power and either considering efficacy or futility. The last element of the vector
  # corresponds to the test-based final analysis.
  # Futility: With an analysis performed at i patients, an observed number k of successes
  # leads to a negative decision (suspension/stopping of the trial for futility or a non-
  # signifcant result in the final analysis) iff k <= v.boundary[i].
  # Efficacy: to a negative decision (suspension/stopping of the trial for efficacy, a
  # significant result in the final analysis) iff k >= v.boundary[i].
  # The critical boundaries calculated for early stopping at interim analyses are based
  # on the critical boundaries (k.crit.fin) for success of the trial in the final analysis.

  # Args:
  # n:            sample size in the final analysis
  # crit:         critical level (%) for probability of a negative/positive judgement
  #                 on the treatment based on predictive power
  # k.crit.fin    critical value for significant results used in the final analysis
  # EF:           type of stopping criterion: 'E'=efficacy, 'F'=futility
  # shape1:       first parameter of the beta prior
  # shape2;       second parameter of the beta prior

  # Returns:
  # Vector v.boundary of of critical boundaries for an analysis at i=1,...,n patients.

  # Function invoked: PredPower
  #--------------------------------------------------------------------------------------------

  if (EF == 'F') {
    v.boundary <- rep(-1, n)
  }
  else {
    v.boundary <- rep(n + 1, n)
  }
  # Initially: no stopping for efficacy or futility

  # Loop over patient numbers 1,...,n-1 (corresponding to number of patients avalaible for
  # interim analyses):

  for (i in 1:(n - 1)) {

    if (EF == 'F') {

      # Loop over the number of successes observed among the i patients
      for (k in 0:i) {

        z <- PredPower(n, i, k, k.crit.fin, shape1F, shape2F)

        if (z >= crit) break
        # k is high enough to warrant continuation of the trial
        # the last value of k is then critical for early stopping
        v.boundary[i] <- k

      }  # end loop over observed successes k

    }  # end if EF=='F'

    if (EF == 'E') {

      # Loop over the number of successes observed among the i patients
      for (k in i:0) {

        z <- PredPower(n, i, k, k.crit.fin, shape1E, shape2E)

        if (z <= crit) break
        # k is low enough to warrant continuation of the trial
        # the last value of k is then critical for early stopping
        v.boundary[i] <- k

      }  # end loop over observed successes k

    }  # end if EF=='E'

  }  # end loop over patient numbers i

  # Final analysis:
  if (EF == 'F') {
    v.boundary[n] <- k.crit.fin - 1
  }
  else {
    v.boundary[n] <- k.crit.fin
  }

  return(v.boundary)

}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
CritBoundaryE5=function (n, pE, critE, shape1=1,shape2=1) {

  # CritBoundaryE

  # Vector of critical boundaries corresponding to interim analyses and final analysis,
  # based on the posterior distribution (=beta distribution) for the response probability r.

  # With an analysis performed at i patients, an observed number k of responses
  # leads to a positive decision (early stopping of the trial for efficacy, a positive result
  # in the final analysis) iff k >= v.boundary[i]. More specifically, k is the smallest
  # integer such that P(r >= pE) >= crit.

  # Args:
  # n:           total number of study patients
  # pE:          response rate corresponding to the futility bound
  # critE        critical level of the posterior probability that r >= pE used for
  #                 early termination

  # Returns:
  # Vector v.boundary of of critical boundaries for an analysis at i=1,...,n patients.

  #-----------------------------------------------------------------------------------------

  v.boundary <- rep(n+1,n)
  # Initially: no stopping for efficacy

  # Loop over patient numbers 1,...,n:
  for (i in 1:n) {
    for (k in 0:i) {
      # k = number of observed responses among i patients

      z <- pbeta(q = pE, shape1 = shape1 + k, shape2 = shape2 + i - k, lower.tail = FALSE)
      # =Posterior probability (calculated from k responses observed among i
      # patients) that r >= pE

      if (z >= critE) {
        v.boundary[i] <- k
        # k responses are sufficient to establish efficacy
        break
      }

      v.boundary[i] <- k
    }
    if(z<critE) v.boundary[i]=i+1
  }

  return(v.boundary)

}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
CritBoundaryF5=function (n, pF, critF, shape1=1,shape2=1) {

  # CritBoundaryF

  # Vector of critical boundaries corresponding to interim analyses and final analysis,
  # based on the posterior distribution (=beta distribution) for the response probability r.

  # With an analysis performed at i patients, an observed number k of responses
  # leads to a negative decision (suspension/early stopping of the trial, a null result
  # in the final analysis) iff k <= v.boundary[i]. More specifically, k is the highest
  # integer such that P(r >= pF) < crit.

  # Args:
  # n:           total number of study patients
  # pF:          response rate corresponding to the futility bound
  # critF:       critical level of the posterior probability that r >= pF used for
  #                 early termination

  # Returns:
  # Vector v.boundary of of critical boundaries for an analysis at i=1,...,n patients.

  #-----------------------------------------------------------------------------------------

  v.boundary <- rep(-1, n)
  # Initially: no stopping for futility

  # Loop over patient numbers 1,...,n:
  for (i in 1:n) {
    for (k in 0:i) {
      # k = number of observed responses among i patients

      z <- pbeta(q = pF, shape1 = shape1 + k, shape2 = shape2 + i - k, lower.tail = FALSE)
      # =Posterior probability (calculated from k responses observed among i
      # patients) that r >= pF

      if (z >= critF) break
      # k is high enough to warrant the continuation of the trial;

      # the last value of k is then the critical one.
      v.boundary[i] <- k
    }
  }

  return(v.boundary)

}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
CritBoundaryPost=function (n, p.hyp, crit,shape1=1,shape2=1) {

  # CritBoundaryPost

  # Vector of critical boundaries corresponding to interim analyses and final analysis,
  # based on the posterior distribution (=beta distribution) for the success probability
  # (rate) psuccess.
  # With an analysis performed at i patients, an observed number k of successes
  # leads to a negative decision (suspension/early stopping of the trial; a null result
  # in the final analysis) iff k <= v.boundary[i]. More specifically, k is the highest
  # integer such that P(psuccess > (=) p.hyp) < crit (if p.hyp corresponds to H0,
  # crit must be high, e.g., = 0.7 or even = 0.95; if p.hyp corresponds to H1, crit
  # must be low, e.g. = 0.1).

  # Args:
  # n:           total number of study patients
  # p.hyp:       success rate (%) corresponding to H0 or H1 (i.e., argument p0 or p1 from
  #                function Futility)
  # crit:        critical level of the posterior probability that psuccess > p0 or psuccess
  #                >= p1 used for early termination

  # Returns:
  # Vector v.boundary of of critical boundaries for an analysis at i=1,...,n patients.

  #-----------------------------------------------------------------------------------------

  v.boundary <- rep(-1, n)
  # Initially: no stopping for futility

  # Loop over patient numbers 1,...,n:
  for (i in 1:n) {
    for (k in 0:i) {
      # k = number of observed successes among i patients

      z <- pbeta(q = p.hyp, shape1 = shape1 + k, shape2 = shape2 + i - k, lower.tail = FALSE)
      # =Posterior probability (calculated from k successes observed among i
      # patients) that psuccess >(=) p.hyp

      if (z >= crit) break
      # k is high enough to warrant the continuation of the trial;

      # the last value of k is then the critical one.
      v.boundary[i] <- k
    }
  }

  return(v.boundary)

}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
CritBoundaryRates=function (n, crit) {

  # CritBoundaryRates

  # Vector of critical boundaries corresponding to interim analyses based on the
  # critical rate (asssumed to be identical for all interim analyses).
  # With an analysis performed at i patients, an observed number k of successes
  # leads to a negative decision (suspension/early stopping of the trial)
  # if k/i <= v.boundary[i].

  # Args:
  # n:           total number of study patients
  # crit:        critical success rate. Futility is declared if the estimated rate is
  #                < crit

  # Returns:
  # Vector v.boundary of of critical boundaries for an analysis at i=1,...,n patients.

  #-----------------------------------------------------------------------------------------

  v.boundary <- rep(-1, n)
  # Initially: no stopping for futility

  # Loop over patient numbers 1,...,n:
  for (i in 1:n) {
    for (k in 0:i) {
      # k = number of observed successes among i patients

      if (k/i >= crit) break
      # k is high enough to warrant the continuation of the trial;
      # the last value of k is then the critical one.

      v.boundary[i] <- k
    }
  }

  return(v.boundary)

}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
CritBoundaryTest=function (n, p, crit, type) {

  # CritBoundaryTest

  # Vector of critical boundaries corresponding to interim analyses based on p values
  # calculated under p.

  # type = 8:  With an analysis performed at i patients, an observed number k of
  # successes leads to a negative decision (suspension/early stopping of the trial)
  # if p-value=p(k):=P(number of successes >=k| success probability=p0) <= v.boundary[i],
  # which is equivalent to p(k) >= pcrit.

  # type=9: Stop/suspend if p'(k):=P(number of successes <=k| success probability=p1)
  # <= v.boundary[i], which is equivalent to p'(k) < pcrit.

  # Args:
  # n:           total number of study patients
  # p:           probability of success used for evaluating the criterion (p=p0 or =p1)
  # crit:        critical alpha level. Futility is declared if the p-value calculated
  #                under success probability p is >= crit (type 8) or < crit (type 9)
  # type:        type of futility analysis (type 8 or 9)

  # Returns:
  # Vector v.boundary of of critical boundaries for an analysis at i=1,...,n-1 patients.

  #-----------------------------------------------------------------------------------------

  v.boundary <- rep(-1, n - 1)
  # Initially: no stopping for futility

  # Loop over patient numbers 1,...,n-1:
  for (i in 1:(n - 1)) {
    for (k in 0:i) {
      # k = number of observed successes among i patients

      if (type == 8) {
        pvalue <- pbinom(q = k - 1, size = i, prob = p, lower.tail = FALSE)
        #print('i,k')
        #print(c(i,k))
        #print(c('pvalue',pvalue))
        if (pvalue < crit) break
        # k is high enough to warrant the continuation of the trial;
        # the last value of k is then the critical one.
        # Note that (exceptionally and due to particular form of futility
        # criterion of type 8), the "<" sign must be used in the if-statement).

        v.boundary[i] <- k
      }
      if (type == 9) {
        pvalue <- pbinom(q = k, size = i, prob = p, lower.tail = TRUE)
        if (pvalue >= crit) break
        # k is high enough to warrant the continuation of the trial;
        # the last value of k is then the critical one.

        v.boundary[i] <- k
      }
    }
  }

  return(v.boundary)

}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
critEF=function (n,  vn.int,  crit,  EF="F", shape1=1, shape2=1, pE=NULL,pF=NULL) {

  # Vector of critical boundaries corresponding to interim analyses and final analysis,
  # based on the posterior distribution (=beta distribution) for the response probability r
  # for either futility or efficacy evaluation.

  # Futility: With an analysis performed at i patients, an observed number k of responses
  # leads to a negative decision (suspension/early stopping of the trial, a null result
  # in the final analysis) iff k <= v.boundary[i]. More specifically, k is the highest
  # integer such that P(r >= pF) < crit.

  # Efficacy: With an analysis performed at i patients, an observed number k of responses
  # leads to a positive decision (early stopping of the trial for efficacy, a positive result
  # in the final analysis) iff k >= v.boundary[i]. More specifically, k is the smallest
  # integer such that P(r >= pE) >= crit.

  # Args:
  # n:           total number of study patients
  # vn.int:      vector of sample sizes at the interim analyses
  # crit        critical level of the posterior probability that r > pE or r > pF used for
  #                 early termination
  # pF:          response rate corresponding to the futility bound
  # pE:          response rate corresponding to the efficacy bound
  # EF:         decides whether calculations are for futility (="F") or efficacy(="E)
  # shape1:     shape1 parameter for prior distribution
  # shape2:     shape2 parameter for prior distribution


  # Returns:
  # Vector v.boundary of of critical boundaries for an analysis at i=1,...,n patients.


  #--------------------------------------------------------------------------------------------
  # Input check:
  if (max(vn.int) >= n) stop('Interim analyses must be performed at patient numbers < n.')
  if (is.unsorted(vn.int) == TRUE) stop('vn.int must be sorted in ascending order')

  # Preliminary commands:
  n.int <- length(vn.int)
  n.an <- n.int + 1
  vn.an <- c(vn.int, n)





  #--------------------------------------------------------------------------------------------
  # 2. Vectors of critical boundaries corresponding to interim analyses at each possible
  #    patient number and final analysis:

  if (EF == 'F') {
    if( is.null(pF)) stop("for type 5 (posterior probabilities), stopping for futility, pF must be specified.")
    v.boundary <- CritBoundaryF5(n, pF, crit,shape1=shape1,shape2=shape2)
  }
  else {
    if( is.null(pE)) stop("for type 5 (posterior probabilities), stopping for efficacy, pE must be specified.")
    v.boundary <- CritBoundaryE5(n, pE, crit,shape1=shape1,shape2=shape2)

  }


  v.crit <- v.boundary[vn.an]
  v.crit
}


