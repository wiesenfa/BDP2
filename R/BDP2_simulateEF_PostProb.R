BDP2_simulateEF_PostProb=function (n, vn.int, p, pF, cF, pE, cE,nsim, eff.stop = TRUE,
                    shape1F =1, shape2F = 1,shape1E=NULL, shape2E=NULL) {

  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F

# EF5

# Determines the operating characteristics of a single-arm trial with a binary
# endpoint (response, success) and interim efficacy and futility analyses.
# Declaration of efficacy and futility (including possibly early stopping)
# is based on the posterior probability that the true response rate
# is at least pE , pF respectively.

# Assumptions:
# Endpoint (response/no response) data available for all study patients.
# Beta-binomial model. Prior distribution = Beta(shape1, shape2). Default=uniform.

# The posterior distribution at interim analysis with n.int  patients and k.int
# successes is Beta(k.int + shape1, n.int + shape2 - k.int)
# Efficacy is declared if the posterior probability P(true response rate >= pE)
# is >= cE (cE must be high, e.g. 70%).
# Futility is declared if the posterior probability P(true success rate >= pF)
# is < cF (cF must be small, e.g. 10%).
# cE, cF translate into futility/efficacy boundaries (maximum number of
# responses leading to early termination for futility/ minimum number of responses
# leading to declaring of, or early termination for, efficacy).

# Args:
# n:           sample size at the final analysis
# pE:          response rate (not in %!) used for the efficacy criterion
# pF:          response rate (not in %!) used for the futility criterion (may be identical to pE)
# p:           true (assumed) response rate (not in %!) used for simulating the trial
# vn.int:      vector of sample sizes at the interim analyses (the vector may be equal
#                to 1:(n-1) = continuous monitoring of futility)
# cE:       critical level of posterior probabilities used for declaring efficacy (not in %!)
# cF:       critical level of posterior probabilities used for declaring futility (not in %!)
# nsim:        number of simulation runs
# eff.stop:    TRUE (yes=default), i.e. the study ends if the efficacy criterion is
#                reached at an interim analysis. Different input: no stop for
#                efficacy; in this case the program merely calculates the probability
#                that the efficicacy criterion is satisfied (possibly triggering a
#                notification of the DMC and the start of the planning of a subsequent
#                trial).
# shape1:      first parameter of the Beta prior (default = 1)
# shape2:      second parameter of the Beta prior (default = 1)

# Output:
# Prints the results of the analysis.

# Functions invoked:
# CritBoundaryE5, CritBoundaryF5, SimEF

#--------------------------------------------------------------------------------------------
# Input check:
  if (max(vn.int) >= n) stop('Interim analyses must be performed at patient numbers < n.')
  if (is.unsorted(vn.int) == TRUE) stop('vn.int must be sorted in ascending order')

# Preliminary commands:

  n.int <- length(vn.int)
  n.an <- n.int + 1
  vn.an <- c(vn.int, n)


# Vectors of critical boundaries corresponding to interim analyses at each possible
# patient number and final analysis:
  v.boundaryE <- CritBoundaryE5(n, pE, cE, shape1E, shape2E)
  v.boundaryF <- CritBoundaryF5(n, pF, cF, shape1F, shape2F)
  v.critE <- v.boundaryE[vn.an]
  v.critF <- v.boundaryF[vn.an]

# Simulations determining the OC
# (Probabilities of delaring efficacy/futility for the assumed response rate p):
  SimResEF <- SimEF(n, p, vn.int, v.boundaryE, v.boundaryF, nsim, eff.stop)

# Data frames for output:
  if (eff.stop) {
    SimResDF <- data.frame(Int.Analysis = 1:n.an, Pat.number = vn.an, Crit.boundaryE = v.critE,
                  P.effic = SimResEF[[1]], P.effic.cum = SimResEF[[2]], Crit.boundaryF = v.critF,
                  P.futil = SimResEF[[3]], P.futil.cum = SimResEF[[4]])
  }
  else {
    SimResDF <- data.frame(Int.Analysis = 1:n.an, Pat.number = vn.an, Crit.boundaryE = v.critE,
                  P.effic = SimResEF[[1]], Crit.boundaryF = v.critF,
                  P.futil = SimResEF[[2]], P.futil.cum = SimResEF[[3]])
  }

# Expected number of patients:
   if (eff.stop) {
    mean.npat <- SimResEF[[5]]
  }
  else {
   mean.npat <- SimResEF[[4]]
  }

# Estimated probability of inconclusive results (%):
   if (eff.stop) {
    prob.inconclus <- SimResEF[[6]]
  }
  else {
    prob.inconclus <- SimResEF[[5]]
  }


# Output:
  cat('#######################################################################\n\n')
  cat('Simulation results (Operating Characteristics of the design):\n\n')
  cat('\n')
  cat('Critical boundaries for declaring efficacy/futility and corresponding\n')
  cat('probabilities.\n\n')
  cat('Efficacy is declared if the number of responses at an interim analysis is \n')
  cat('>= the critical boundary Crit.BoundaryE.\n')
  cat('Futility is declared if the number of responses at an interim analysis is \n')
  cat('<= the critical boundary Crit.BoundaryF.\n')
  cat('The last row corresponds to the final analysis\n\n')
  print(SimResDF)
  cat('\n')
  cat('Expected number of patients enrolled in the trial   :', mean.npat, '\n\n')
  cat('Estimated probability of inconclusive results   :', prob.inconclus, '\n\n')

}
