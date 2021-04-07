BDP2_simulateEF_PredPow=function (n, vn.int, p, p0, p1, cF, cE, alpha, nsim, eff.stop = TRUE,
              shape1F=1, shape2F=1,shape1E=NULL, shape2E=NULL) {

  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F

# EF1

# Determines the operating characteristics of a single-arm trial with a binary
# endpoint (response, success) and interim efficacy and futility analyses.

# Assumptions:
# Endpoint (success/no success) data available for all study patients.
# Prior distribution = Beta(shape1, shape2); default=uniform, i.e., shape1=
# shape2=1.
# One-sided testing in the final analysis.

# Some methodological details:

# Given the results of the interim analysis, the predictive power at the final analysis
# (n patients, critical number of successes k.crit) is P(X >= k.crit - k.int), where
# X  follows a beta-binomial distribution with parameters n'= n - n.int, a = k.int + shape1,
# and b = n.int - k.int + shape2.

# Efficacy is declared if the predictive power is >= cE (cE must be high, e.g. 70%).
# Futility is declared if the predictive power is < cF (cF must be small, e.g. 10%).
# cE, cF translate into futility/efficacy boundaries (maximum number of
# responses leading to early termination for futility/ minimum number of responses
# leading to declaring of, or early termination for, efficacy).

# Args:
# n:           sample size at the final analysis
# p0:          response rate (not in %!) corresponding to H0
# p1:          response rate (not in %!) corresponding to H1 (p1 > p0)
# p:           true (assumed) response rate (%) used for simulating the trial
# vn.int:      vector of sample sizes at the interim analyses (the vector may be equal
#                to 1:(n-1) = continuous monitoring of futility)
# cE:       critical level of posterior probabilities (not in %!) used for declaring efficacy
# cF:       critical level of posterior probabilities (not in %!) used for declaring futility
# nsim:        number of simulation runs
# eff.stop:    TRUE (yes=default), i.e. the study ends if the efficacy criterion is
#                reached at an interim analysis. Different input: no stop for
#                efficacy; in this case the program merely calculates the probability
#                that the efficicacy criterion is satisfied (possibly triggering a
#                notification of the DMC and the start of the planning of a subsequent
#                trial).
# shape1:       first parameter of the Beta prior (default = 1)
# shape2:       second parameter of the Beta prior (default = 1)

# Output:
# Prints the results of the analysis.

# Functions invoked:
# CritBoundary1, SimEF

#--------------------------------------------------------------------------------------------
# Input check:
  if (p0 > p1) stop('p1 must be > p0.')
  if (max(vn.int) >= n) stop('Interim analyses must be performed at patient numbers < n.')
  if (is.unsorted(vn.int) == TRUE) stop('vn.int must be sorted in ascending order')

# Preliminary commands:

  n.int <- length(vn.int)
  n.an <- n.int + 1
  vn.an <- c(vn.int, n)


# Critical values for testing used in the final analysis (one-sided binomial test):

  k.crit.fin <- qbinom(size = n, prob = p0, p = alpha, lower.tail = FALSE) + 1

# Power of the final test:

  power.fin <- pbinom(size = n, q = k.crit.fin - 1, prob = p1, lower.tail = FALSE)
#  power.fin <- round(100*power.fin, digits=1)

  cat('\n')
#  cat('Exact power of the test used in the final analysis (%):  ', power.fin, '\n\n')
  cat('Exact power of the test used in the final analysis:  ', power.fin, '\n\n')

  if (power.fin <= cF) {
    stop('        The critical predictive power used for early stopping for futility
          must be lower than the power of the final test.')
  }
  if (power.fin > cE) {
    stop('        The critical predictive power used for declaring efficacy in interim
          analyses must be at least as high as the power of the final test.')
  }

# Vectors of critical boundaries corresponding to interim analyses at each possible
# patient number and final analysis:

  v.boundaryE <- CritBoundary1 (n, cE, k.crit.fin, 'E', shape1E, shape2E)
  v.boundaryF <- CritBoundary1 (n, cF, k.crit.fin, 'F', shape1F, shape2F)

  v.critE <- v.boundaryE[vn.an]
  v.critF <- v.boundaryF[vn.an]

# Simulations determining the OC
# (Probabilities of delaring efficacy/futility for the assumed response rate p):
  SimResEF <- SimEF (n, p, vn.int, v.boundaryE, v.boundaryF, nsim, eff.stop)

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

# Output:
  cat('#######################################################################\n\n')
  cat('Simulation results (Operating Characteristics of the design):\n\n')
  cat('\n')
  cat('Critical boundaries for declaring efficacy/futility and corresponding\n')
  cat('probabilites.\n\n')
  cat('Efficacy is declared if the number of responses at an interim analysis is\n')
  cat('>= the critical boundary Crit.BoundaryE.\n')
  cat('Futility is declared if the number of responses at an interim analysis is\n')
  cat('<= the critical boundary Crit.BoundaryF.\n')
  cat('The last row corresponds to the final analysis.\n\n')
  print(SimResDF)
  cat('\n')
  cat('Expected number of patients enrolled in the trial   :', mean.npat, '\n\n')
}
