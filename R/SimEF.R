SimEF=function (n, p, vn.int, v.boundaryE, v.boundaryF, nsim, eff.stop) {

# SimEF

# Simulation program determining the probability of early stopping for futility
# or declaring/stopping for efficacy (or obtaining the corresponding results in
# the final analysis);

# Args:
# n:            sample size at the final analysis
# p:            true (assumed) response rate (not in %!) used for simulating the trial
# vn.int:       vector of sample sizes at the interim analyses
# v.boundaryE:  vector of critical boundaries for declaring efficacy.
#                 The i-th component (0<i<=n) contains the critical boundary
#                 (critical number of observed responses for early stopping/a
#                 positive result in the final analysis (respectively) when
#                 evaluated at the i-th patient.
#                 Stop/declare efficacy if at least v.boundary[i] many
#                 responses have been observed up until the i-th patient.
# v.boundaryF:  vector of critical boundaries for declaring futility.
#                 The i-th component (0<i<=n) contains the critical boundary
#                 (critical number of observed responses for early stopping/a
#                 negative result in the final analysis (respectively) when
#                 evaluated at the i-th patient.
#                 Stop/declare futility if at most v.boundary[i] many
#                 responses have been observed up until the i-th patient.
# nsim:         number of simulation runs
# eff.stop:     'y': the study ends if the efficacy criterion is
#                reached at an interim analysis. Values different form 'y': no
#                stop for efficacy.

# Returns:
# A list containing six elements: vectors probstopE, probstop.sumE giving
# the probabilities/cumulative probabilites of stopping for efficacy (positive result
# in the final analysis, respectively) for each analysis;
# Analogously: vectors probstopF, probstop.sumF (for futility);
# the fifth element is the mean number of patients  enrolled.
# the sixth element is the proportion of inconclusive results.

#-----------------------------------------------------------------------------------

# Preliminary commands:
  vn.an <- c(vn.int, n)
  n.an <- length(vn.an)
  nposE <- nstopE <- rep(0, n.an)
    # nposE and nstopE are vectors of length n.an whose m-th component contains the number
    # of simulation runs in which the efficalcy criterion is triggered (and the accrual stopped,
    # respectively) at the m-th analysis.
  nstopF <- rep(0, n.an)
    # nstopF is a vector of length n.an whose m-th component contains the number
    # of simulation runs in which the trial is stopped for futility at the m-th
    # interim analysis, or (if no early stopping occurs) a negative result is declared
    # in the final analysis.
  n.inconclus <- 0
  v.critE <- v.boundaryE[vn.an]
  v.critF <- v.boundaryF[vn.an]
  vsim.npat <- rep(n, nsim)
    # vector of length nsim whose m-th component contains the number of patients
    # actually enrolled in the m-th simulated trial. Initially, n patient are
    # assumed to be enrolled

# Simulation:

  for (i in 1:nsim) {

    stopF <- stopE <- 0
     # Variables indicating early stopping for efficacy or futility

    posE.last <- 0
     # This variable indicates a positive result obtained in the final analysis

    # Binomial random draw for endpoint response for all patients:
      response <- rbinom(size = 1, prob = p, n = n)

    # Cumulative number of responses:
      response.cum <- cumsum(response)

    # Loop over analyses:

      for (j in 1:n.an) {

        npat <- vn.an[j]

       # Positive result if number of responses = critical value for efficacy:
           if (response.cum[npat] >= (v.critE[j] - 0.5))  {
             nposE[j] <- nposE[j] + 1
             if (j == n.an) posE.last <- 1

             if (eff.stop == 'y') {
               vsim.npat[i] <- npat
               stopE <- 1
               nstopE[j] <- nstopE[j] + 1
                  # if stop for efficacy is planned then accrual stops here
             }
           }

       # Stop (or null result at final analysis) if number of successes
       # <= critical value for futility:
           if (response.cum[npat] <= (v.critF[j] + 0.5))  {
             stopF <- 1
             vsim.npat[i] <- npat
             nstopF[j] <- nstopF[j] + 1
           }

       if ((j == n.an) && (posE.last + stopF == 0)) n.inconclus <- n.inconclus + 1
             # No declaration of efficacy or futility up until the final analysis:

       if (((stopE == 1) && (eff.stop == 'y')) || (stopF == 1)) break
           # No further analyses are carried out

       } # End loop over analyses

  }  # End loop over the simulation runs

# Proportion of studies with efficacy found at the j-th analysis (j=1,..,n.an):
  prob.eff <- round(nposE/nsim, digits = 4)

# Proportion of studies stopping at the j-th analysis (j=1,..,n.an):
  prob.effstop <- round(nstopE/nsim, digits = 4)
  prob.futilstop <- round(nstopF/nsim, digits = 4)

# Vector of cumulative probabilites of stopping:
  prob.effstop.cum <- round(cumsum(prob.effstop), digits = 4)
  prob.futilstop.cum <- round(cumsum(prob.futilstop), digits = 4)

# Mean patient number enrolled:
  mean.npat <- round(mean(vsim.npat), digits = 1)

# Estimated probability of inconclusive results:
  prob.inconclus <- round(n.inconclus/nsim, digits = 4)

# Output:
  if (eff.stop == 'y') {
    return(list(prob.effstop, prob.effstop.cum, prob.futilstop, prob.futilstop.cum, mean.npat, prob.inconclus))
  }
  else {
    return(list(prob.eff, prob.futilstop, prob.futilstop.cum,  mean.npat, prob.inconclus))
  }

}
