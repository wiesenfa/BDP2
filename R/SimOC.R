SimOC=function (n, p, vn.int, v.boundary, nsim) {

# SimOC

# Simulation program determining the probability of early stopping (or, 
# of a null result obaitained in the final analysis);
# In case of statistical testing used in the final analysis: type 1 error
# and power are also determined (accounting for interim analyses). 

# Args:
# n:           sample size at the final analysis
# p:           true (assumed) success rate (not in %!) used for simulating the trial
# vn.int:      vector of sample sizes at the interim analyses 
# v.boundary:  vector of critical boundaries for curtailment . 
#                 The i-th component (0<i<=n) contains the critical boundary 
#                 (critical number of observed successes for early stopping/a null
#                 result in the final analysis (respectively) when evaluated at 
#                 the ith patient. 
#                 Stop/declare the result non-significant if at most v.boundary[i] 
#                 successes have been observed up until the i-th patient.
# nsim:        number of simulation runs

# Returns:
# A list containing three elements: two vectors probstop, probstop.sum giving 
# the probabilities/cumulative probabilites of stopping (null result in the final 
# analysis, respectively) for each analysis; the third element is the mean patient 
# number enrolled.

#-----------------------------------------------------------------------------------

# Preliminary commands:
  vn.an <- c(vn.int, n)
  n.an <- length(vn.an)
  nstop <- rep(0, n.an)
    # nstop is a vector of length n.an whose m-th component contains the number
    # of simulation runs in which the trial is stopped at the m-th interim
    # analysis, or (if no early stopping occurs) a null result is obtained in the 
    # final analysis. 
  v.crit <- v.boundary[vn.an]
  vsim.npat <- rep(n, nsim)
    # vector of length nsim whose m-th component contains the number of patients 
    # actually enrolled in the m-th simulated trial. Initially, n patient are
    # assumed to be enrolled


# Simulation:

  for (i in 1:nsim) {

    stop <- 0

    # Binomial random draw for endpoint success for all patients: 
      success <- rbinom(size = 1, prob = p, n = n)

    # Cumulative number of successes:
      success.cum <- cumsum(success)
    
    # Loop over analyses:

      for (j in 1:n.an) {
      
        npat <- vn.an[j]

        # Stop (or null result at final analysis) if number of successes 
        # <= critical value:             
          if (success.cum[npat] <= (v.crit[j] + 0.5))  {
            stop <- 1  
            vsim.npat[i] <- npat
            nstop[j] <- nstop[j] + 1
          } 

        if (stop == 1) break
          # No further analyses are carried out 
        
       } # End loop over analyses
 
  }  # End loop over the simulation runs

# Percentage of studies stopping at the j-th analysis (j=1,..,n.an):
  probstop <- round(nstop/nsim, digits = 4) 

# Vector of cumulative probabilites of stopping:
  probstop.cum <- round(cumsum(probstop), digits = 4)

# Mean patient number enrolled:
  mean.npat <- round(mean(vsim.npat), digits = 1)
  
# Output:
  return(list(probstop, probstop.cum, mean.npat))

}
