PredPower=function (n, n.int, k.int, k.crit, shape1, shape2) {

# PredPower

# Calculates the predictive power at an interim analysis.

# Assumptions: 
# Beta-binomial model. 
# Binary variable success-failure. One-sided testing.
# Prior distribution: Beta(shape1, shape2).
 
# The posterior distribution at interim  analysis with n.int patients and k.int 
# successes is then equal to Beta(k.int + shape1, n + shape2 - k.int) and - given
# the results of the interim analysis - the predictive power for a significant
# result in the final analysis (n patients, critical number of successes k.crit) 
# is P(X >= k.crit - k.int), where X  follows a beta-binomial distribution 
# with parameters n'= n - n.int, a = k.int + shape1, and b = n.int - k.int + shape2.

# Args:
# n:           sample size at the final analysis
# n.int:       sample size at the interim analsis
# k.int:       number of successes observed up until the interim analysis
# k.crit:      critical number of successes at the final analysis (lower tail)   
# shape1;      first parameter of the Beta prior
# shape2:      second parameter of the Beta prior

# Returns:
# The predictive power.

# Function invoked: BetaBinom
#-----------------------------------------------------------------------------
# Predictive power, upper tail:

  pred.power <- BetaBinom(n - n.int, k.crit - k.int, k.int + shape1, 
                        n.int + shape2 - k.int)[2]

return(pred.power)

}
