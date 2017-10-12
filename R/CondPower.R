CondPower=function (n, n.int, k.int, k.crit, p) {

# CondPower

# Calculates the conditional power for a binominal endpoint (success) at an 
# interim analysis. Assumption: True success rate = p. 

# Args:
# n:           sample size at the final analysis
# n.int:       sample size at the interim analsis
# k.int:       number of events observed up until the interim analysis
# k.crit:      critical number of events at the final analysis (lower tail)   
# p:           success rate 

# Returns:
# The conditional power.

#---------------------------------------------------------------------------------
# Conditional power, upper tail:

  cond.power <- pbinom(size = n - n.int, q = (k.crit - k.int) - 1, prob = p, 
                          lower.tail = FALSE)

return(cond.power)

}
