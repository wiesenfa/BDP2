BetaBinom=function (n, x, a, b) {

# BetaBinom

# Calculates and returns P(X<=x) and P(X>=x), where X follows a beta-binomial 
# distribution with parameters n, a, b.
#----------------------------------------------------------------------------
v.low <- 0:x
v.up <- x:n

P.low <- sum(choose(n, v.low)*beta(a + v.low, b + n - v.low)/beta(a, b))
P.up <- sum(choose(n, v.up)*beta(a + v.up, b + n - v.up)/beta(a, b))

if (x > n) {
  P.low <- 1
  P.up <- 0
}

return(c(P.low, P.up))
}
