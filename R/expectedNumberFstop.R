expectedNumberFstop=function(n, p, pF, vn.int, cF, shape1=1,shape2=1){
  v.crit <-crit(n=n, pThreshold=pF, vn.int=vn.int, crit=cF, shape1 = shape1, shape2 = shape2)
  probstop_prop <- c(as.numeric(pFstop(p,c(vn.int,n),v.crit)[,4]))

  #calculations for the frequency distribution of #patients

  cum.probstop_prop <- cumsum(probstop_prop)
  probabilities <- c(probstop_prop[-length(probstop_prop)],1-cum.probstop_prop[length(probstop_prop)-1])
 #vector probabilites gives the frequency distribution

  expected.n <- sum(c(vn.int,n)*probabilities)
  invisible(expected.n)
}