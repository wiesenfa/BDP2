expectedNumberFstopEstop=function(n, p, pF, vn.int, cF, pE=NULL, cE,
                                  shape1F, shape2F,shape1E=NULL, shape2E=NULL){
  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F
  v.critE <- critEF(n=n,  vn.int=vn.int, crit=cE, pE=pE,pF=pF,  EF="E", shape1=shape1E, shape2=shape2E)
  v.critF <- critEF(n=n,  vn.int=vn.int, crit=cF, pE=pE,pF=pF,  EF="F", shape1=shape1F, shape2=shape2F)

  my.pFstopEstop <- pFstopEstop(p,c(vn.int,n),v.critE,v.critF)
  
  print(my.pFstopEstop)
  
  probstop_prop <- c(as.numeric(my.pFstopEstop[,4])+my.pFstopEstop[,7])
  
  #calculations for the frequency distribution of #patients
  cum.probstop_prop <- cumsum(probstop_prop)
  probabilities <- c(probstop_prop[-length(probstop_prop)],1-cum.probstop_prop[length(probstop_prop)-1])
  
  
  expected.n <- sum(c(vn.int,n)*probabilities)
  invisible(expected.n)
}
