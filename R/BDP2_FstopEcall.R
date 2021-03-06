BDP2_FstopEcall=function(n, vn.int, p, pF=NULL,cF,pE=NULL,  cE,
                                       shape1F, shape2F,shape1E=NULL, shape2E=NULL){
  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F

  v.critE <- critEF(n=n,  vn.int=vn.int, crit=cE, pE=pE,pF=pF,  EF="E", shape1=shape1E, shape2=shape2E)
  v.critF <- critEF(n=n,  vn.int=vn.int, crit=cF, pE=pE,pF=pF,  EF="F", shape1=shape1F, shape2=shape2F)

  my.pFstopEcall <- pFstopEcall(p,c(vn.int,n),v.critE,v.critF)

  print(my.pFstopEcall)

  probstop_prop <- c(as.numeric(my.pFstopEcall[,6]))

  #calculations for the frequency distribution of #patients
    cum.probstop_prop <- cumsum(probstop_prop)
    probabilities <- c(probstop_prop[-length(probstop_prop)],1-cum.probstop_prop[length(probstop_prop)-1])
  #vector probabilites gives the frequency distribution
    expected.n <- sum(c(vn.int,n)*probabilities)

  cat("\nExpected number of patients enrolled in the trial:", round(expected.n, digits = 2))
  Pinconcl <- 1 - my.pFstopEcall[length(c(vn.int,n)),4] -my.pFstopEcall[length(c(vn.int,n)),7]
  cat("\nEstimated probability of inconclusive results                   :", round(Pinconcl, digits = 4))
  invisible(list(summary=my.pFstopEcall,expectedNumber=expected.n,probInconclusive=Pinconcl))
}


