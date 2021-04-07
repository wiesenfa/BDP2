
BDP2_Fstop=function(n, vn.int, p, pF, cF, shape1=1,shape2=1){

  v.crit <-crit(n=n, pThreshold=pF,  vn.int=vn.int, crit=cF, shape1 = shape1, shape2=shape2)
  my.pstop <- pFstop(p,c(vn.int,n),v.crit)
  my.expectedNumber <- expectedNumberFstop(n, p,pF, vn.int, cF, shape1,shape2)

  print(my.pstop)
  cat('\n')
  cat("Expected number of patients enrolled in the trial:", round(my.expectedNumber, digits = 2))
  invisible(list(summary=my.pstop,expectedNumber=my.expectedNumber))
}


