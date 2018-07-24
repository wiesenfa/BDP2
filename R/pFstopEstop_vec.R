pFstopEstop_vec=function(n, p, pE,pF, vn.int, cEvec, cF, shape1F,shape2F,shape1E=NULL,shape2E=NULL,printit=FALSE){

  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F
  
  out=lapply(cEvec, function(cE) {
    v.critE <- critEF(n=n,  vn.int=vn.int, crit=cE, pE=pE,pF=pF,  EF="E", 
                      shape1=shape1E,shape2=shape2E)
    v.critF <- critEF(n=n,  vn.int=vn.int, crit=cF, pE=pE,pF=pF,  EF="F", 
                      shape1=shape1F,shape2=shape2F)
    res <- pFstopEstop(p,c(vn.int,n),v.critE,v.critF)
    res
  }
  )
  names(out)=cEvec

#  summ=sapply(out, function(x) x$P.stop.cum[c(1,length(vn.int)+1)])
  summ=sapply(out, function(x) x$P.effic.cum)
  rownames(summ)= paste0("Pat.",c(vn.int,n))
  if (printit==TRUE) print(summ)
  res=list(summary=summ, all=out)
  invisible(res)
}


