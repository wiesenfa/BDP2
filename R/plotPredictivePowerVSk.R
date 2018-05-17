plotPredictivePowerVSk=function(n,vn.int,pF,cF,pE,cE,shape1F,shape2F,shape1E=NULL,shape2E=NULL,show=TRUE,...){
  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F
  
  if (is.unsorted(vn.int) == TRUE) stop('vn.int must be sorted in ascending order')

  if (n < min(vn.int)) stop('at least one interim necessary for calculation')
  
  vn.int1=vn.int[which(vn.int < n)]
  v.critE <- critEF(n=n,  vn.int=vn.int1, crit=cE, pE=pE,pF=pF,  EF="E",
                    shape1=shape1E,shape2=shape2E)
  v.critF <- critEF(n=n,  vn.int=vn.int1, crit=cF, pE=pE,pF=pF,  EF="F",
                    shape1=shape1F,shape2=shape2F)
  
  n.evaluations=length(vn.int1)
  
  plot(c(1,vn.int1[n.evaluations]),c(0,1), type="n",xlab="k (successes at interim)", ylab="Predictive Power",
       sub=paste0("Interim analyses at (",paste(vn.int1,collapse=", "),"), final analysis at ",n,", pF=",pF,", cF=",cF,", pE =",pE,", cE=",cE))
  
  for (nint in c(1:n.evaluations)){
    for (k.int in c(0:vn.int1[nint])){
      n.int=vn.int1[nint]
      y=predictivePower(n.int,k.int,vn.int=c(vn.int1,n),v.critE,v.critF,shape1E,shape2E)
      points(k.int,y[1],col=nint,pch=nint)
    }
  }

  legend("bottomright",legend=paste(rep("interim",n.evaluations),c(1:n.evaluations)),col = c(1:n.evaluations),pch=c(1:n.evaluations))
  
#  class(res)="n_vs_pEcall_p0_p1"
#  invisible(res)
}
