plotExpectedNumberVSptrue_Estop=function(n,vn.int,pF,cF,pE,cE,pvec,shape1F,shape2F,shape1E=NULL,shape2E=NULL, col=1,add=FALSE,show=TRUE,...){

  
  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F
  
  v.critE <- critEF(n=n,  vn.int=vn.int, crit=cE, pE=pE,pF=pF,  EF="E", shape1=shape1E, shape2=shape2E)
  v.critF <- critEF(n=n,  vn.int=vn.int, crit=cF, pE=pE,pF=pF,  EF="F", shape1=shape1F, shape2=shape2F)
  
  out=lapply(pvec, function(p) {expectedNumberFstopEstop(n, p, pF, vn.int, cF, pE, cE,
                                                         shape1F, shape2F,shape1E, shape2E)})

  if (show){
    if (!add)   plot(pvec,out,xlab=expression(p["true"]),ylab="Expected number of patients",type="n",
                       sub=paste("Interim analyses at",paste(vn.int,collapse=", "),", pF =",pF,", cF=",cF,", pE =",pE,", cE=",cE ),
                     las=1,...)
    lines(pvec,out,lwd=2,col=col)
  }
  class(out)="ptrue_vs_expectedNumber"
  invisible(out)

}
