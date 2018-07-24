plotExpectedNumberVSptrue=function(n,vn.int,pF,cF,pE,cE,pvec,shape1F,shape2F, col=1,add=FALSE,show=TRUE,...){

  out=lapply(pvec, function(p) {expectedNumberFstop(n, p, pF, vn.int, cF, shape1=shape1F,shape2=shape2F)})

  if (show){
    if (!add)   plot(pvec,out,xlab=expression(p["true"]),ylab="Expected number of patients",type="n",
                       sub=paste("Interim analyses at",paste(vn.int,collapse=", "),", pF =",pF,", cF=",cF,", pE =",pE,", cE=",cE ),
                     las=1,...)
    lines(pvec,out,lwd=2,col=col)
  }
  class(out)="ptrue_vs_expectedNumber"
  invisible(out)

}
