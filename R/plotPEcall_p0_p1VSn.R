plotPEcall_p0_p1VSn=function(nvec,vn.int,pF,cF,pE,cE,p0,p1,shape1F,shape2F,shape1E=NULL,shape2E=NULL,col1="green",col2="red",cex.legend=1,show=TRUE,...){
  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F

  if (is.unsorted(vn.int) == TRUE) stop('vn.int must be sorted in ascending order')

  out=lapply(nvec, function(n) {
    if (n < min(vn.int)) stop('at least one interim necessary for calculation')
    vn.int1=vn.int[which(vn.int < n)]
    v.critE <- critEF(n=n,  vn.int=vn.int1, crit=cE, pE=pE,pF=pF,  EF="E",
                      shape1=shape1E,shape2=shape2E)
    v.critF <- critEF(n=n,  vn.int=vn.int1, crit=cF, pE=pE,pF=pF,  EF="F",
                      shape1=shape1F,shape2=shape2F)

#    res=c(pFstopEcall(p,c(vn.int1,n),v.critE,v.critF)$P.effic[length(vn.int1)+1],
#          pFstopEcall(p,c(vn.int1,n),v.critE,v.critF)$P.futil.cum[length(vn.int1)+1])

    res=c(pFstopEcall(p0,c(vn.int1,n),v.critE,v.critF)$P.effic[length(vn.int1)+1],
          pFstopEcall(p1,c(vn.int1,n),v.critE,v.critF)$P.effic[length(vn.int1)+1])


     res

  }
  )
  names(out)=nvec

  effp0=sapply(out, function(x) x[1])
  effp1=sapply(out, function(x) x[2])
  res=list(effp0=effp0,effp1=effp1)

  plot(as.numeric(names(effp0)),effp0,xlab="n",
       ylab="Probability",
       main=paste("Interim analyses at",paste(vn.int,collapse=", "),", final at n \npF =",pF,", cF=",cF,", pE =",pE,", cE=",cE ),
       ylim=c(0,1),pch=20,col=1,xlim=c(min(nvec),max(nvec)), type="n",las=1,...)
  lines(as.numeric(names(effp0)),effp0,lwd=2, col=col2)
  lines(as.numeric(names(effp1)),effp1,lwd=2, col=col1)
  legend("topleft",legend=c("Efficacy at final for ptrue=p0","Efficacy at final for ptrue=p1" ),lty=1,col = c(col2,col1),cex=par()$cex*cex.legend)

  class(res)="n_vs_pEcall_p0_p1"
  invisible(res)
}
