plotPFstopEcallVSn=function(nvec,vn.int,pF,cF,pE,cE,p,shape1F,shape2F,shape1E=NULL,shape2E=NULL,col1="green",col2="red",progress=FALSE,cex.legend=1,show=TRUE,...){
  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F

  if (is.unsorted(vn.int) == TRUE) stop('vn.int must be sorted in ascending order')

  if (progress)  withProgress(message = 'Compute Operating Characteristics', value = 0,{
      out=lapply(nvec, function(n) {
      if (n < min(vn.int)) stop('at least one interim necessary for calculation')
      incProgress(1/nvec, detail = paste("for n=", n))
      vn.int1=vn.int[which(vn.int < n)]
      v.critE <- critEF(n=n,  vn.int=vn.int1, crit=cE, pE=pE,pF=pF,  EF="E",
                        shape1=shape1E,shape2=shape2E)
      v.critF <- critEF(n=n,  vn.int=vn.int1, crit=cF, pE=pE,pF=pF,  EF="F",
                        shape1=shape1F,shape2=shape2F)
      
      temp=pFstopEcall(p,c(vn.int1,n),v.critE,v.critF)
      res=c(temp$P.effic[length(vn.int1)+1],
            temp$P.futil.cum[length(vn.int1)+1])
      res

    })

  }) else   out=lapply(nvec, function(n) {
      if (n < min(vn.int)) stop('at least one interim necessary for calculation')
       vn.int1=vn.int[which(vn.int < n)]
      v.critE <- critEF(n=n,  vn.int=vn.int1, crit=cE, pE=pE,pF=pF,  EF="E",
                        shape1=shape1E,shape2=shape2E)
      v.critF <- critEF(n=n,  vn.int=vn.int1, crit=cF, pE=pE,pF=pF,  EF="F",
                        shape1=shape1F,shape2=shape2F)
      
      temp=pFstopEcall(p,c(vn.int1,n),v.critE,v.critF)
      res=c(temp$P.effic[length(vn.int1)+1],
            temp$P.futil.cum[length(vn.int1)+1])
      res

    })
  names(out)=nvec



  eff=sapply(out, function(x) x[1])
  fut=sapply(out, function(x) x[2])
  res=list(eff=eff,fut=fut)

  if (show){
    plot(as.numeric(names(eff)),eff,xlab="final sample size (n)",
         ylab="Probability",
         sub=paste0("Interim analyses at (",paste(vn.int,collapse=", "),"),  ptrue=",p,", pF=",pF,", cF=",cF,", pE =",pE,", cE=",cE ),
         ylim=c(0,1),xlim=c(min(nvec),max(nvec)), type="n",las=1,...)
#    abline(v=min(nvec):max(nvec),lty=3,lwd=1,col="grey")
    points(as.numeric(names(eff)),eff,pch=4, col=col1)
    points(as.numeric(names(fut)),fut,pch=3, col=col2)
    legend("topleft",legend=c("Efficacy at Final","Cumulative Futility at Final" ),pch=c(4,3),col = c(col1,col2),cex=par()$cex*cex.legend)
  }

  class(res)="n_vs_pFstopEcall"
  invisible(res)
}
