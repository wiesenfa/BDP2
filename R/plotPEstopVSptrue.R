plotPEstopVSptrue=function(n,vn.int,pF,cF,pE,cE,pvec,shape1F,shape2F,shape1E=NULL,shape2E=NULL, col=1, add=FALSE,show=TRUE,...){
  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F

  v.critE <- critEF(n=n,  vn.int=vn.int, crit=cE, pE=pE,pF=pF,  EF="E",
                    shape1=shape1E,shape2=shape2E)
  v.critF <- critEF(n=n,  vn.int=vn.int, crit=cF, pE=pE,pF=pF,  EF="F",
                    shape1=shape1F,shape2=shape2F)

  out=lapply(pvec, function(p) {pFstopEstop(p,c(vn.int,n),v.critE,v.critF)})

  names(out)=pvec

  summ=sapply(out, function(x) x$P.effic.cum)
  rownames(summ)= paste0("Pat.",c(vn.int,n))
  res=list(summary=summ, all=out)

  if (!add) plot(as.numeric(colnames(res$summary)),res$summary[length(vn.int)+1,],xlab=expression(p["true"]),
                 ylab="Cumulative Probability of Efficacy at Final",
#                 main=paste("Analyses at",paste(c(vn.int,n),collapse=", "),"\npF =",pF,", cF=",cF,", pE =",pE,", cE=",cE ),
                     sub=paste("Interim analyses at",paste(vn.int,collapse=", "),", pF =",pF,", cF=",cF,", pE =",pE,", cE=",cE ),
               ylim=c(0,1),pch=20,xlim=c(min(pvec),max(pvec)), type="n",las=1,...) #final
  lines(as.numeric(colnames(res$summary)),res$summary[length(vn.int)+1,],lwd=2,col=col)

  class(res)="ptrue_vs_pEstop"
  invisible(res)


}
