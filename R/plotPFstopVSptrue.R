plotPFstopVSptrue=function(n,vn.int,pF,cF,pvec,shape1,shape2, col=1, add=FALSE,show=TRUE,...){
  v.crit <-crit(n=n, pThreshold=pF, vn.int=vn.int, crit=cF, shape1 = shape1, shape2 = shape2)

  out=lapply(pvec, function(p) {pFstop(p,c(vn.int,n),v.crit)})

  names(out)=pvec

  summ=sapply(out, function(x) x$P.stop.cum)

  rownames(summ)= paste0("Pat.",c(vn.int,n))
  res=list(summary=summ, all=out)

  if (show){
     if (!add) plot(as.numeric(colnames(res$summary)),res$summary[length(vn.int)+1,],xlab=expression(p["true"]),
                   ylab="Cumulative Stopping Probability at Final",
                   main=paste("Analyses at",paste(c(vn.int,n),collapse=", "),"\npF =",pF,", cF=",cF ),
                   ylim=c(0,1),xlim=c(min(pvec),max(pvec)), type="n",las=1,...) #final
    lines(as.numeric(colnames(res$summary)),res$summary[length(vn.int)+1,],lwd=2,col=col)
   
  }
  class(res)="ptrue_vs_pFstop"
  invisible(res)
}
