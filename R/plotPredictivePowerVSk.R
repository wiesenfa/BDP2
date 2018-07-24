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

  for (nint in c(1:n.evaluations)){
    n.int=vn.int1[nint]
    plot(c(1,n.int),c(0,1), type="n",xlab="", ylab="Predictive Power",xaxt="n")
    axis(1, c(0:n.int),  labels=c(0:n.int))
    axis(1, c(0:n.int),  labels=c(0:n.int)/n.int, line=1, lty=0, cex.axis=par()$cex.axis-.3)
    mtext('k (successes at interim)', 1, line=3)
    mtext('proportion of successes at interim', 1, line=4, cex=par()$cex.lab-.3)
    title(main=paste0("Interim analyses at (",paste(vn.int1,collapse=", "),"), final analysis at ",n,", pF=",pF,", cF=",cF,", pE =",pE,", cE=",cE),
          cex.main=par()$cex.main-.2,font.main=1)
        title(main=paste0("Interim at ",vn.int1[nint]),
          cex.main=par()$cex.main-.2,font.main=2,line=.7)
  
    y=sapply(0:n.int, function(k.int) predictivePower(n.int,k.int,vn.int=c(vn.int1,n),v.critE,v.critF,shape1E,shape2E)[,"PredPow_wFut"])
    points(0:n.int,y,col="black",pch=4)
  }

 
#  class(res)="n_vs_pEcall_p0_p1"
#  invisible(res)
}
