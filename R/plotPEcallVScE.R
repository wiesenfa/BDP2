plotPEcallVScE=function(n,vn.int,pF,cF,pE,cEvec,p0=NULL,p1=NULL,shape1F,shape2F,shape1E=NULL,shape2E=NULL,col0="red",col1="green",...){
  if (is.null(shape1E)) shape1E=shape1F
  if (is.null(shape2E)) shape2E=shape2F

  vn.int=vn.int[vn.int<n]

  pEcall_p0 <-pFstopEcall_vec(n=n, p=p0, pE=pE, pF=pF, vn.int=vn.int, cEvec, cF,
                              shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E)
  pEcall_p1 <-pFstopEcall_vec(n=n, p=p1, pE=pE, pF=pF, vn.int=vn.int, cEvec, cF,
                              shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E)

  plot(as.numeric(colnames(pEcall_p0$summary)),pEcall_p0$summary[length(vn.int)+1,],xlab=expression(c[E]),
       ylab="Probability of Efficacy at Final",
  #     sub=paste("Analyses at",paste(c(vn.int,n),collapse=", ") ),
        sub=paste0("Interim analyses at (",paste(vn.int,collapse=", "),"), final analysis at ",n,", pF=",pF,", cF=",cF,", pE =",pE),
       ylim=c(0,1),pch=20,col=1,xlim=c(min(cEvec),max(cEvec)), type="n",las=1,...) #final
  lines(as.numeric(colnames(pEcall_p0$summary)),pEcall_p0$summary[length(vn.int)+1,],
        col=col0,lwd=2)
  lines(as.numeric(colnames(pEcall_p1$summary)),pEcall_p1$summary[length(vn.int)+1,],
        col=col1,lwd=2)
  legend("bottomleft",legend=c(paste("ptrue=",p0),paste("ptrue=",p1)),lty=1,col = c(col0,col1))
  invisible(list(pEcall_p0=pEcall_p0,pEcall_p1=pEcall_p1,
                 x.p0=as.numeric(colnames(pEcall_p0$summary)),y.p0=pEcall_p0$summary[length(vn.int)+1,],
                 x.p1=as.numeric(colnames(pEcall_p1$summary)),y.p1=pEcall_p1$summary[length(vn.int)+1,]
                 ))
}



