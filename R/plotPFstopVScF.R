plotPFstopVScF=function(n,vn.int,pF,cFvec,p0=NULL,p1=NULL,shape1,shape2,col1="green",col2="red",cex.legend=1,...){
  pFstop_p0 <-pFstop_vec(n=n, p=p0, pF=pF, vn.int=vn.int, cFvec=cFvec, shape1=shape1,shape2=shape2)
  pFstop_p1 <-pFstop_vec(n=n, p=p1, pF=pF, vn.int=vn.int, cFvec=cFvec, shape1=shape1,shape2=shape2)

  plot(as.numeric(colnames(pFstop_p0$summary)),pFstop_p0$summary[length(vn.int)+1,],xlab=expression(c[F]),
       ylab="Cumulative Stopping Probability at Final",main=paste("Analyses at",paste(c(vn.int,n),collapse=", ") ),
       ylim=c(0,1),xlim=c(min(cFvec),max(cFvec)), type="n",las=1,...) #final
  lines(as.numeric(colnames(pFstop_p0$summary)),pFstop_p0$summary[length(vn.int)+1,],
        col=col1,lwd=2)
  lines(as.numeric(colnames(pFstop_p1$summary)),pFstop_p1$summary[length(vn.int)+1,],
        col=col2,lwd=2)
  legend("bottomright",legend=c(paste("ptrue=",p0),paste("ptrue=",p1)),lty=1,col = c(col1,col2),cex=par()$cex*cex.legend)
}

