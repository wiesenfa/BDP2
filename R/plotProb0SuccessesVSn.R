plotProb0SuccessesVSn=function(nmin,nmax,p0,p1,cex.legend=1,...){
  plot(c(nmin,nmax), c(0,1), type="n",xlab="first interim at n=",ylab="Probability of 0 successes out of n",
       xaxt="n",main=paste("Probability of 0 successes out of n for p0=",p0," and p1=",p1),yaxs="i",las=1,...)
  axis(side=1,at=c(nmin:nmax))
  abline(h=c(0.01,0.05,0.1), col="grey")
  abline(v=c(nmin+1:nmax-1), col="grey",lty=2)
  lines(c(nmin:nmax),dbinom(0,c(nmin:nmax),p1),col="red",lwd=3,lty=1)
  lines(c(nmin:nmax),dbinom(0,c(nmin:nmax),p0),col="green",lwd=3,lty=1)

  legend("topright",legend=c(paste("ptrue=",p0),paste("ptrue=",p1)), lty=1, lwd=3,col=c("green","red"),cex=par()$cex*cex.legend)

}
