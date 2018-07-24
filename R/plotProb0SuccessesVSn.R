plotProb0SuccessesVSn=function(nmin,nmax,p0,p1,cex.legend=1,...){
  plot(c(nmin,nmax), c(0,1), type="n",xlab="first interim at n=",ylab="Probability of 0 successes out of n",
       xaxt="n",main=paste("Probability of 0 successes out of n for p0=",p0," and p1=",p1),yaxs="i",las=1,...)
  axis(side=1,at=c(nmin:nmax))
  abline(h=c(0.01,0.05,0.1), col="grey")
#  abline(v=c(nmin+1:nmax-1), col="grey",lty=2)
  points(c(nmin:nmax),dbinom(0,c(nmin:nmax),p1),col="red",pch=4)
  points(c(nmin:nmax),dbinom(0,c(nmin:nmax),p0),col="green",pch=3)

  legend("topright",legend=c(paste("ptrue=",p0),paste("ptrue=",p1)), col=c("green","red"),pch=c(3:4),cex=par()$cex*cex.legend)

}
