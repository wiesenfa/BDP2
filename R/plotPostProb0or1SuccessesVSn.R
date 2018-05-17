plotPostProb0or1SuccessesVSn=function(nmin,nmax,pF,shape1F,shape2F,cex.legend=1,...){
  ymax=max(1-pbeta(pF,shape1F,shape2F+nmin),1-pbeta(pF,shape1F+1,shape2F+nmin-1))
  plot(c(nmin,nmax),c(0,ymax), type="n",xlab="first interim at n=",ylab="P(p>pF|x successes out of n) (=cF)",
       xaxt="n",main=paste("Posterior probability of 0 or 1 successes out of n for pF=",pF),yaxs="i",las=1,...)
  axis(side=1,at=c(nmin:nmax))
  abline(h=c(0.01,0.05,0.1), col="grey")
#  abline(v=c(nmin+1:nmax-1), col="grey",lty=2)
  function0 = function(x){1-pbeta(pF,shape1F,shape2F+x)}
  function1 = function(x){1-pbeta(pF,shape1F+1,shape2F+x-1)}
  points(c(nmin:nmax),function0(nmin:nmax),col="red",pch=4)
  points(c(nmin:nmax),function1(nmin:nmax),col="green",pch=3)
#  plot(function(x) 1-pbeta(pF,shape1F,shape2F+x),nmin,nmax,add=T,col="red",lwd=3,lty=1)
#  plot(function(x) 1-pbeta(pF,shape1F+1,shape2F+x-1),nmin,nmax,add=T,col="green",lwd=3,lty=1)
  
  legend("topright",legend=c("1 success out of n","no successes out of n"), pch=c(3,4), col=c("green","red"),cex=par()$cex*cex.legend)

}
