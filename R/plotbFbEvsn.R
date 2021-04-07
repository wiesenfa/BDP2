plotbFbEVSn=function(n,pF,cF,pE,cE, shape1F,shape2F,shape1E,shape2E, colF="red",colE="green",...){
  pF=pF
  pE=pE
  cF=cF
  cE=cE
  bF=CritBoundaryF5(n=n, pF=pF, critF=cF, shape1=shape1F,shape2=shape2F)
  bE=CritBoundaryE5(n=n, pE=pE, critE=cE, shape1=shape1E,shape2=shape2E)

  t_col <- function(color, percent = 50, name = NULL) {
    #	  color = color name
    #	percent = % transparency
    #	   name = an optional name for the color
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 maxColorValue =  255,
                 alpha = (100-percent)*255/100,
                 names = name)
    ## Save the color
    invisible(t.col)
     }

  mycolF <- t_col(colF, perc = 50, name = "lt.colF")
  mycolE <- t_col(colE, perc = 50, name = "lt.colE")

  plot(c(0:n),c(0:n), ylim=c(-0.5,max(c(bF,bE))),yaxt="n",xaxs="i",yaxs="i",type="n",
       xlab="Number of enrolled patients",ylab="Decision boundaries",
       main=paste0("Boundaries for futility (",colF,") and efficacy (",colE,")"),
       sub=paste0("pF =",pF,", cF=",cF,", pE=",pE,", cE=",cE ),... )
  axis(side=2,at=c(0:max(c(bF,bE))))
  lines(c(1:n),bF, col=colF,lwd=3)
  lines(c(1:n),bE, col=colE,lwd=3)
  abline(v=0:n,lty=3,lwd=1,col="grey")
  abline(h=c(0:max(c(bF,bE))), lty=3, lwd=1,col="grey")
  polygon(c(1,c(1:n),n),c(-1,bF,-1), border = NA,col=mycolF)
  polygon(c(1,c(1:n),n),c(max(c(bF,bE)),bE,max(c(bF,bE))), border = NA,col=mycolE)
}



