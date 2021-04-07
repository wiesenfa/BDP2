plotBDP2=function(
                    x=c("n","k", "ptrue","cE","cF"),
                    y=c("Prob0Successes","PostProb0or1Successes","bFbE",
                        "PEcall_p0_p1","PEstop_p0_p1",
                        "PFstopEcall","PFstopEstop",  #if x=="n"
                        "PEcall","PEstop","PFstop","PFstopEstop",
                        "ExpectedNumber",
                        "PredictivePower"),   #if x in ptrue, cE, cF
                    n,   #could be vector or scalar
                    interim.at,
                     ptrue, #could be vector or scalar
                   pF, cF,
                    pE, cE,  #cF, cE could be vector or scalar
                    p0,p1,
                   Estop=FALSE,
                    shape1F,shape2F,shape1E=NULL,shape2E=NULL,
                    col=c("green","red"), #could be vector or scalar
                    cex.legend=1,add=FALSE,show=TRUE,progress=FALSE,...){

  args=as.list(match.call(expand.dots = FALSE))
  x=match.arg(x)
  y=match.arg(y)
  
  switch(x,
         n= switch(y,
                  Prob0Successes=  plotProb0SuccessesVSn(nmin=n[1],nmax=n[length(n)],
                                                          p0=p0,p1=p1,cex.legend=cex.legend,...),

                  PostProb0or1Successes=  plotPostProb0or1SuccessesVSn(nmin=n[1],nmax=n[length(n)],
                                                                      pF=pF,
                                                                      shape1F=shape1F,shape2F=shape2F,cex.legend=cex.legend,...),
                  bFbE=   plotbFbEVSn(n=n,
                                       pF=pF,cF=cF,pE=pE,cE=cE,
                                       shape1F,shape2F,shape1E,shape2E,
                                       colF=col[1],colE=col[2],...),
                  PEcall_p0_p1=  plotPEcall_p0_p1VSn(nvec=n,vn.int=interim.at,
                                                     pF=pF,cF=cF,pE=pE,cE=cE,
                                                     p0=p0,p1=p1,
                                                     shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                                     col1=col[1],col2=col[2],cex.legend=cex.legend,show=show,...),
                 PEstop_p0_p1=  plotPEstop_p0_p1VSn(nvec=n,vn.int=interim.at,
                                                     pF=pF,cF=cF,pE=pE,cE=cE,
                                                     p0=p0,p1=p1,
                                                     shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                                     col1=col[1],col2=col[2],cex.legend=cex.legend,show=show,...),
                  
                  
                  PFstopEcall=  plotPFstopEcallVSn(nvec=n,vn.int=interim.at,
                                                   pF=pF,cF=cF,pE=pE,cE=cE,
                                                   p=ptrue,
                                                   shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                                   col1=col[1],col2=col[2],progress=progress,cex.legend=cex.legend,show=show,...),
                 PFstopEstop=  plotPFstopEstopVSn(nvec=n,vn.int=interim.at,
                                                   pF=pF,cF=cF,pE=pE,cE=cE,
                                                   p=ptrue,
                                                   shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                                   col1=col[1],col2=col[2],progress=progress,cex.legend=cex.legend,show=show,...)
               ),
         k=switch(y,
                    PredictivePower=plotPredictivePowerVSk(n=n,vn.int=interim.at,
                                                           pF=pF,cF=cF,pE=pE,cE=cE,
                                                           shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                                           show=show,...)
                  ),
         ptrue= switch(y,
                         PEcall=  plotPEcallVSptrue(n=n,vn.int=interim.at,    #vielleicht ProbCallEfficacy statt PEcall
                                               pF=pF,cF=cF,pE=pE,cE=cE,
                                               pvec=ptrue,
                                               shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                               col=col,add=add,show=show,...),
                         PEstop=  plotPEstopVSptrue(n=n,vn.int=interim.at,    #vielleicht ProbCallEfficacy statt PEcall
                                               pF=pF,cF=cF,pE=pE,cE=cE,
                                               pvec=ptrue,
                                               shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                               col=col,add=add,show=show,...),
                        ExpectedNumber= ifelse(Estop==TRUE,
                                               no = plotExpectedNumberVSptrue(n=n,vn.int=interim.at,
                                                   pF=pF,cF=cF,pE=pE,cE=cE,
                                                   pvec=ptrue,
                                                   shape1F=shape1F,shape2F=shape2F,
                                                   col=col,add=add,show=show,...),
                                               yes=plotExpectedNumberVSptrue_Estop(n=n,vn.int=interim.at,
                                                   pF=pF,cF=cF,pE=pE,cE=cE,
                                                   pvec=ptrue,
                                                   shape1F=shape1F,shape2F=shape2F,
                                                   col=col,add=add,show=show,...)),
                        PFstop= plotPFstopVSptrue(n=n,vn.int=interim.at,
                                             pF=pF,cF=cF,
                                             pvec=ptrue,
                                             shape1=shape1F,shape2=shape2F,
                                             col=col,add=add,show=show,...)
                      ),
          cE= switch(y,
                    PEcall= do.call("plotPEcallVScE",c(list(n=n,vn.int=interim.at,   #vielleicht ProbCallEfficacy statt PEcall
                                                       pF=pF,cF=cF,pE=pE,cEvec=cE,
                                                       p0=p0,p1=p1,
                                                       shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                                       col0=col[1],col1=col[2],cex.legend=cex.legend,show=show),args$...)),
                   PEstop= do.call("plotPEstopVScE",c(list(n=n,vn.int=interim.at,   #vielleicht ProbCallEfficacy statt PEcall
                                   pF=pF,cF=cF,pE=pE,cEvec=cE,
                                   p0=p0,p1=p1,
                                   shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                   col0=col[1],col1=col[2],cex.legend=cex.legend,show=show),args$...))

                    ),
          cF= switch(y,
                    PFstop= plotPFstopVScF(n=n,vn.int=interim.at,
                                           pF=pF,cFvec=cF,
                                           p0=p0,p1=p1,
                                           shape1=shape1F,shape2=shape2F,
                                           col1=col[1],col2=col[2]),cex.legend=cex.legend,...)
  )



}




