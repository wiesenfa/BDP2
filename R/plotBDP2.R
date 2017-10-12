plotBDP2=function(
                    x=c("n","ptrue","cE","cF"),
                    y=c("Prob0Successes","PostProb0or1Successes","bFbE","PEcall_p0_p1","PFstopEcall",  #if x=="n"
                        "PEcall","PFstop","ExpectedNumber"),   #if x in ptrue, cE, cF
                    n,   #could be vector or scalar
                    interim.at,
                     ptrue, #could be vector or scalar
                   pF, cF,
                    pE, cE,  #cF, cE could be vector or scalar
                    p0,p1,
                    shape1F,shape2F,shape1E=NULL,shape2E=NULL,
                    col=c("green","red"), #could be vector or scalar
                    add=FALSE,progress=FALSE,...){


  x=match.arg(x)
  y=match.arg(y)
  switch(x,
         n= switch(y,
                  Prob0Successes=  plotProb0SuccessesVSn(nmin=n[1],nmax=n[length(n)],
                                                          p0=p0,p1=p1,...),

                  PostProb0or1Successes=  plotPostProb0or1SuccessesVSn(nmin=n[1],nmax=n[length(n)],
                                                                      pF=pF,
                                                                      shape1F=shape1F,shape2F=shape2F,...),
                  bFbE=   plotbFbEVSn(n=n,
                                       pF=pF,cF=cF,pE=pE,cE=cE,
                                       shape1F,shape2F,shape1E,shape2E,
                                       colF=col[1],colE=col[2],...),
                  PEcall_p0_p1=  plotPEcall_p0_p1VSn(nvec=n,vn.int=interim.at,
                                                     pF=pF,cF=cF,pE=pE,cE=cE,
                                                     p0=p0,p1=p1,
                                                     shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                                     col1=col[1],col2=col[2],...),
                  PFstopEcall=  plotPFstopEcallVSn(nvec=n,vn.int=interim.at,
                                                   pF=pF,cF=cF,pE=pE,cE=cE,
                                                   p=ptrue,
                                                   shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                                   col1=col[1],col2=col[2],progress=progress,...)
               ),
         ptrue= switch(y,
                         PEcall=  plotPEcallVSptrue(n=n,vn.int=interim.at,    #vielleicht ProbCallEfficacy statt PEcall
                                               pF=pF,cF=cF,pE=pE,cE=cE,
                                               pvec=ptrue,
                                               shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                               col=col,add=add,...),
                        ExpectedNumber= plotExpectedNumberVSptrue(n=n,vn.int=interim.at,
                                                   pF=pF,cF=cF,pE=pE,cE=cE,
                                                   pvec=ptrue,
                                                   shape1F=shape1F,shape2F=shape2F,
                                                   col=col,add=add,...),
                        PFstop= plotPFstopVSptrue(n=n,vn.int=interim.at,
                                             pF=pF,cF=cF,
                                             pvec=ptrue,
                                             shape1=shape1F,shape2=shape2F,
                                             col=col,add=add,...)
                      ),
          cE= switch(y,
                    PEcall= plotPEcallVScE(n=n,vn.int=interim.at,   #vielleicht ProbCallEfficacy statt PEcall
                                                       pF=pF,cF=cF,pE=pE,cEvec=cE,
                                                       p0=p0,p1=p1,
                                                       shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,
                                                       col0=col[1],col1=col[2],...)
                    ),
          cF= switch(y,
                    PFstop= plotPFstopVScF(n=n,vn.int=interim.at,
                                           pF=pF,cFvec=cF,
                                           p0=p0,p1=p1,
                                           shape1=shape1F,shape2=shape2F,
                                           col1=col[1],col2=col[2]),...)
  )



}



