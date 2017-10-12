BDP2=function(n, interim.at, ptrue, eff.stop=FALSE, pF,cF,pE=NULL,  cE=NULL,
                                       shape1F, shape2F,shape1E=NULL, shape2E=NULL,simulate=FALSE,nsim=10000){



  if (!simulate){
    if (is.na(eff.stop) | eff.stop==FALSE)  BDP2_Fstop(n=n, vn.int=interim.at, p=ptrue,
                                  pF=pF, cF=cF,
                                  shape1=shape1F,shape2=shape1F)
    else if (eff.stop=="call") BDP2_FstopEcall(n=n, vn.int=interim.at, p=ptrue,
                                               pF=pF,cF=cF,pE=pE,cE=cE,
                                               shape1F=shape1F, shape2F=shape2F,shape1E=NULL, shape2E=NULL)
    else if (eff.stop=="stop") BDP2_FstopEstop(n=n, vn.int=interim.at, p=ptrue,
                                               pF=pF,cF=cF,pE=pE,cE=cE,
                                               shape1F=shape1F, shape2F=shape2F,shape1E=NULL, shape2E=NULL)

  } else {
      BDP2_simulateEF_PostProb(n=n, vn.int=interim.at, p=ptrue,
                                                     pF=pF,cF=cF,pE=pE, cE=cE,
                                                    nsim=nsim, eff.stop = (eff.stop=="stop"),
                    shape1F=shape1F, shape2F=shape2F,shape1E=NULL, shape2E=NULL)

  }

}





