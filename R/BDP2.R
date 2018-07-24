BDP2=function(n, interim.at, ptrue, eff.stop=FALSE, pF,cF,pE=NULL,  cE=NULL,
              type="PostProb", alpha=0.05,
                                       shape1F, shape2F,shape1E=NULL, shape2E=NULL,simulate=FALSE,nsim=10000){

#####################
  type
  alpha
  
############

  if (!simulate){
    if (type!="PostProb") stop("type must be 'PostProb', decisions based on predictive power only available for simulate==TRUE")
    if (is.na(eff.stop) | eff.stop==FALSE)  BDP2_Fstop(n=n, vn.int=interim.at, p=ptrue,
                                  pF=pF, cF=cF,
                                  shape1=shape1F,shape2=shape1F)
    else if (eff.stop=="call") BDP2_FstopEcall(n=n, vn.int=interim.at, p=ptrue,
                                               pF=pF,cF=cF,pE=pE,cE=cE,
                                               shape1F=shape1F, shape2F=shape2F,shape1E=shape1E, shape2E=shape2E)
    else if (eff.stop=="stop") BDP2_FstopEstop(n=n, vn.int=interim.at, p=ptrue,
                                               pF=pF,cF=cF,pE=pE,cE=cE,
                                               shape1F=shape1F, shape2F=shape2F,shape1E=shape1E, shape2E=shape2E)

  } else {
    if (type=="PostProb")  BDP2_simulateEF_PostProb(n=n, vn.int=interim.at, p=ptrue,
                                                     pF=pF,cF=cF,pE=pE, cE=cE,
                                                    nsim=nsim, eff.stop = (eff.stop=="stop"),
                    shape1F=shape1F, shape2F=shape2F,shape1E=shape1E, shape2E=shape2E)
    else if (type=="PredictivePower") BDP2_simulateEF_PredPow(n=n, vn.int=interim.at, p=ptrue, 
                            p0=pF,cF=cF,p1=pE, cE=cE, 
                            alpha=alpha, 
                            nsim=nsim, eff.stop = (eff.stop=="stop"),
              shape1F=shape1F, shape2F=shape2F,shape1E=shape1E, shape2E=shape2E) 
    else stop("type must be either 'PostProb' or 'PredictivePower'")

    

  }

}





