pinconsistent=function(p,vn.int,v.critE,v.critF){
  
  # Probability of inconsistent calls in the case of stopping for futility and calling efficacy at interims
  # a call is ionconsistent if a trial is called efficacious at an interim 
  #     and stopped for futility at a later interim
  
  # Args:
  # p:           true response rate
  # vn.int:      vector of sample sizes at the interim analyses 
  # v.critE:     vector of critical boundaries for calling efficacy at the interim analyses 
  # v.critF:     vector of critical boundaries for futility stopping at the interim analyses 
  
  
  n.interim=length(vn.int)
  if (n.interim>5) stop("max 5 interims including final")
  pincon <- rep(0,n.interim)
  pincon_3_1 <- 0
  pincon_3_2 <- 0
  pincon_4_1 <- 0
  pincon_4_2 <- 0
  pincon_4_3 <- 0
  pincon_5_1 <- 0
  pincon_5_2 <- 0
  pincon_5_3 <- 0
  pincon_5_4 <- 0
 
  pincon[1] <- 0
  
  #Eff in [1], Fut in [2]   
  if (n.interim>1)  
    for (my.i1 in v.critE[1]:vn.int[1]) {
      for (my.i2 in 0:(v.critF[2] - my.i1)) { 
        pincon[2] <- pincon[2] + 
          dbinom(size=vn.int[1],x=my.i1,p)*
          (v.critF[2] - my.i1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)
      }
    }

##3 interims
 #Eff in [1], keine Fut in [2], Fut in [3] 
  if (n.interim>2)  
    for (my.i1 in v.critE[1] :vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2] - vn.int[1])) {
        for (my.i3 in 0:(v.critF[3] - my.i1 - my.i2 )) {   
          pincon_3_1 <- pincon_3_1 + 
            
            dbinom(size=vn.int[1],x=my.i1,p)*
            (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (v.critF[3] - my.i1 - my.i2  >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)
          
        }
      }
    }
 
 #keine Fut keine Eff in [1], Eff in [2], Fut in [3] 
  if (n.interim>2)  
    for (my.i1 in (v.critF[1] + 1) :(v.critE[1] - 1)) {
      for (my.i2 in (v.critE[2] - my.i1):(vn.int[2] - vn.int[1]  )) {
        for (my.i3 in 0:(v.critF[3] - my.i1 - my.i2 )) {   
          pincon_3_2 <- pincon_3_2 + 
            
            dbinom(size=vn.int[1],x=my.i1,p)*
            (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (v.critF[3] - my.i1 - my.i2  >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)
          
        }
      }
    }
if (n.interim>2)  pincon[3] = pincon_3_1 + pincon_3_2
  
##4 interims
#Eff in [1], keine Fut in [2], keine Fut in [3], Fut in [4]
if (n.interim>3)  
  for (my.i1 in v.critE[1] :vn.int[1]) {
    for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2] - vn.int[1])) {
      for (my.i3 in (v.critF[3]+1 - my.i1):(vn.int[3] - vn.int[2])) {
        for (my.i4 in 0:(v.critF[4] - my.i1 - my.i2 - my.i3)) {   
        pincon_4_1 <- pincon_4_1 + 
          
          dbinom(size=vn.int[1],x=my.i1,p)*
          (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
          (my.i3 >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
          (v.critF[4] - my.i1 - my.i2 - my.i3 >= 0) * dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)
        }       
      }
    }
  }

#keine Fut keine Eff in [1], Eff in [2], keine Fut in [3],  Fut in [4] 
if (n.interim>3)  
  for (my.i1 in (v.critF[1] + 1) :(v.critE[1] - 1)) {
    for (my.i2 in (v.critE[2] - my.i1):(vn.int[2] - vn.int[1]  )) {
      for (my.i3 in (v.critF[3] + 1 - my.i1 - my.i2):(vn.int[3] - vn.int[2]  )) {
        for (my.i4 in 0:(v.critF[4] - my.i1 - my.i2 - my.i3)) {   
          pincon_4_2 <- pincon_4_2 + 
            
            (v.critE[1] - 1 - v.critF[1] - 1 >= 0) * dbinom(size=vn.int[1],x=my.i1,p)*
            (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (my.i3 >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
            (v.critF[4] - my.i1 - my.i2 - my.i3 >= 0) * dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)
        }        
      }
    }
  }

#keine Fut keine Eff in [1], keine Fut keine Eff in [2], Eff in [3],  Fut in [4] 
if (n.interim>3)  
  for (my.i1 in (v.critF[1] + 1) :(v.critE[1] - 1)) {
    for (my.i2 in (v.critF[2] + 1 - my.i1) :(v.critE[2] - 1 - my.i1)) {
      for (my.i3 in (v.critE[3]  - my.i1 - my.i2):(vn.int[3] - vn.int[2]  )) {
        for (my.i4 in 0:(v.critF[4] - my.i1 - my.i2 - my.i3)) {   
          pincon_4_2 <- pincon_4_2 + 
            
            (v.critE[1] - 1 - v.critF[1] - 1 >= 0) * dbinom(size=vn.int[1],x=my.i1,p)*
            (v.critE[2] - 1 - v.critF[2] - 1 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (my.i3 >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
            (v.critF[4] - my.i1 - my.i2 - my.i3 >= 0) * dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)
        }        
      }
    }
  }

if (n.interim>3)  pincon[4] = pincon_4_1 + pincon_4_2 + pincon_4_3

##5 interims
#Eff in [1], keine Fut in [2], keine Fut in [3], keine Fut in [4], Fut in [5]
if (n.interim>4)  
  for (my.i1 in v.critE[1] :vn.int[1]) {
    for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2] - vn.int[1])) {
      for (my.i3 in (v.critF[3]+1 - my.i1):(vn.int[3] - vn.int[2])) {
        for (my.i4 in (v.critF[4]+1 - my.i1):(vn.int[4] - vn.int[3])) {
          for (my.i5 in 0:(v.critF[5] - my.i1 - my.i2 - my.i3 -my.i4)) {   
          pincon_5_1 <- pincon_5_1 + 
            
            dbinom(size=vn.int[1],x=my.i1,p)*
            (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (my.i3 >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
            (my.i4 >= 0) * dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
            (v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4>= 0) * dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)
          }
        }       
      }
    }
  }
#keine Fut keine Eff in [1], Eff in [2], keine Fut in [3], keine Fut in [4],  Fut in [5] 
  if (n.interim>4)  
    for (my.i1 in (v.critF[1] + 1) :(v.critE[1] - 1)) {
      for (my.i2 in (v.critE[2] - my.i1):(vn.int[2] - vn.int[1]  )) {
        for (my.i3 in (v.critF[3] + 1 - my.i1 - my.i2):(vn.int[3] - vn.int[2]  )) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(vn.int[4] - vn.int[3]  )) {
            for (my.i5 in 0:(v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4)) {   
            pincon_5_2 <- pincon_5_2 + 
              
              (v.critE[1] - 1 - v.critF[1] - 1 >= 0) * dbinom(size=vn.int[1],x=my.i1,p)*
              (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
              (my.i3 >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
              (my.i4 >= 0) * dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
              (v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4>= 0) * dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)
            }
          }        
        }
      }
    }
#keine Fut keine Eff in [1], keine Fut keine Eff in [2], Eff in [3],  keine Fut in [4],  Fut in [5]  
if (n.interim>4)  
  for (my.i1 in (v.critF[1] + 1) :(v.critE[1] - 1)) {
    for (my.i2 in (v.critF[2] + 1 - my.i1) :(v.critE[2] - 1 - my.i1)) {
      for (my.i3 in (v.critE[3]  - my.i1 - my.i2):(vn.int[3] - vn.int[2]  )) {
        for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(vn.int[4] - vn.int[3]  )) {
          for (my.i5 in 0:(v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4)) {   
          pincon_5_3 <- pincon_5_3 + 
            
            (v.critE[1] - 1 - v.critF[1] - 1 >= 0) * dbinom(size=vn.int[1],x=my.i1,p)*
            (v.critE[2] - 1 - v.critF[2] - 1 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (my.i3 >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
            (my.i4 >= 0) * dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
            (v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4 >= 0) * dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)
          }
        }        
      }
    }
  }
  #keine Fut keine Eff in [1], keine Fut keine Eff in [2],keine Fut keine Eff in [3], Eff in [4],  Fut in [5]  
  if (n.interim>4)  
    for (my.i1 in (v.critF[1] + 1) :(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1) :(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1 - my.i1 - my.i2) :(v.critE[3] - 1 - my.i1 - my.i3)) {
            for (my.i4 in (v.critE[4] - my.i1 - my.i2 - my.i3):(vn.int[4] - vn.int[3]  )) {
              for (my.i5 in 0:(v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4)) {   
              pincon_5_4 <- pincon_5_4 + 
                
                (v.critE[1] - 1 - v.critF[1] - 1 >= 0) * dbinom(size=vn.int[1],x=my.i1,p)*
                (v.critE[2] - 1 - v.critF[2] - 1 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                (v.critE[3] - 1 - v.critF[3] - 1 >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                (my.i4 >= 0) * dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                (v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4 >= 0) * dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)
            }        
          }
        }
      }
    }   

if (n.interim>4)  pincon[5] = pincon_5_1 + pincon_5_2 + pincon_5_3 + pincon_5_4




 
probcall <- rep(0,n.interim)
probstop <- rep(0,n.interim)

probstop[1] <- pbinom(size=vn.int[1],q=v.critF[1] ,p, lower.tail=TRUE)

if (n.interim>1)  
  for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    for (my.i2 in 0:(v.critF[2] - my.i1)) { 
      probstop[2] <- probstop[2] + 
        dbinom(size=vn.int[1],x=my.i1,p)*
        (v.critF[2] -my.i1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)
    }
  }

if (n.interim>2)  
  for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
      for (my.i3 in 0:(v.critF[3] - my.i1 - my.i2 )) {   
        probstop[3] <- probstop[3] + 
          
          dbinom(size=vn.int[1],x=my.i1,p)*
          (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
          (v.critF[3] - my.i1 - my.i2  >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)
        
      }
    }
  }

if (n.interim>3)  
  for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
      for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
        for (my.i4 in 0:(v.critF[4]- my.i1 - my.i2 - my.i3 )) {
          probstop[4] <- probstop[4] + 
            
            dbinom(size=vn.int[1],x=my.i1,p)*
            (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
            (v.critF[4] - my.i3 - my.i2 - my.i1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) 
          
        }
      }
    }
  }

if (n.interim>4)  
  for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
      for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
        for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
          for (my.i5 in 0:(v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4)) {
            
            probstop[5] <- probstop[5] + 
              
              dbinom(size=vn.int[1],x=my.i1,p)*
              (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
              (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
              (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
              (v.critF[5]- my.i4 -my.i3 - my.i2 - my.i1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) 
          }
        }
      }
    }
  }



probcall[1] <- pbinom(size=vn.int[1],q=v.critE[1]-1 ,p, lower.tail=FALSE)

if (n.interim>1)  
  for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    for (my.i2 in (v.critE[2] - my.i1):(vn.int[2]-vn.int[1])) { 
      probcall[2] <- probcall[2] + 
        dbinom(size=vn.int[1],x=my.i1,p)*
        (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)
    }
  }

if (n.interim>2)  
  for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
      for (my.i3 in (v.critE[3] - my.i1 - my.i2):(vn.int[3]-vn.int[2])) {   
        probcall[3] <- probcall[3] + 
          
          dbinom(size=vn.int[1],x=my.i1,p)*
          (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
          (my.i3  >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)
        
      }
    }
  }

if (n.interim>3)  
  for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
      for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
        for (my.i4 in (v.critE[4] - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
          probcall[4] <- probcall[4] + 
            
            dbinom(size=vn.int[1],x=my.i1,p)*
            (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
            (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) 
          
        }
      }
    }
  }

if (n.interim>4)  
  for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
      for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
        for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
          for (my.i5 in (v.critE[5] - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
            
            probcall[5] <- probcall[5] + 
              
              dbinom(size=vn.int[1],x=my.i1,p)*
              (my.i2 >= 0) * dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
              (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
              (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
              (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) 
          }
        }
      }
    }
  }


  
  pinconsistent=data.frame(Int.Analysis=1:length(vn.int),Pat.number=vn.int,Crit.boundaryE=v.critE,
                            P.effic=round(probcall, digits = 4),
                            Crit.boundaryF=v.critF,P.futil=round(probstop, digits = 4),
                            P.futil.cum=round(cumsum(probstop), digits = 4),
                            P.inconsistent=round(pincon, digits = 4),
                            P.inconsistent.cum=round(cumsum(pincon), digits = 4))

  pinconsistent
  
}




