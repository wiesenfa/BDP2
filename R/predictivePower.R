predictivePower=function(n.int,k.int,vn.int,v.critE,v.critF,shape1E,shape2E){
  # Futility stopping probability at every interim: Evaluate at vn.int[i] and stop if #successes <= v.critF[i]
  # Probability of Calling efficacy at final: Evaluate at vn.int[i] and call if #successes >= v.critE[length(vn.int)]

  # Args:
  # n.int:       sample size at the interim analsis
  # k.int:       number of successes observed up until the interim analysis
  # vn.int:      vector of sample sizes at the interim analyses
  # v.critE:     vector of critical boundaries for calling efficacy at the interim analyses
  # v.critF:     vector of critical boundaries for futility stopping at the interim analyses

  #use nomenclature from paper: final is the number of the final analysis in terms of number of interim analyses
  final=length(vn.int)
  if (final>7) stop("max 7 interims including final")
  
  #identify which interim corresponds to n.int
  if (identical(which(vn.int==n.int),integer(0))) stop("interim analysis not planned")
  current.interim=which(vn.int==n.int)
  
  my.dbetabinom=function(x,size,a,b){
    if (x < 0 | x > size | a < 0 | b < 0){
      xxx=0
     } else{
        xxx=choose(size, x)*beta(a + x, b + size - x)/beta(a, b)
     }
    return(xxx)
  }
  
  
  # if k.int was observed in final analysis
  if (current.interim == final){
    pp=ifelse(k.int >= v.critE(final),1,0)
  }
  
  if (current.interim == final-1){
    pp=0
    if (k.int > v.critF[final-1]){
    a1=shape1E+k.int
    b1=shape2E+n.int-k.int
    for (my.i in (v.critE[final] - k.int):(vn.int[final]-vn.int[final-1])) {
        pp <- pp +
          (my.i >= 0)*my.dbetabinom(x=my.i,size=vn.int[final]-vn.int[final-1],a=a1,b=b1)
     
      }
    }
  }
  
  if (current.interim == final-2){
    pp=0
    #my.i1 for (final-2 to) final-1
    #my.i for (final-1 to) final
    if (k.int > v.critF[final-2]){
    a1=shape1E+k.int
    b1=shape2E+n.int-k.int
   for (my.i1 in (v.critF[final-1] + 1 - k.int):(vn.int[final-1]-vn.int[final-2])) {
     a=shape1E+k.int + my.i1
     b=shape2E+n.int-k.int + vn.int[final-1]-vn.int[final-2] - my.i1
       for (my.i in (v.critE[final] - my.i1 - k.int):(vn.int[final]-vn.int[final-1])) {
         pp <- pp +
           (my.i1 >= 0)*
           my.dbetabinom(x=my.i1,size=vn.int[final-1]-vn.int[final-2],a=a1,b=b1)*
          (my.i >= 0)*
           my.dbetabinom(x=my.i,size=vn.int[final]-vn.int[final-1],a=a,b=b)
         }
      }
    }
  }
  
  if (current.interim == final-3){
    pp=0
    #my.i2 for (final-3 to) final-2
    #my.i1 for (final-2 to) final-1
    #my.i for (final-1 to) final
    if (k.int > v.critF[final-3]){
    a2=shape1E+k.int
    b2=shape2E+n.int-k.int
    
    for (my.i2 in (v.critF[final-2] + 1 - k.int):(vn.int[final-2]-vn.int[final-3])) {
      a1=shape1E+k.int + my.i2
      b1=shape2E+n.int-k.int + vn.int[final-2]-vn.int[final-3] - my.i2
      
      for (my.i1 in (v.critF[final-1] + 1 - my.i2 - k.int):(vn.int[final-1]-vn.int[final-2])) {
        a=shape1E+k.int + my.i2 + my.i1 
        b=shape2E+n.int-k.int + vn.int[final-1]-vn.int[final-3] - my.i2 - my.i1
        
        for (my.i in (v.critE[final] - my.i1 - my.i2 - k.int):(vn.int[final]-vn.int[final-1])) {
          pp <- pp +
            (my.i2 >= 0)*
            my.dbetabinom(x=my.i2,size=vn.int[final-2]-vn.int[final-3],a=a2,b=b2)*
            (my.i1 >= 0)*
            my.dbetabinom(x=my.i1,size=vn.int[final-1]-vn.int[final-2],a=a1,b=b1)*
            (my.i >= 0)*
            my.dbetabinom(x=my.i,size=vn.int[final]-vn.int[final-1],a=a,b=b)
          }
        }
      }
    }
  }
  
  if (current.interim == final-4){
    pp=0
    #my.i3 for (final-4 to) final-3
    #my.i2 for (final-3 to) final-2
    #my.i1 for (final-2 to) final-1
    #my.i for (final-1 to) final
    if (k.int > v.critF[final-4]){
    a3=shape1E+k.int
    b3=shape2E+n.int-k.int
    for (my.i3 in (v.critF[final-3] + 1 - k.int):(vn.int[final-3]-vn.int[final-4])) {
      a2=shape1E+k.int + my.i3
      b2=shape2E+n.int-k.int + vn.int[final-3]-vn.int[final-4] - my.i3
    
    for (my.i2 in (v.critF[final-2] + 1 -my.i3 - k.int):(vn.int[final-2]-vn.int[final-3])) {
      a1=shape1E+k.int + my.i3 + my.i2
      b1=shape2E+n.int-k.int + vn.int[final-2]-vn.int[final-4] - my.i3 - my.i2
      
      for (my.i1 in (v.critF[final-1] + 1 - my.i3 - my.i2 - k.int):(vn.int[final-1]-vn.int[final-2])) {
        a=shape1E+k.int + my.i3 + my.i2 + my.i1 
        b=shape2E+n.int-k.int + vn.int[final-1]-vn.int[final-4] - my.i3 - my.i2 - my.i1
        
        for (my.i in (v.critE[final] - my.i1 - my.i2 - my.i3 - k.int):(vn.int[final]-vn.int[final-1])) {
          pp <- pp +
            (my.i3 >= 0)*
            my.dbetabinom(x=my.i3,size=vn.int[final-3]-vn.int[final-4],a=a3,b=b3)*
            (my.i2 >= 0)*
            my.dbetabinom(x=my.i2,size=vn.int[final-2]-vn.int[final-3],a=a2,b=b2)*
            (my.i1 >= 0)*
            my.dbetabinom(x=my.i1,size=vn.int[final-1]-vn.int[final-2],a=a1,b=b1)*
            (my.i >= 0)*
            my.dbetabinom(x=my.i,size=vn.int[final]-vn.int[final-1],a=a,b=b)
            }
          }
        }
      }
    }
  }
  
  if (current.interim == final-5){
    pp=0
    #my.i4 for (final-5 to) final-4
    #my.i3 for (final-4 to) final-3
    #my.i2 for (final-3 to) final-2
    #my.i1 for (final-2 to) final-1
    #my.i for (final-1 to) final
    if (k.int > v.critF[final-5]){
      a4=shape1E+k.int
      b4=shape2E+n.int-k.int
      for (my.i4 in (v.critF[final-4] + 1 - k.int):(vn.int[final-4]-vn.int[final-5])) {
        a3=shape1E+k.int + my.i4
        b3=shape2E+n.int-k.int + vn.int[final-4]-vn.int[final-5] - my.i4
        
        for (my.i3 in (v.critF[final-3] + 1 - my.i4 - k.int):(vn.int[final-3]-vn.int[final-4])) {
          a2=shape1E+k.int + my.i4 + my.i3
          b2=shape2E+n.int-k.int + vn.int[final-3]-vn.int[final-5] - my.i4 - my.i3
        
          for (my.i2 in (v.critF[final-2] + 1 -my.i4 - my.i3 - k.int):(vn.int[final-2]-vn.int[final-3])) {
            a1=shape1E+k.int + my.i4 + my.i3 + my.i2
            b1=shape2E+n.int-k.int + vn.int[final-2]-vn.int[final-5] - my.i4 - my.i3 - my.i2
          
            for (my.i1 in (v.critF[final-1] + 1 - my.i4 - my.i3 - my.i2 - k.int):(vn.int[final-1]-vn.int[final-2])) {
              a=shape1E+k.int + my.i4 + my.i3 + my.i2 + my.i1 
              b=shape2E+n.int-k.int + vn.int[final-1]-vn.int[final-5] - my.i4 - my.i3 - my.i2 - my.i1
            
              for (my.i in (v.critE[final] - my.i1 - my.i2 - my.i3 - my.i4 - k.int):(vn.int[final]-vn.int[final-1])) {
                pp <- pp +
                  (my.i4 >= 0)*
                  my.dbetabinom(x=my.i4,size=vn.int[final-4]-vn.int[final-5],a=a4,b=b4)*
                  (my.i3 >= 0)*
                  my.dbetabinom(x=my.i3,size=vn.int[final-3]-vn.int[final-4],a=a3,b=b3)*
                  (my.i2 >= 0)*
                  my.dbetabinom(x=my.i2,size=vn.int[final-2]-vn.int[final-3],a=a2,b=b2)*
                  (my.i1 >= 0)*
                  my.dbetabinom(x=my.i1,size=vn.int[final-1]-vn.int[final-2],a=a1,b=b1)*
                  (my.i >= 0)*
                  my.dbetabinom(x=my.i,size=vn.int[final]-vn.int[final-1],a=a,b=b)
              }
            }
          }
        }
      }
    }
  }
  
  if (current.interim == final-6){
    pp=0
    #my.i5 for (final-6 to) final-5
    #my.i4 for (final-5 to) final-4
    #my.i3 for (final-4 to) final-3
    #my.i2 for (final-3 to) final-2
    #my.i1 for (final-2 to) final-1
    #my.i for (final-1 to) final
    if (k.int > v.critF[final-6]){
      a5=shape1E+k.int
      b5=shape2E+n.int-k.int
      for (my.i5 in (v.critF[final-5] + 1 - k.int):(vn.int[final-5]-vn.int[final-6])) {
        a4=shape1E+k.int + my.i5
        b4=shape2E+n.int-k.int + vn.int[final-5]-vn.int[final-6] - my.i5
        
        for (my.i4 in (v.critF[final-4] + 1 - my.i5 - k.int):(vn.int[final-4]-vn.int[final-5])) {
          a3=shape1E+k.int + my.i5 + my.i4
          b3=shape2E+n.int-k.int + vn.int[final-4]-vn.int[final-6] - my.i5 - my.i4
        
          for (my.i3 in (v.critF[final-3] + 1 - my.i5 - my.i4 - k.int):(vn.int[final-3]-vn.int[final-4])) {
            a2=shape1E+k.int + my.i5 + my.i4 + my.i3
            b2=shape2E+n.int-k.int + vn.int[final-3]-vn.int[final-6] - my.i5 - my.i4 - my.i3
          
            for (my.i2 in (v.critF[final-2] + 1 - my.i5 - my.i4 - my.i3 - k.int):(vn.int[final-2]-vn.int[final-3])) {
              a1=shape1E+k.int + my.i5 + my.i4 + my.i3 + my.i2
              b1=shape2E+n.int-k.int + vn.int[final-2]-vn.int[final-6] - my.i5 - my.i4 - my.i3 - my.i2
            
              for (my.i1 in (v.critF[final-1] + 1 - my.i5 - my.i4 - my.i3 - my.i2 - k.int):(vn.int[final-1]-vn.int[final-2])) {
                a=shape1E+k.int + my.i5 + my.i4 + my.i3 + my.i2 + my.i1 
                b=shape2E+n.int-k.int + vn.int[final-1]-vn.int[final-6] - my.i5 - my.i4 - my.i3 - my.i2 - my.i1
              
                for (my.i in (v.critE[final] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - k.int):(vn.int[final]-vn.int[final-1])) {
                  pp <- pp +
                    (my.i5 >= 0)*
                    my.dbetabinom(x=my.i5,size=vn.int[final-5]-vn.int[final-6],a=a5,b=b5)*
                    (my.i4 >= 0)*
                    my.dbetabinom(x=my.i4,size=vn.int[final-4]-vn.int[final-5],a=a4,b=b4)*
                    (my.i3 >= 0)*
                    my.dbetabinom(x=my.i3,size=vn.int[final-3]-vn.int[final-4],a=a3,b=b3)*
                    (my.i2 >= 0)*
                    my.dbetabinom(x=my.i2,size=vn.int[final-2]-vn.int[final-3],a=a2,b=b2)*
                    (my.i1 >= 0)*
                    my.dbetabinom(x=my.i1,size=vn.int[final-1]-vn.int[final-2],a=a1,b=b1)*
                    (my.i >= 0)*
                    my.dbetabinom(x=my.i,size=vn.int[final]-vn.int[final-1],a=a,b=b)
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  pp_withoutFutility=PredPower(n=vn.int[final], n.int=n.int,k.int=k.int, k.crit=v.critE[final], shape1=shape1E, shape2=shape2E)
  

  
  
  


predictivePower=data.frame(PredPow_wFut=pp,PredPow_woFut=pp_withoutFutility)
predictivePower
}




