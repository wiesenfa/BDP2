pFstopEcall=function(p,vn.int,v.critE,v.critF){
  # Futility stopping probability at every interim: Evaluate at vn.int[i] and stop if #successes <= v.critF[i]
  # Probability of Calling efficacy at every interim: Evaluate at vn.int[i] and call if #successes >= v.critE[i]

  # Args:
  # p:           true response rate
  # vn.int:      vector of sample sizes at the interim analyses
  # v.critE:     vector of critical boundaries for calling efficacy at the interim analyses
  # v.critF:     vector of critical boundaries for futility stopping at the interim analyses

  n.interim=length(vn.int)
  if (n.interim>10) stop("max 10 interims including final")
  probcall <- rep(0,n.interim)
  probstop <- rep(0,n.interim)

  probstop[1] <- pbinom(size=vn.int[1],q=v.critF[1] ,p, lower.tail=TRUE)

  dbin.i1=dbinom(size=vn.int[1],x=(v.critF[1] + 1):vn.int[1],p)

  if (n.interim>1){
    i1=0
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      i1=i1+1
      dbin.i2=dbinom(size=vn.int[2]-vn.int[1],x=0:(v.critF[2] - my.i1),p)
      i2=0
      for (my.i2 in 0:(v.critF[2] - my.i1)) {
        i2=i2+1
        probstop[2] <- probstop[2] +
          dbin.i1[i1]*#dbinom(size=vn.int[1],x=my.i1,p)* #
          (v.critF[2] - my.i1 >= 0)*dbin.i2[i2]#dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)#
      }
    }
  }
  if (n.interim>2){
    i1=0
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      i1=i1+1
      dbin.i2=dbinom(size=vn.int[2]-vn.int[1],x=(v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] ),p)
      i2=0
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        i2=i2+1
        dbin.i3=dbinom(size=vn.int[3]-vn.int[2],x= 0:(v.critF[3] - my.i1 - my.i2 ),p)
        i3=0
        for (my.i3 in 0:(v.critF[3] - my.i1 - my.i2 )) {
          i3=i3+1
          probstop[3] <- probstop[3] +
            dbin.i1[i1]*# dbinom(size=vn.int[1],x=my.i1,p)*
            (my.i2 >= 0)*dbin.i2[i2]*#dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (v.critF[3] - my.i1 - my.i2  >= 0) *dbin.i3[i3]#* dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)

        }
      }
    }
  }
  if (n.interim>3){

    i1=0
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      i1=i1+1
      dbin.i2=dbinom(size=vn.int[2]-vn.int[1],x= (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] ),p)
      i2=0
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        i2=i2+1
        dbin.i3=dbinom(size=vn.int[3]-vn.int[2],x= (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] ),p)
        i3=0
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          i3=i3+1
          dbin.i4=dbinom(size=vn.int[4]-vn.int[3],x=  0:(v.critF[4]- my.i1 - my.i2 - my.i3 ),p)
          i4=0
          for (my.i4 in 0:(v.critF[4]- my.i1 - my.i2 - my.i3 )) {
            i4=i4+1
            probstop[4] <- probstop[4] +
              dbin.i1[i1]*#dbinom(size=vn.int[1],x=my.i1,p)*
              (my.i2 >= 0)*dbin.i2[i2]*#dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
              (my.i3 >= 0)*dbin.i3[i3]*#dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
              (v.critF[4] - my.i3 - my.i2 - my.i1 >= 0)*dbin.i4[i4]#dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p)
          }
        }
      }
    }
}
  if (n.interim>4){
    i1=0
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      i1=i1+1
      dbin.i2=dbinom(size=vn.int[2]-vn.int[1],x= (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] ),p)
      i2=0
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        i2=i2+1
        dbin.i3=dbinom(size=vn.int[3]-vn.int[2],x= (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] ),p)
        i3=0
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          i3=i3+1
          dbin.i4=dbinom(size=vn.int[4]-vn.int[3],x=   (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] ),p)
          i4=0
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            i4=i4+1
            dbin.i5=dbinom(size=vn.int[5]-vn.int[4],x=0:(v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4) ,p)
            i5=0
            for (my.i5 in 0:(v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4)) {
              i5=i5+1
              probstop[5] <- probstop[5] +
                dbin.i1[i1]*#dbinom(size=vn.int[1],x=my.i1,p)*
                (my.i2 >= 0)*dbin.i2[i2]*#dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                (my.i3 >= 0)*dbin.i3[i3]*#dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                (my.i4 >= 0)*dbin.i4[i4]*
                (v.critF[5]- my.i4 -my.i3 - my.i2 - my.i1 >= 0)*dbin.i5[i5]
          }
        }
      }
      }
    }


    # for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    #   for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
    #     for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
    #       for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
    #         for (my.i5 in 0:(v.critF[5] - my.i1 - my.i2 - my.i3 - my.i4)) {
    #
    #           probstop[5] <- probstop[5] +
    #
    #             dbinom(size=vn.int[1],x=my.i1,p)*
    #             (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
    #             (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
    #             (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
    #             (v.critF[5]- my.i4 -my.i3 - my.i2 - my.i1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p)
    #         }
    #       }
    #     }
    #   }
    # }
}
  if (n.interim>5){
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in 0:(v.critF[6] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5)) {

                probstop[6] <- probstop[6] +

                  dbinom(size=vn.int[1],x=my.i1,p)*
                  (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                  (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                  (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                  (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                  (v.critF[6]- my.i5 - my.i4 -my.i3 - my.i2 - my.i1 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p)
              }
            }
          }
        }
      }
    }
  }

  if (n.interim>6){
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in (v.critF[6]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(vn.int[6]-vn.int[5] )) {
                for (my.i7 in 0:(v.critF[7] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6)) {

                  probstop[7] <- probstop[7] +

                    dbinom(size=vn.int[1],x=my.i1,p)*
                    (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                    (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                    (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                    (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                    (my.i6 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p) *
                    (v.critF[7]- my.i6 - my.i5 - my.i4 -my.i3 - my.i2 - my.i1 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7 ,p)
                }
              }
            }
          }
        }
      }
    }
  }
  if (n.interim>7){
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in (v.critF[6]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(vn.int[6]-vn.int[5] )) {
                for (my.i7 in (v.critF[7]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6):(vn.int[7]-vn.int[6] )) {
                  for (my.i8 in 0:(v.critF[8] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6-my.i7)) {

                    probstop[8] <- probstop[8] +

                      dbinom(size=vn.int[1],x=my.i1,p)*
                      (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2, p) *
                      (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                      (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                      (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                      (my.i6 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p) *
                      (my.i7 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7 ,p) *
                      (v.critF[8]- my.i7 -my.i6 - my.i5 - my.i4 -my.i3 - my.i2 - my.i1 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8 ,p)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (n.interim>8){
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in (v.critF[6]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(vn.int[6]-vn.int[5] )) {
                for (my.i7 in (v.critF[7]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6):(vn.int[7]-vn.int[6] )) {
                  for (my.i8 in (v.critF[8]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6-my.i7):(vn.int[8]-vn.int[7] )) {
                    for (my.i9 in 0:(v.critF[9] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6-my.i7-my.i8)) {

                      probstop[9] <- probstop[9] +

                        dbinom(size=vn.int[1],x=my.i1,p)*
                        (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                        (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                        (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                        (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                        (my.i6 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p) *
                        (my.i7 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7 ,p) *
                        (my.i8 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8 ,p) *
                        (v.critF[9]- my.i8 - my.i7 -my.i6 - my.i5 - my.i4 -my.i3 - my.i2 - my.i1 >= 0)*dbinom(size=vn.int[9]-vn.int[8],x=my.i9
                                                                                                            ,p)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (n.interim>9){
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in (v.critF[6]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(vn.int[6]-vn.int[5] )) {
                for (my.i7 in (v.critF[7]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6):(vn.int[7]-vn.int[6] )) {
                  for (my.i8 in (v.critF[8]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7):(vn.int[8]-vn.int[7] )) {
                    for (my.i9 in (v.critF[9]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8):(vn.int[9]-vn.int[8] )) {
                      for (my.i10 in 0:(v.critF[10] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8 - my.i9)) {

                        probstop[10] <- probstop[10] +

                          dbinom(size=vn.int[1],x=my.i1,p)*
                          (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                          (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                          (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                          (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                          (my.i6 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p) *
                          (my.i7 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7 ,p) *
                          (my.i8 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8 ,p) *
                          (my.i9 >= 0)*dbinom(size=vn.int[9]-vn.int[8],x=my.i9 ,p) *
                          (v.critF[10] - my.i9 - my.i8 - my.i7 - my.i6 - my.i5 - my.i4 -my.i3 - my.i2 - my.i1 >= 0)*dbinom(size=vn.int[10]-vn.int
                                                                                                                      [9],x=my.i10 ,p)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
#  print("eff")
  probcall[1] <- pbinom(size=vn.int[1],q=v.critE[1]-1 ,p, lower.tail=FALSE)

  dbin.i1=dbinom(size=vn.int[1],x= (v.critF[1] + 1):vn.int[1],p)

  if (n.interim>1){
    i1=0
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      i1=i1+1
      dbin.i2=dbinom(size=vn.int[2]-vn.int[1],x= (v.critE[2] - my.i1):(vn.int[2]-vn.int[1]),p)
      i2=0
      for (my.i2 in (v.critE[2] - my.i1):(vn.int[2]-vn.int[1])) {
        i2=i2+1
        probcall[2] <- probcall[2] +
                dbin.i1[i1]*#dbinom(size=vn.int[1],x=my.i1,p)*
                (my.i2 >= 0)*dbin.i2[i2]#dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
       }
    }

    # for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    #   for (my.i2 in (v.critE[2] - my.i1):(vn.int[2]-vn.int[1])) {
    #     probcall[2] <- probcall[2] +
    #       dbinom(size=vn.int[1],x=my.i1,p)*
    #       (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)
    #   }
    # }

  }

  if (n.interim>2){
    i1=0
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      i1=i1+1
      dbin.i2=dbinom(size=vn.int[2]-vn.int[1],x= (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] ),p)
      i2=0
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        i2=i2+1
        dbin.i3=dbinom(size=vn.int[3]-vn.int[2],x= (v.critE[3] - my.i1 - my.i2):(vn.int[3]-vn.int[2]),p)
        i3=0
        for (my.i3 in (v.critE[3] - my.i1 - my.i2):(vn.int[3]-vn.int[2])) {
          i3=i3+1
          probcall[3] <- probcall[3] +
                dbin.i1[i1]*#dbinom(size=vn.int[1],x=my.i1,p)*
                (my.i2 >= 0)*dbin.i2[i2]*#dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                (my.i3 >= 0)*dbin.i3[i3]#dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
       }
      }
    }



    #  for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    #   for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
    #     for (my.i3 in (v.critE[3] - my.i1 - my.i2):(vn.int[3]-vn.int[2])) {
    #       probcall[3] <- probcall[3] +
    #
    #         dbinom(size=vn.int[1],x=my.i1,p)*
    #         (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
    #         (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)
    #
    #     }
    #   }
    # }

  }

  if (n.interim>3){
    i1=0
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      i1=i1+1
      dbin.i2=dbinom(size=vn.int[2]-vn.int[1],x= (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] ),p)
      i2=0
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        i2=i2+1
        dbin.i3=dbinom(size=vn.int[3]-vn.int[2],x= (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] ),p)
        i3=0
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          i3=i3+1
          dbin.i4=dbinom(size=vn.int[4]-vn.int[3],x=    (v.critE[4] - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] ),p)
          i4=0
          for (my.i4 in (v.critE[4] - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            i4=i4+1
            probcall[4] <- probcall[4] +
                dbin.i1[i1]*#dbinom(size=vn.int[1],x=my.i1,p)*
                (my.i2 >= 0)*dbin.i2[i2]*#dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                (my.i3 >= 0)*dbin.i3[i3]*#dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                (my.i4 >= 0)*dbin.i4[i4]#dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p)

        }
      }
      }
    }






 #    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
 #      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
 #        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
 #          for (my.i4 in (v.critE[4] - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
 #            probcall[4] <- probcall[4] +
 #
 #              dbinom(size=vn.int[1],x=my.i1,p)*
 #              (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p)*
 #              (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p)*
 #              (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p)
 #
 #          }
 #        }
 #      }
 #    }
  }

  if (n.interim>4){
    i1=0
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      i1=i1+1
      dbin.i2=dbinom(size=vn.int[2]-vn.int[1],x= (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] ),p)
      i2=0
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        i2=i2+1
        dbin.i3=dbinom(size=vn.int[3]-vn.int[2],x= (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] ),p)
        i3=0
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          i3=i3+1
          dbin.i4=dbinom(size=vn.int[4]-vn.int[3],x=   (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] ),p)
          i4=0
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            i4=i4+1
            dbin.i5=dbinom(size=vn.int[5]-vn.int[4],x=(v.critE[5] - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] ) ,p)
            i5=0
            for (my.i5 in (v.critE[5] - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              i5=i5+1
              probcall[5] <- probcall[5] +
                dbin.i1[i1]*#dbinom(size=vn.int[1],x=my.i1,p)*
                (my.i2 >= 0)*dbin.i2[i2]*#dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                (my.i3 >= 0)*dbin.i3[i3]*#dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                (my.i4 >= 0)*dbin.i4[i4]*#dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p)
                (my.i5 >= 0)*dbin.i5[i5]
          }
        }
      }
      }
    }








    # for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
    #   for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
    #     for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
    #       for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
    #         for (my.i5 in (v.critE[5] - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
    #
    #           probcall[5] <- probcall[5] +
    #
    #             dbinom(size=vn.int[1],x=my.i1,p)*
    #             (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p)*
    #             (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p)*
    #             (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p)*
    #             (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p)
    #         }
    #       }
    #     }
    #   }
    # }

  }

  if (n.interim>5)
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in (v.critE[6] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(vn.int[6]-vn.int[5] )) {

                probcall[6] <- probcall[6] +

                  dbinom(size=vn.int[1],x=my.i1,p)*
                  (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                  (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                  (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                  (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                  (my.i6 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p)
              }
            }
          }
        }
      }
    }

  if (n.interim>6)
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in (v.critF[6]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(vn.int[6]-vn.int[5] )) {
                for (my.i7 in (v.critE[7] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6):(vn.int[7]-vn.int[6] )) {

                  probcall[7] <- probcall[7] +

                    dbinom(size=vn.int[1],x=my.i1,p)*
                    (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                    (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                    (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                    (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                    (my.i6 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p) *
                    (my.i7 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7 ,p)
                }
              }
            }
          }
        }
      }
    }

  if (n.interim>7)
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in (v.critF[6]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(vn.int[6]-vn.int[5] )) {
                for (my.i7 in (v.critF[7]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6):(vn.int[7]-vn.int[6] )) {
                  for (my.i8 in (v.critE[8] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7):(vn.int[8]-vn.int[7] )) {

                    probcall[8] <- probcall[8] +

                      dbinom(size=vn.int[1],x=my.i1,p)*
                      (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                      (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                      (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                      (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                      (my.i6 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p) *
                      (my.i7 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7 ,p) *
                      (my.i8 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8 ,p)
                  }
                }
              }
            }
          }
        }
      }
    }

  if (n.interim>8)
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in (v.critF[6]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(vn.int[6]-vn.int[5] )) {
                for (my.i7 in (v.critF[7]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6):(vn.int[7]-vn.int[6] )) {
                  for (my.i8 in (v.critF[8]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7):(vn.int[8]-vn.int[7] )) {
                    for (my.i9 in (v.critE[9] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8):(vn.int[9]-vn.int[8] )) {

                      probcall[9] <- probcall[9] +

                        dbinom(size=vn.int[1],x=my.i1,p)*
                        (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p) *
                        (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                        (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                        (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                        (my.i6 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p) *
                        (my.i7 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7 ,p) *
                        (my.i8 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8 ,p) *
                        (my.i9 >= 0)*dbinom(size=vn.int[9]-vn.int[8],x=my.i9 ,p)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

  if (n.interim>9)
    for (my.i1 in (v.critF[1] + 1):vn.int[1]) {
      for (my.i2 in (v.critF[2]+1 - my.i1):(vn.int[2]-vn.int[1] )) {
        for (my.i3 in (v.critF[3]+1 - my.i1 - my.i2):(vn.int[3]-vn.int[2] )) {
          for (my.i4 in (v.critF[4]+1 - my.i1 - my.i2 - my.i3):(vn.int[4]-vn.int[3] )) {
            for (my.i5 in (v.critF[5]+1 - my.i1 - my.i2 - my.i3 - my.i4):(vn.int[5]-vn.int[4] )) {
              for (my.i6 in (v.critF[6]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(vn.int[6]-vn.int[5] )) {
                for (my.i7 in (v.critF[7]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6):(vn.int[7]-vn.int[6] )) {
                  for (my.i8 in (v.critF[8]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7):(vn.int[8]-vn.int[7] )) {
                    for (my.i9 in (v.critF[9]+1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8):(vn.int[9]-vn.int[8] )) {
                      for (my.i10 in (v.critE[10] - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8 - my.i9):(vn.int[10]-vn.int[9] )) {

                        probcall[10] <- probcall[10] +

                          dbinom(size=vn.int[1],x=my.i1,p)*
                          (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2 ,p)*
                          (my.i3 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3 ,p) *
                          (my.i4 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4 ,p) *
                          (my.i5 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5 ,p) *
                          (my.i6 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6 ,p) *
                          (my.i7 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7 ,p) *
                          (my.i8 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8 ,p) *
                          (my.i9 >= 0)*dbinom(size=vn.int[9]-vn.int[8],x=my.i9 ,p) *
                          (my.i10 >= 0)*dbinom(size=vn.int[10]-vn.int[9],x=my.i10 ,p)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }


probFstopEcall=data.frame(Int.Analysis=1:length(vn.int),Pat.number=vn.int,
                            Crit.boundaryE=v.critE,
                            P.effic=round(probcall, digits = 4),
  #                          P.effic.cum=pmax(100,round(100*cumsum(probcall), digits = 2)),
                            Crit.boundaryF=v.critF,
                            P.futil=round(probstop, digits = 4),
                            P.futil.cum=round(cumsum(probstop), digits = 4) )
  probFstopEcall
}




