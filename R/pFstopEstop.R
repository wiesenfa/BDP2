pFstopEstop=function(p,vn.int,v.critE,v.critF){

#### Stopping probability at every interim: Evaluate at vn.int[i]
#### stop for futility if #successes <= v.critF[i]
#### stop for efficacy if #successes >= v.critE[i]

  n.interim=length(vn.int)
  if (n.interim>10) stop("max 10 interims including final")
  probstopF <- rep(0,n.interim)
  probstopE <- rep(0,n.interim)

  probstopF[1] <- pbinom(size=vn.int[1],q=v.critF[1] ,p, lower.tail=TRUE)

  if (n.interim>1)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in 0:(v.critF[2] - my.i1)) {
        probstopF[2] <- probstopF[2] +
          (v.critE[1] - 1 - v.critF[1] - 1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
          (v.critF[2] - my.i1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)

      }
    }

  if (n.interim>2)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in 0:(v.critF[3] - my.i1 - my.i2 )) {
        probstopF[3] <- probstopF[3] +
          (v.critE[1] - 1 - v.critF[1] - 1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
          (v.critE[2] - 1 - v.critF[2] - 1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
          (v.critF[3] - my.i1 - my.i2  >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)
        }
      }
    }

  if (n.interim>3)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1 -my.i1 - my.i2):(v.critE[3] - 1-my.i1 - my.i2)) {
          for (my.i4 in 0:(v.critF[4] - my.i1 - my.i2 - my.i3)) {
          probstopF[4] <- probstopF[4] +
            (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
            (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
            (v.critF[4] - my.i1 - my.i2 - my.i3 >= 0) * dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)
          }
        }
      }
    }


  if (n.interim>4)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 -my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1-my.i1 -my.i2):(v.critE[3] - 1-my.i1 -my.i2)) {
          for (my.i4 in (v.critF[4] + 1-my.i1 -my.i2 - my.i3):(v.critE[4] - 1-my.i1 - my.i2 - my.i3)) {
            for (my.i5 in 0:(v.critF[5] - my.i1 - my.i2 -my.i3 - my.i4)) {
            probstopF[5] <- probstopF[5] +
              (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
              (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
              (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
              (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
              (v.critF[5] - my.i1 - my.i2 -my.i3  - my.i4>= 0) * dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)
            }
          }
        }
      }
    }

  if (n.interim>5)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in 0:(v.critF[6] - my.i1 - my.i2 -my.i3 - my.i4 - my.i5)) {
                probstopF[6] <- probstopF[6] +
                  (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                  (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                  (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                  (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                  (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                  (v.critF[6] - my.i1 - my.i2 - my.i3  - my.i4 - my.i5>= 0) * dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)
              }
            }
          }
        }
      }
    }

  if (n.interim>6)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in (v.critF[6] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(v.critE[6] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5)) {
                for (my.i7 in 0:(v.critF[7] - my.i1 - my.i2 -my.i3 - my.i4 - my.i5 - my.i6)) {
                probstopF[7] <- probstopF[7] +
                  (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                  (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                  (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                  (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                  (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                  (v.critE[6] - 1 -v.critF[6] -1 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)*
                  (v.critF[7] - my.i1 - my.i2 - my.i3  - my.i4 - my.i5 - my.i6>= 0) * dbinom(size=vn.int[7]-vn.int[6],x=my.i7,p)
                }
              }
            }
          }
        }
      }
    }

  if (n.interim>7)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in (v.critF[6] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(v.critE[6] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5)) {
                for (my.i7 in (v.critF[7] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 -my.i6):(v.critE[7] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6)) {
                  for (my.i8 in 0:(v.critF[8] - my.i1 - my.i2 -my.i3 - my.i4 - my.i5 - my.i6 - my.i7)) {
                  probstopF[8] <- probstopF[8] +
                    (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                    (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                    (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                    (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                    (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                    (v.critE[6] - 1 -v.critF[6] -1 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)*
                    (v.critE[7] - 1 -v.critF[7] -1 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7,p)*
                    (v.critF[8] - my.i1 - my.i2 - my.i3  - my.i4 - my.i5 - my.i6 -my.i7>= 0) * dbinom(size=vn.int[8]-vn.int[7],x=my.i8,p)
                  }
                }
              }
            }
          }
        }
      }
    }

  if (n.interim>8)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in (v.critF[6] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(v.critE[6] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5)) {
                for (my.i7 in (v.critF[7] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 -my.i6):(v.critE[7] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6)) {
                  for (my.i8 in (v.critF[8] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7):(v.critE[8] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7)) {
                    for (my.i9 in 0:(v.critF[9] - my.i1 - my.i2 -my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8)) {
                    probstopF[9] <- probstopF[9] +
                      (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                      (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                      (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                      (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                      (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                      (v.critE[6] - 1 -v.critF[6] -1 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)*
                      (v.critE[7] - 1 -v.critF[7] -1 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7,p)*
                      (v.critE[8] - 1 -v.critF[8] -1 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8,p)*
                      (v.critF[9] - my.i1 - my.i2 - my.i3  - my.i4 - my.i5 - my.i6 - my.i7 - my.i8>= 0) * dbinom(size=vn.int[9]-vn.int[8],x=my.i9,p)
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
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in (v.critF[6] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(v.critE[6] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5)) {
                for (my.i7 in (v.critF[7] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 -my.i6):(v.critE[7] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6)) {
                  for (my.i8 in (v.critF[8] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7):(v.critE[8] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7)) {
                    for (my.i9 in (v.critF[9] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8):(v.critE[9] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8)) {
                      for (my.i10 in 0:(v.critF[10] - my.i1 - my.i2 -my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8 - my.i9)) {
                      probstopF[10] <- probstopF[10] +
                        (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                        (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                        (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                        (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                        (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                        (v.critE[6] - 1 -v.critF[6] -1 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)*
                        (v.critE[7] - 1 -v.critF[7] -1 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7,p)*
                        (v.critE[8] - 1 -v.critF[8] -1 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8,p)*
                        (v.critE[9] - 1 -v.critF[9] -1 >= 0)*dbinom(size=vn.int[9]-vn.int[8],x=my.i9,p)*
                        (v.critF[10] - my.i1 - my.i2 - my.i3  - my.i4 - my.i5 - my.i6 - my.i7 - my.i8 - my.i9>= 0) * dbinom(size=vn.int[10]-vn.int[9],x=my.i10,p)
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







  probstopE[1] <- pbinom(size=vn.int[1],q=v.critE[1]-1 ,p, lower.tail=FALSE)

  if (n.interim>1)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critE[2] -my.i1):(vn.int[2]-vn.int[1])) {
        probstopE[2] <- probstopE[2] +
          (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
          (my.i2 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)
      }
    }

  if (n.interim>2)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1- my.i1):(v.critE[2] - 1- my.i1)) {
        for (my.i3 in (v.critE[3] - my.i1 -my.i2):(vn.int[3]-vn.int[2])) {
          probstopE[3] <- probstopE[3] +
            (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
            (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
            (my.i3  >= 0) * dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)
        }
      }
    }

  if (n.interim>3)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 -my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critE[4] - my.i1 -my.i2 - my.i3):(vn.int[4]-vn.int[3])) {
            probstopE[4] <- probstopE[4] +
              (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
              (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
              (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
              (my.i4 >= 0) * dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)
          }
        }
      }
    }

  if (n.interim>4)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 -my.i1):(v.critE[2] - 1 -my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 -my.i2):(v.critE[3] - 1- my.i1 -my.i2)) {
          for (my.i4 in (v.critF[4] + 1- my.i1 -my.i2 - my.i3):(v.critE[4] - 1- my.i1 - my.i2 - my.i3)) {
            for (my.i5 in (v.critE[5] - my.i1 -my.i2 - my.i3 -my.i4):(vn.int[5]-vn.int[4])) {
              probstopE[5] <- probstopE[5] +
                (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                (my.i5>= 0) * dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)
            }
          }
        }
      }
    }

  if (n.interim>5)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in (v.critE[6] - my.i1 -my.i2 - my.i3 -my.i4 - my.i5):(vn.int[6]-vn.int[5])) {
                probstopE[6] <- probstopE[6] +
                  (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                  (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                  (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                  (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                  (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                  (my.i6>= 0) * dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)
              }
            }
          }
        }
      }
    }

  if (n.interim>6)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in (v.critF[6] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(v.critE[6] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5)) {
                for (my.i7 in (v.critE[7] - my.i1 -my.i2 - my.i3 -my.i4 - my.i5 - my.i6):(vn.int[7]-vn.int[6])) {
                  probstopE[7] <- probstopE[7] +
                    (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                    (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                    (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                    (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                    (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                    (v.critE[6] - 1 -v.critF[6] -1 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)*
                    (my.i7>= 0) * dbinom(size=vn.int[7]-vn.int[6],x=my.i7,p)
                }
              }
            }
          }
        }
      }
    }

  if (n.interim>7)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in (v.critF[6] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(v.critE[6] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5)) {
                for (my.i7 in (v.critF[7] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 -my.i6):(v.critE[7] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6)) {
                  for (my.i8 in (v.critE[8] - my.i1 -my.i2 - my.i3 -my.i4 - my.i5 - my.i6 - my.i7):(vn.int[8]-vn.int[7])) {
                    probstopE[8] <- probstopE[8] +
                      (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                      (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                      (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                      (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                      (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                      (v.critE[6] - 1 -v.critF[6] -1 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)*
                      (v.critE[7] - 1 -v.critF[7] -1 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7,p)*
                      (my.i8>= 0) * dbinom(size=vn.int[8]-vn.int[7],x=my.i8,p)
                  }
                }
              }
            }
          }
        }
      }
    }

  if (n.interim>8)
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in (v.critF[6] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(v.critE[6] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5)) {
                for (my.i7 in (v.critF[7] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 -my.i6):(v.critE[7] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6)) {
                  for (my.i8 in (v.critF[8] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7):(v.critE[8] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7)) {
                    for (my.i9 in (v.critE[9] - my.i1 -my.i2 - my.i3 -my.i4 - my.i5 - my.i6 - my.i7 - my.i8):(vn.int[9]-vn.int[8])) {
                      probstopE[9] <- probstopE[9] +
                        (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                        (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                        (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                        (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                        (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                        (v.critE[6] - 1 -v.critF[6] -1 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)*
                        (v.critE[7] - 1 -v.critF[7] -1 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7,p)*
                        (v.critE[8] - 1 -v.critF[8] -1 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8,p)*
                        (my.i9>= 0) * dbinom(size=vn.int[9]-vn.int[8],x=my.i9,p)
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
    for (my.i1 in (v.critF[1] + 1):(v.critE[1] - 1)) {
      for (my.i2 in (v.critF[2] + 1 - my.i1):(v.critE[2] - 1 - my.i1)) {
        for (my.i3 in (v.critF[3] + 1- my.i1 - my.i2):(v.critE[3] - 1 - my.i1 - my.i2)) {
          for (my.i4 in (v.critF[4] + 1 - my.i1 - my.i2 - my.i3):(v.critE[4] - 1 - my.i1 -my.i2 - my.i3)) {
            for (my.i5 in (v.critF[5] + 1 - my.i1 - my.i2 - my.i3 - my.i4):(v.critE[5] - 1 - my.i1 - my.i2 - my.i3 - my.i4)) {
              for (my.i6 in (v.critF[6] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5):(v.critE[6] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5)) {
                for (my.i7 in (v.critF[7] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 -my.i6):(v.critE[7] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6)) {
                  for (my.i8 in (v.critF[8] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7):(v.critE[8] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7)) {
                    for (my.i9 in (v.critF[9] + 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8):(v.critE[9] - 1 - my.i1 - my.i2 - my.i3 - my.i4 - my.i5 - my.i6 - my.i7 - my.i8)) {
                      for (my.i10 in (v.critE[10] - my.i1 -my.i2 - my.i3 -my.i4 - my.i5 - my.i6 - my.i7 - my.i8 - my.i9):(vn.int[10]-vn.int[9])) {
                        probstopE[10] <- probstopE[10] +
                          (v.critE[1] - 1 -v.critF[1] -1 >= 0)*dbinom(size=vn.int[1],x=my.i1,p)*
                          (v.critE[2] - 1 -v.critF[2] -1 >= 0)*dbinom(size=vn.int[2]-vn.int[1],x=my.i2,p)*
                          (v.critE[3] - 1 -v.critF[3] -1 >= 0)*dbinom(size=vn.int[3]-vn.int[2],x=my.i3,p)*
                          (v.critE[4] - 1 -v.critF[4] -1 >= 0)*dbinom(size=vn.int[4]-vn.int[3],x=my.i4,p)*
                          (v.critE[5] - 1 -v.critF[5] -1 >= 0)*dbinom(size=vn.int[5]-vn.int[4],x=my.i5,p)*
                          (v.critE[6] - 1 -v.critF[6] -1 >= 0)*dbinom(size=vn.int[6]-vn.int[5],x=my.i6,p)*
                          (v.critE[7] - 1 -v.critF[7] -1 >= 0)*dbinom(size=vn.int[7]-vn.int[6],x=my.i7,p)*
                          (v.critE[8] - 1 -v.critF[8] -1 >= 0)*dbinom(size=vn.int[8]-vn.int[7],x=my.i8,p)*
                          (v.critE[9] - 1 -v.critF[9] -1 >= 0)*dbinom(size=vn.int[9]-vn.int[8],x=my.i9,p)*
                          (my.i10>= 0) * dbinom(size=vn.int[10]-vn.int[9],x=my.i10,p)
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



  probstopFstopE=data.frame(Int.Analysis=1:length(vn.int),Pat.number=vn.int,
                            Crit.boundaryE=v.critE,P.effic=round(probstopE, digits = 4),P.effic.cum=round(cumsum(probstopE), digits = 4),
                            Crit.boundaryF=v.critF,P.futil=round(probstopF, digits = 4),P.futil.cum=round(cumsum(probstopF), digits = 4) )
  probstopFstopE
}




