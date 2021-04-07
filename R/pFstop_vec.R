pFstop_vec=function(n, p, pF, vn.int, cFvec, shape1=1,shape2=1, printit=FALSE){

  out=lapply(cFvec, function(cF) {
    v.crit <-crit(n=n, pThreshold=pF, vn.int=vn.int, crit=cF, shape1 = shape1, shape2 = shape2)
    res=    pFstop(p,c(vn.int,n),v.crit)
    res
  }
  )
  names(out)=cFvec

  summ=sapply(out, function(x) x$P.stop.cum)
  rownames(summ)= paste0("Pat.",c(vn.int,n))
  if (printit==TRUE) print(summ)
  res=list(summary=summ, all=out)
  invisible(res)
}


