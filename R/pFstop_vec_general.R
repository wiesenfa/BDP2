pFstop_vec_general=function(n, p, p0, p1, vn.int, cFvec, shape1=1,shape2=1, type=5, alpha=0.05,printit=FALSE){

  out=lapply(cFvec, function(cF) {
    v.crit <-crit_general(n=n, p0=p0, p1=p1, vn.int=vn.int, alpha=alpha, crit=cF, type=type, shape1 = shape1, shape2 = shape2)
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


