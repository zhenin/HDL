#
# R wrapper for the banded multiplic function 
# pout = plink output in vector form
# M = number of markers
# nb = bandwidth above the diagonal
#
# input: x(M)
# output: out(M) = R%*% x
Rbmult = function(M,nb,pout,x)
{
## 1. extend plink output to length=(M-1)*nb
##    with padded zeroes at the tailend
## 2. extended x with dummy zeros
## tailmat is read colwise
tailmat = matrix(0,nb,nb-1)
pick = lower.tri(tailmat)[nb:1,]
tailmat[pick] =pout[((M-nb)*nb+1): length(pout)]
# ttailmat = t(tailmat)  ## to check= original row-wise form
#
plkvec = c(pout[1:((M-nb)*nb)], c(tailmat))
x0 = c(x,rep(0,nb))
out = rep(0,M)
call1 = .Fortran("bmult", 
            as.integer(M),
            as.integer(nb),
            as.double(plkvec),
            as.double(x0),
            out=as.double(out))
 return(call1$out)
}

## for eigs(.) function in RSpectra
Rxfun= function(x,args){
  M = length(x)
  nb = args$nb
  pout= args$pout
  out=Rbmult(M,nb,pout,x)
  return(out)
}

## computing LDscore from plink output
Rldscore = function(M,nb,pout)
{
## similar to Rbmult(.)
tailmat = matrix(0,nb,nb-1)
  pick = lower.tri(tailmat)[nb:1,]
tailmat[pick] =pout[((M-nb)*nb+1): length(pout)]
#
plkvec = c(pout[1:((M-nb)*nb)], c(tailmat))
score = rep(0,M)
call1 = .Fortran("ldscore", 
            as.integer(M),
            as.integer(nb),
            as.double(plkvec),
            score=as.double(score))
 return(call1$score)
}

