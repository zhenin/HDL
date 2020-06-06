llfun.gcov.part.2 = function(param,h11,h22,rho12, M, N1, N2, N0, Nref,
                             lam0, lam1, lam2, bstar1, bstar2, lim=exp(-10)){
  h12 = param[1]
  int = param[2]
  ## sample fractions
  p1 = N0/N1; p2= N0/N2
  ## must follow the formula for lamh2 used in llfun4
  lam11 =   h11[1]/M*lam1^2 - h11[1]*lam1/Nref + h11[2]*lam1/N1
  lam11 = ifelse(lam11<lim, lim,lam11)
  lam22 = h22[1]/M*lam2^2 - h22[1]*lam2/Nref + h22[2]*lam2/N2
  lam22 = ifelse(lam22<lim, lim,lam22)
  #lam12 = h12/M*lam1*lam2 - p1*p2*h12*lam1/Nref + p1*p2*int*lam1/N0
  if (N0>0) lam12 = h12/M*lam1*lam2 + p1*p2*int*lam1/N0  ## key change here
  if (N0==0) lam12 = h12/M*lam1*lam2
  ##  resid of bstar2 ~bstar1
  ustar = bstar2 - lam12/lam11*bstar1  ## see note
  lam22.1 = lam22 - lam12^2/lam11
  lam22.1 = ifelse(lam22.1<lim, lim,lam22.1)
  ll = sum(log(lam22.1)) + sum(ustar^2/(lam22.1))
  return(ll)
}