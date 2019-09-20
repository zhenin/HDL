llfun.rg.full.likelihood <-
function(param, M, N1, N2, N0, Nref, 
                                    lam1, lam2, bstar1, bstar2, lim=exp(-10)){
  h11 = param[1]
  int.h11 = param[2]
  h22 = param[3]
  int.h22 = param[4]
  rg = param[5]
  int.h12 = param[6]
  ## sample fractions
  p1 = N0/N1; p2= N0/N2
  ## must follow the formula for lamh2 used in llfun
  lam11 =   h11/M*lam1^2 + int.h11*lam1/N1
  lam11 = ifelse(lam11<lim, lim,lam11)
  lam22 = h22/M*lam2^2 + int.h22*lam2/N2
  lam22 = ifelse(lam22<lim, lim,lam22)
  #lam12 = h12/M*lam1*lam2 - p1*p2*h12*lam1/Nref + p1*p2*int*lam1/N0
  if (N0>0) lam12 = rg*sqrt(h11*h22)/M*lam1*lam2 + p1*p2*int.h12*lam1/N0  ## key change here
  if (N0==0) lam12 = rg*sqrt(h11*h22)/M*lam1*lam2
  ## Full likelihood
  lam11.inv <- 1/lam11+(lam12/lam11)^2/(lam22-lam12^2/lam11)
  lam12.inv <- -lam12/lam11/(lam22-lam12^2/lam11)
  lam22.inv <- 1/(lam22-lam12^2/lam11)
  if(any(lam11*lam22-lam12^2<0))
    ll = 1000000
  else
    ll = sum(log(lam11*lam22-lam12^2)) + sum(bstar1^2*lam11.inv) + 2*sum(bstar1*bstar2*lam12.inv) + sum(bstar2^2*lam22.inv)
  return(ll)
}
