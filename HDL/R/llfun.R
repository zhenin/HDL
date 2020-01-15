llfun <-
function(param, N, M,Nref=1000, lam, bstar, lim=exp(-10)){
  h2 = param[1]
  int= param[2]
  lamh2 = h2/M*lam^2 - h2*lam/Nref + int*lam/N
  lamh2 = ifelse(lamh2<lim, lim,lamh2)
  ll = sum(log(lamh2)) + sum(bstar^2/(lamh2))
  return(ll)
}
