hist_match = function(y, x, n1, n, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  dat = list(y = y, x = x, n1 = n1, n = n)
  thismodel = jags.model(file = "bugs/match.txt", 
                             data = dat, 
                             inits = list(z = y,
                                          w = rep(1, n -  n1),
                                          u0 = rep(0, 7),
                                          u1 = rep(0, 7)),
                             n.adapt = n.adapt)
  update(thismodel, n.burn)
  res.bugs = jags.samples(thismodel, 
                              variable.names = c("w", "u0", "u1"),
                              n.iter = n.iter)
  w_post = matrix(res.bugs$w, n - n1, n.iter)
  u0_post = matrix(res.bugs$u0, 7, n.iter)
  u1_post = matrix(res.bugs$u1, 7, n.iter)
  return (list("w_post" = w_post, "u0_post" = u0_post, "u1_post" = u1_post))
}
