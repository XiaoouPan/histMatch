hist_match = function(y, x, n1, n, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  dat = list(y = y, x = x, n1 = n1, n = n)
  thismodel = try(jags.model(file = "bugs/match.txt", 
                             data = dat, 
                             inits = list(z = y,
                                          w = rep(1, n -  n1),
                                          u0 = rep(0, 6),
                                          u1 = rep(0, 6)),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("w", "u0", "u1"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  w_rec = matrix(res.bugs$w, n, n.iter)
  return (list("w_rec" = w_rec))
}
