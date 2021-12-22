post_sas = function(response, activity, N, ninter, p0, mu1_h0, a0, mu2_h0, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  p1 = 0.25 + 0.5 * as.numeric(rowMeans(response) > p0 + 0.15)
  p2 = 0.25 + 0.5 * as.numeric(rowMeans(activity) > a0 + 0.15)
  dat = list(response = response,
             activity = activity,
             N = N,
             ninter = ninter,
             p1 = p1,
             p2 = p2,
             mu1_h0 = mu1_h0, 
             mu2_h0 = mu2_h0)
  Z = array(0, dim = c(N, ninter, 2))
  for (i in 1:N) {
    Z[i, , ] = mvrnorm(ninter, c(mu1_h0[i], mu2_h0[i]), matrix(c(1, 0.5, 0.5, 1), 2, 2))
    for (j in 1:ninter) {
      Z[i, j, 1] = ifelse(dat$response[i, j] == 1, abs(Z[i, j, 1]), -abs(Z[i, j, 1]))
      Z[i, j, 2] = ifelse(dat$activity[i, j] == 1, abs(Z[i, j, 2]), -abs(Z[i, j, 2]))
    }
  }
  ss10 = as.numeric(rowMeans(response) > p0 + 0.2)
  ss20 = as.numeric(rowMeans(activity) > a0 + 0.2)
  thismodel = try(jags.model(file = "bugs/match.txt", 
                             data = dat, 
                             inits = list(Z = Z,
                                          diff1 = 1,
                                          diff2 = 1,
                                          ss1 = ss10,
                                          ss2 = ss20,
                                          rho = rep(0.5, dat$N),
                                          tau1 = 0.001,
                                          tau2 = 0.001),
                             n.adapt = n.adapt, quiet = TRUE), silent = TRUE)
  try(update(thismodel, n.burn, progress.bar = "none"), silent = TRUE)
  res.bugs = try(jags.samples(thismodel, 
                              variable.names = c("mu1", "mu2", "rho"),
                              n.iter = n.iter, progress.bar = "none"), silent = TRUE)
  mu1_rec = matrix(res.bugs$mu1, N, n.iter)
  mu2_rec = matrix(res.bugs$mu2, N, n.iter)
  rho = matrix(res.bugs$rho, N, n.iter)
  return (list("mu1_rec" = mu1_rec, "mu2_rec" = mu2_rec, "rho" = rho))
}
