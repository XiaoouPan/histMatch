hist_match = function(y, x, size, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  dat = list(y = y, x = x, size = size)
  thismodel = jags.model(file = "bugs/match.txt", 
                             data = dat, 
                             inits = list(z = y,
                                          w = rep(1, 4),
                                          u0 = rep(0, 7),
                                          u1 = rep(0, 7),
                                          u2 = rep(0, 7),
                                          u3 = rep(0, 7),
                                          u4 = rep(0, 7),
                                          tau = 0.001),
                             n.adapt = n.adapt)
  update(thismodel, n.burn)
  res.bugs = jags.samples(thismodel, 
                              variable.names = c("w", "u0", "u1", "u2", "u3", "u4", "tau"),
                              n.iter = n.iter)
  w_post = matrix(res.bugs$w, 4, n.iter)
  u0_post = matrix(res.bugs$u0, 7, n.iter)
  u1_post = matrix(res.bugs$u1, 7, n.iter)
  u2_post = matrix(res.bugs$u2, 7, n.iter)
  u3_post = matrix(res.bugs$u3, 7, n.iter)
  u4_post = matrix(res.bugs$u4, 7, n.iter)
  tau_post = as.numeric(res.bugs$tau)
  return (list("w_post" = w_post, "u0_post" = u0_post, "u1_post" = u1_post, "u2_post" = u2_post, "u3_post" = u3_post, 
               "u4_post" = u4_post, "tau_post" = tau_post))
}

mix_effect = function(y, year, omega, size, mu, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  dat = list(y = y, year = year, omega = omega, size = size, mu = mu)
  thismodel = jags.model(file = "bugs/effects.txt", 
                         data = dat, 
                         inits = list(alpha = rep(0, 5),
                                      beta = 0,
                                      effect = rep(0, 5),
                                      tau = 0.001),
                         n.adapt = n.adapt)
  update(thismodel, n.burn)
  res.bugs = jags.samples(thismodel, 
                          variable.names = c("alpha", "beta", "effect", "tau"),
                          n.iter = n.iter)
  alpha_post = matrix(res.bugs$alpha, 5, n.iter)
  beta_post = as.numeric(res.bugs$beta)
  effect_post = matrix(res.bugs$effect, 5, n.iter)
  tau_post = as.numeric(res.bugs$tau)
  return (list("alpha_post" = alpha_post, "beta_post" = beta_post, "effect_post" = effect_post, "tau_post" = tau_post))
}

naive = function(y, size, n.adapt = 5000, n.burn = 5000, n.iter = 10000) {
  dat = list(y = y, size = size)
  thismodel = jags.model(file = "bugs/naive.txt", 
                         data = dat, 
                         inits = list(alpha = 0),
                         n.adapt = n.adapt)
  update(thismodel, n.burn)
  res.bugs = jags.samples(thismodel, 
                          variable.names = c("alpha"),
                          n.iter = n.iter)
  alpha_post = as.numeric(res.bugs$alpha)
  return (alpha_post)
}
