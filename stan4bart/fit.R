library(stan4bart)
`%notin%` <- Negate(`%in%`)

setwd('~/Dropbox/dbarts_family/stan4bart/')
source('draw.R')
source('pull_estimates.R')

sim <- draw_dat()


m1 <- stan4bart::stan4bart(y ~ bart(. -schoolid -ID) + (1|schoolid),
                           data = sim[['worldA']][["data"]],
                           treatment = z,
                           cores = 1,
                           chains = 10,
                           iter = 4000
)


sim[['worldA']][["satt"]] <- cbind(sim[['worldA']][["satt"]], pull_estimates(fit = m1)[['satt_est']])

sim[["worldA"]][["gsatt"]] <- cbind(sim[["worldA"]][["gsatt"]],
                                    merge(sim[['worldA']][["gsatt"]], pull_estimates(fit = m1)[['gsatt_est']]))

sim[["worldA"]][["icatt"]] <- cbind(sim[["worldA"]][["icatt"]], 
                                    merge(sim[['worldA']][["icatt"]], pull_estimates(fit = m1)[['icatt_est']]))

sim[["worldA"]][['rhat']] <- rbind(sim[["worldA"]][['rhat']], 
                                   pull_convergence(fit = m1))


m2 <- stan4bart::stan4bart(y ~ . -schoolid -ID + bart(. -schoolid -ID) + (1|schoolid),
                           data = sim[['worldA']][["data"]],
                           treatment = z,
                           cores = 1,
                           chains = 10,
                           iter = 4000
)

sim[['worldA']][["satt"]] <- cbind(sim[['worldA']][["satt"]], pull_estimates(fit = m2)[['satt_est']])

sim[["worldA"]][["gsatt"]] <- cbind(sim[["worldA"]][["gsatt"]],
                                    merge(sim[['worldA']][["gsatt"]], pull_estimates(fit = m2)[['gsatt_est']]))

sim[["worldA"]][["icatt"]] <- cbind(sim[["worldA"]][["icatt"]], 
                                    merge(sim[['worldA']][["icatt"]], pull_estimates(fit = m2)[['icatt_est']]))

sim[["worldA"]][['rhat']] <- rbind(sim[["worldA"]][['rhat']], 
                                   pull_convergence(fit = m2))



m3 <- stan4bart::stan4bart(y ~ z + bart(. -schoolid -ID) + (1|schoolid),
                           data = sim[['worldA']][["data"]],
                           treatment = z,
                           cores = 1,
                           chains = 10,
                           iter = 4000
)

sim[['worldA']][["satt"]] <- cbind(sim[['worldA']][["satt"]], pull_estimates(fit = m3)[['satt_est']])

sim[["worldA"]][["gsatt"]] <- cbind(sim[["worldA"]][["gsatt"]],
                                    merge(sim[['worldA']][["gsatt"]], pull_estimates(fit = m3)[['gsatt_est']]))

sim[["worldA"]][["icatt"]] <- cbind(sim[["worldA"]][["icatt"]], 
                                    merge(sim[['worldA']][["icatt"]], pull_estimates(fit = m3)[['icatt_est']]))

sim[["worldA"]][['rhat']] <- rbind(sim[["worldA"]][['rhat']], 
                                   pull_convergence(fit = m3))



m4 <- stan4bart::stan4bart(y ~ bart(. -schoolid -ID) + (z|schoolid),
                           data = sim[['worldA']][["data"]],
                           treatment = z,
                           cores = 1,
                           chains = 10,
                           iter = 4000
)

sim[['worldA']][["satt"]] <- cbind(sim[['worldA']][["satt"]], pull_estimates(fit = m4)[['satt_est']])

sim[["worldA"]][["gsatt"]] <- cbind(sim[["worldA"]][["gsatt"]],
                                    merge(sim[['worldA']][["gsatt"]], pull_estimates(fit = m4)[['gsatt_est']]))

sim[["worldA"]][["icatt"]] <- cbind(sim[["worldA"]][["icatt"]], 
                                    merge(sim[['worldA']][["icatt"]], pull_estimates(fit = m4)[['icatt_est']]))

sim[["worldA"]][['rhat']] <- rbind(sim[["worldA"]][['rhat']], 
                                   pull_convergence(fit = m4))


m5 <- stan4bart::stan4bart(y ~ . -schoolid -ID + bart(. -schoolid -ID) + (z|schoolid),
                           data = sim[['worldA']][["data"]],
                           treatment = z,
                           cores = 1,
                           chains = 10,
                           iter = 4000
)

sim[['worldA']][["satt"]] <- cbind(sim[['worldA']][["satt"]], pull_estimates(fit = m5)[['satt_est']])

sim[["worldA"]][["gsatt"]] <- cbind(sim[["worldA"]][["gsatt"]],
                                    merge(sim[['worldA']][["gsatt"]], pull_estimates(fit = m5)[['gsatt_est']]))

sim[["worldA"]][["icatt"]] <- cbind(sim[["worldA"]][["icatt"]], 
                                    merge(sim[['worldA']][["icatt"]], pull_estimates(fit = m5)[['icatt_est']]))

sim[["worldA"]][['rhat']] <- rbind(sim[["worldA"]][['rhat']], 
                                   pull_convergence(fit = m5))



m6 <- stan4bart::stan4bart(y ~ z + bart(. -schoolid -ID) + (z|schoolid),
                           data = sim[['worldA']][["data"]],
                           treatment = z,
                           cores = 1,
                           chains = 10,
                           iter = 4000
)

sim[['worldA']][["satt"]] <- cbind(sim[['worldA']][["satt"]], pull_estimates(fit = m6)[['satt_est']])

sim[["worldA"]][["gsatt"]] <- cbind(sim[["worldA"]][["gsatt"]],
                                    merge(sim[['worldA']][["gsatt"]], pull_estimates(fit = m6)[['gsatt_est']]))

sim[["worldA"]][["icatt"]] <- cbind(sim[["worldA"]][["icatt"]], 
                                    merge(sim[['worldA']][["icatt"]], pull_estimates(fit = m6)[['icatt_est']]))

sim[["worldA"]][['rhat']] <- rbind(sim[["worldA"]][['rhat']], 
                                   pull_convergence(fit = m6))


m7 <- stan4bart::stan4bart(y ~ z + bart(. -schoolid -ID -z) + (z|schoolid),
                           data = sim[['worldA']][["data"]],
                           treatment = z,
                           cores = 1,
                           chains = 10,
                           iter = 4000
)

sim[['worldA']][["satt"]] <- cbind(sim[['worldA']][["satt"]], pull_estimates(fit = m7)[['satt_est']])

sim[["worldA"]][["gsatt"]] <- cbind(sim[["worldA"]][["gsatt"]],
                                    merge(sim[['worldA']][["gsatt"]], pull_estimates(fit = m7)[['gsatt_est']]))

sim[["worldA"]][["icatt"]] <- cbind(sim[["worldA"]][["icatt"]], 
                                    merge(sim[['worldA']][["icatt"]], pull_estimates(fit = m7)[['icatt_est']]))

sim[["worldA"]][['rhat']] <- rbind(sim[["worldA"]][['rhat']], 
                                   pull_convergence(fit = m7))

m8 <- stan4bart::stan4bart(y ~ z + bart(. -schoolid -ID -z) + (1|schoolid),
                           data = sim[['worldA']][["data"]],
                           treatment = z,
                           cores = 1,
                           chains = 10,
                           iter = 4000
)

sim[['worldA']][["satt"]] <- cbind(sim[['worldA']][["satt"]], pull_estimates(fit = m8)[['satt_est']])

sim[["worldA"]][["gsatt"]] <- cbind(sim[["worldA"]][["gsatt"]],
                                    merge(sim[['worldA']][["gsatt"]], pull_estimates(fit = m8)[['gsatt_est']]))

sim[["worldA"]][["icatt"]] <- cbind(sim[["worldA"]][["icatt"]], 
                                    merge(sim[['worldA']][["icatt"]], pull_estimates(fit = m8)[['icatt_est']]))

sim[["worldA"]][['rhat']] <- rbind(sim[["worldA"]][['rhat']], 
                                   pull_convergence(fit = m8))




write_rds(dat, glue::glue('results/results_iter{iter}.rds'), compress = 'gz')















