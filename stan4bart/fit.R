library(stan4bart)
`%notin%` <- Negate(`%in%`)

setwd('~/Dropbox/dbarts_family/stan4bart/')
source('draw.R')
source('pull_estimates.R')
source('pull_convergence .R')


iteration <- Sys.getenv('SLURM_ARRAY_TASK_ID')

set.seed(iteration)

sim <- draw_dat()

worlds <- c('worldA', 'worldB', 'worldC')

results <- lapply(worlds, function(world){

    m1 <- stan4bart::stan4bart(y ~ bart(. -schoolid -ID) + (1|schoolid),
                             data = sim[[world]][["data"]],
                             treatment = z,
                             cores = 1,
                             chains = 10,
                             iter = 4000
  )
  
  
  sim[[world]][["satt"]] <<- cbind(sim[[world]][["satt"]], pull_estimates(fit = m1)[['satt_est']])
  
  sim[[world]][["gsatt"]] <<- merge(sim[[world]][["gsatt"]], pull_estimates(fit = m1)[['gsatt_est']])
  
  sim[[world]][["icatt"]] <<- merge(sim[[world]][["icatt"]], pull_estimates(fit = m1)[['icatt_est']])
  
  sim[[world]][['rhat']] <<- rbind(sim[[world]][['rhat']], 
                                     pull_convergence(fit = m1))
  
  rm(m2)
  gc()
  m2 <- stan4bart::stan4bart(y ~ . -schoolid -ID + bart(. -schoolid -ID) + (1|schoolid),
                             data = sim[[world]][["data"]],
                             treatment = z,
                             cores = 1,
                             chains = 10,
                             iter = 4000
  )
  
  sim[[world]][["satt"]] <<- cbind(sim[[world]][["satt"]], pull_estimates(fit = m2)[['satt_est']])
  
  sim[[world]][["gsatt"]] <<- merge(sim[[world]][["gsatt"]], pull_estimates(fit = m1)[['gsatt_est']])
  
  sim[[world]][["icatt"]] <<- merge(sim[[world]][["icatt"]], pull_estimates(fit = m1)[['icatt_est']])
  
  
  sim[[world]][['rhat']] <<- rbind(sim[[world]][['rhat']], 
                                     pull_convergence(fit = m2))
  
  
  rm(m2)
  gc()
  m3 <- stan4bart::stan4bart(y ~ z + bart(. -schoolid -ID) + (1|schoolid),
                             data = sim[[world]][["data"]],
                             treatment = z,
                             cores = 1,
                             chains = 10,
                             iter = 4000
  )
  
  sim[[world]][["satt"]] <<- cbind(sim[[world]][["satt"]], pull_estimates(fit = m3)[['satt_est']])
  
  sim[[world]][["gsatt"]] <<- merge(sim[[world]][["gsatt"]], pull_estimates(fit = m1)[['gsatt_est']])
  
  sim[[world]][["icatt"]] <<- merge(sim[[world]][["icatt"]], pull_estimates(fit = m1)[['icatt_est']])
  
  sim[[world]][['rhat']] <<- rbind(sim[[world]][['rhat']], 
                                     pull_convergence(fit = m3))
  
  
  rm(m3)
  gc()
  m4 <- stan4bart::stan4bart(y ~ bart(. -schoolid -ID) + (z|schoolid),
                             data = sim[[world]][["data"]],
                             treatment = z,
                             cores = 1,
                             chains = 10,
                             iter = 4000
  )
  
  sim[[world]][["satt"]] <<- cbind(sim[[world]][["satt"]], pull_estimates(fit = m4)[['satt_est']])
  
  sim[[world]][["gsatt"]] <<- merge(sim[[world]][["gsatt"]], pull_estimates(fit = m1)[['gsatt_est']])
  
  sim[[world]][["icatt"]] <<- merge(sim[[world]][["icatt"]], pull_estimates(fit = m1)[['icatt_est']])
  
  sim[[world]][['rhat']] <<- rbind(sim[[world]][['rhat']], 
                                     pull_convergence(fit = m4))
  
  rm(m4)
  gc()
  m5 <- stan4bart::stan4bart(y ~ . -schoolid -ID + bart(. -schoolid -ID) + (z|schoolid),
                             data = sim[[world]][["data"]],
                             treatment = z,
                             cores = 1,
                             chains = 10,
                             iter = 4000
  )
  
  sim[[world]][["satt"]] <<- cbind(sim[[world]][["satt"]], pull_estimates(fit = m5)[['satt_est']])
  
  sim[[world]][["gsatt"]] <<- merge(sim[[world]][["gsatt"]], pull_estimates(fit = m1)[['gsatt_est']])
  
  sim[[world]][["icatt"]] <<- merge(sim[[world]][["icatt"]], pull_estimates(fit = m1)[['icatt_est']])
  
  sim[[world]][['rhat']] <<- rbind(sim[[world]][['rhat']], 
                                     pull_convergence(fit = m5))
  
  
  rm(m5)
  gc()
  m6 <- stan4bart::stan4bart(y ~ z + bart(. -schoolid -ID) + (z|schoolid),
                             data = sim[[world]][["data"]],
                             treatment = z,
                             cores = 1,
                             chains = 10,
                             iter = 4000
  )
  
  sim[[world]][["satt"]] <<- cbind(sim[[world]][["satt"]], pull_estimates(fit = m6)[['satt_est']])
  
  sim[[world]][["gsatt"]] <<- merge(sim[[world]][["gsatt"]], pull_estimates(fit = m1)[['gsatt_est']])
  
  sim[[world]][["icatt"]] <<- merge(sim[[world]][["icatt"]], pull_estimates(fit = m1)[['icatt_est']])
  
  sim[[world]][['rhat']] <<- rbind(sim[[world]][['rhat']], 
                                     pull_convergence(fit = m6))
  
  rm(m6)
  gc()
  m7 <- stan4bart::stan4bart(y ~ z + bart(. -schoolid -ID -z) + (z|schoolid),
                             data = sim[[world]][["data"]],
                             treatment = z,
                             cores = 1,
                             chains = 10,
                             iter = 4000
  )
  
  sim[[world]][["satt"]] <<- cbind(sim[[world]][["satt"]], pull_estimates(fit = m7)[['satt_est']])
  
  sim[[world]][["gsatt"]] <<- merge(sim[[world]][["gsatt"]], pull_estimates(fit = m1)[['gsatt_est']])
  
  sim[[world]][["icatt"]] <<- merge(sim[[world]][["icatt"]], pull_estimates(fit = m1)[['icatt_est']])
  
  sim[[world]][['rhat']] <<- rbind(sim[[world]][['rhat']], 
                                     pull_convergence(fit = m7))
  
  rm(m7)
  gc()
  
  m8 <- stan4bart::stan4bart(y ~ z + bart(. -schoolid -ID -z) + (1|schoolid),
                             data = sim[[world]][["data"]],
                             treatment = z,
                             cores = 1,
                             chains = 10,
                             iter = 4000
  )
  
  sim[[world]][["satt"]] <<- cbind(sim[[world]][["satt"]], pull_estimates(fit = m8)[['satt_est']])
  
  sim[[world]][["gsatt"]] <<- merge(sim[[world]][["gsatt"]], pull_estimates(fit = m1)[['gsatt_est']])
  
  sim[[world]][["icatt"]] <<- merge(sim[[world]][["icatt"]], pull_estimates(fit = m1)[['icatt_est']])
  
  sim[[world]][['rhat']] <<- rbind(sim[[world]][['rhat']], 
                                     pull_convergence(fit = m8))


})


write_rds(results, glue::glue('results/results_iter{iteration}.rds'), compress = 'gz')















