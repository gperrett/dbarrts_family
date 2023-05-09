pull_convergence <- function(fit){
  model <- deparse(substitute(fit)) # this will help us keep track of whats what
  
  ppd.cf <- stan4bart:::extract.stan4bartFit(fit, 
                                             "ppd", 
                                             combine_chains = F, 
                                             sample = 'test')
  
  trt <- fit$frame$z
  
  samples.ite <- (fit$frame$y - ppd.cf) * (2 * trt - 1)
  samples.ite <- samples.ite[trt == 1,,]
  rhat <- rstan::monitor(
    array(
      apply(samples.ite, c(2, 3), mean), 
      dim = c(fit$call$iter/2, fit$call$chains, 1)),
    warmup = 0)
  
  rhat$model <- model
  
  return(rhat)
  
}

