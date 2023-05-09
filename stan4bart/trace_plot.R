trace <- function(fit, estimand){

  ppd.cf <- stan4bart:::extract.stan4bartFit(fit, 
                                             "ppd", 
                                             combine_chains = F, 
                                             sample = 'test')
  
  trt <- fit$frame$z
  
  samples.ite <- (fit$frame$y - ppd.cf) * (2 * trt - 1)
  
  samples.ite <- samples.ite[trt == 1,,]
  
  array(
    apply(samples.ite, c(2, 3), mean), 
    dim = c(fit$call$iter, fit$call$chains, 1)) %>% 
    as_tibble() %>% 
    mutate(iteration = 1:fit$call$iter/2) %>% 
    pivot_longer(1:fit$call$chains) %>% 
    ggplot(aes(iteration, value, col = name)) +
    geom_line() + 
    labs(title = paste0(deparse(formula(fit)))) + 
    theme_bw()
  

}
