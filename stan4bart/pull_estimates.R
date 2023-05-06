pull_estimates <- function(fit){
  
  model <- deparse(substitute(fit)) # this will help us keep track of whats what
  trt <- fit$frame$z
  # icates
  mu.train <- stan4bart:::extract.stan4bartFit(fit, sample = 'train', type = 'ev')
  mu.test <- stan4bart:::extract.stan4bartFit(fit, sample = 'test', type = 'ev')
  
  icate <- (mu.train - mu.test) * (2 * trt - 1)
  icatt <- icate[trt == 1,]
  icatt_est <- apply(icatt, 1, mean)
  icatt_est <- cbind(icatt_est, ID = fit$frame$ID[trt == 1])
  names(icatt_est)[1] <- paste(names(icatt_est)[1], model, sep = '_')
  
  # satt 
  ppd.cf <- stan4bart:::extract.stan4bartFit(fit, type = "ppd", sample = "test")
  
  # Individual sample treatment effects
  samples.ite <- (fit$frame$y - ppd.cf) * (2 * trt - 1)
  
  samples.ite <- samples.ite[trt == 1, ]
  
  # Sample average treatment effect
  samples.satt <- apply(samples.ite, 2, mean)
  # pull estimates
  satt.est <- mean(samples.satt)
  satt.lci <- unname(quantile(samples.satt, .025))
  satt.uci <- unname(quantile(samples.satt, .975))
  
  # save satt
  satt_est <- cbind.data.frame(satt.est, satt.lci, satt.uci)
  names(satt_est) <- paste(names(satt_est), model, sep = '_')
  
  gsatt <- cbind.data.frame(schoolid = fit$frame[trt == 1, 'schoolid'], 
                            samples.ite) 
  
  names(gsatt) <- c('schoolid', paste0('sample', 1:ncol(samples.ite))) 
  
  gsatt_est <- gsatt %>% 
    pivot_longer(cols = contains('sample')) %>% 
    group_by(schoolid) %>% 
    summarise(gsatt.est = mean(value), 
              gsatt.lci = quantile(value, .025), 
              gsatt.uci = quantile(value, .975)) 
  
  names(gsatt_est)[2:4] <- paste(names(gsatt_est)[2:4], model, sep = '_')
  
  
  estimates <- list('satt_est' = satt_est, 
                    'gsatt_est' = gsatt_est, 
                    'icatt_est' = icatt_est)

  return(estimates)
  
}