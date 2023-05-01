rm(list = ls())


####Read in Libraries####
library(tidyverse)
library(MASS)


####Read in Data####
dat_raw <- readr::read_csv('data/classroom.csv') %>% 
  na.omit() %>% 
  as.data.frame() 


####Start Function####
satt_function <- function(group_vars = c("schoolid"), random_slope = TRUE, sd1 = 1, sd2 = 2, cov1 = 0.3, cov2 = 0.2){
  
  ####Set Up for Creation of Response Surfaces####
  # filter out groups and save for later
  grp_names <- group_vars
  total_group_cols <- c("schoolid", "classid", "childid")
  
  groups <- 
    dat_raw %>% 
    dplyr::select(all_of(grp_names))
  
  dat <- 
    dat_raw %>% 
    dplyr::select(-all_of(total_group_cols))
  
  # standardize continuous vars
  continuous <- 
    dat %>% 
    dplyr::select(c('mathkind', 'mathgain', 'ses', 'mathknow', 'housepov', 'mathprep', 'yearstea')) %>% 
    apply(2, scale)
  
  # get categorical vars 
  cat <- dat[, c('minority', 'sex')] 
  
  # combine to one matrix for some diagnostics
  X <- cbind(cat, continuous)
  
  # things become a lot easier if everything is positvily correlated 
  # so we will recode things so they are in the same direction 
  cor.mat <- cor(X, method = 'spearman')
  cor.mat 
  
  X[, which(cor.mat[, 'housepov'] < 0)] <- X[, which(cor.mat[, 'housepov'] < 0)]*-1
  cor(X, method = 'spearman') # double check recode worked...
  
  
  # Identidy a tight area of the covariate space for treated units
  treated_groups <- cbind(
    ifelse(continuous[, 'ses'] < -1 , 1, 0),
    ifelse(continuous[, 'mathkind'] < -1,1, 0),
    ifelse(continuous[, 'mathgain'] < -1,1, 0),
    ifelse(continuous[, 'mathknow'] < -1 , 1, 0),
    ifelse(continuous[, 'housepov'] < -1 , 1, 0),
    ifelse(continuous[, 'mathprep'] < -1 , 1, 0),
    ifelse(continuous[, 'yearstea'] < -1 , 1, 0),
    ifelse(cat[, 'minority'] == 1 , 1, 0),
    ifelse(cat[, 'minority'] == 1 & cat[, 'sex'] == 1 , 1, 0)
  )
  
  p.score <- dplyr::case_when(
    rowSums(treated_groups) == 0 ~ 0, 
    rowSums(treated_groups) > 0 & rowSums(treated_groups) < 3 ~ .55,
    rowSums(treated_groups) ==  3  ~ .65, 
    T ~ .45
  )
  
  # now we'll refine a bit further: this is like an interaction terms
  p.score <- ifelse(continuous[, 'ses'] < .5, p.score, 0) 
  p.score <- ifelse(continuous[, 'housepov'] < .5, p.score, 0) 
  p.score <- ifelse(continuous[, 'mathknow'] < .5, p.score, 0) 
  p.score <- ifelse(continuous[, 'mathgain'] < .5, p.score, 0) 
  p.score <- ifelse(continuous[, 'mathkind'] < 1, p.score, 0) 
  p.score <- ifelse(continuous[, 'mathkind'] < 1.5, p.score, 0) 
  # hist(p.score)
  
  z <- rbinom(nrow(dat), 1, p.score)
  
  ##Create Group Error Structure (Just School)
  Sigma_school <- matrix(c(sd1, cov1, cov1, sd1), 2, 2)
  Sigma_class <- matrix(c(sd2, cov2, cov2, sd2), 2, 2) 
  school_error <- 
    mvrnorm(n = length(unique(groups$schoolid)), mu = c(0, 0), Sigma = Sigma_school) %>% 
    as.data.frame() %>% 
    rename("school_e" = "V1",
           "school_z" = "V2") %>%
    mutate(schoolid = unique(groups$schoolid)) %>% 
    relocate(schoolid, .before = "school_e") 
  
  if("classid" %in% group_vars){
    school_error_new <-
      groups %>%
      unique() %>%
      left_join(school_error, by = "schoolid") %>%
      group_by(schoolid) %>%
      unique()
    
    class_error_list <- list()
    for(i in seq_len(nrow(school_error))){
      school <- school_error$schoolid[i]
      class_groups <-
        groups %>%
        filter(schoolid == school) %>%
        unique()
      
      if(nrow(class_groups) == 1){
        df <-
          class_groups %>%
          bind_cols(mvrnorm(nrow(class_groups), mu = c(0, 0), Sigma = Sigma_class) %>%
                      t() %>%
                      as.data.frame() %>%
                      rename("class_e" = "V1",
                             "class_z" = "V2"))
      }else{
        df <-
          class_groups %>%
          bind_cols(mvrnorm(nrow(class_groups), mu = c(0, 0), Sigma = Sigma_class) %>%
                      as.data.frame() %>%
                      rename("class_e" = "V1",
                             "class_z" = "V2"))
      }
      class_error_list[[i]] <- df
    }
    class_error_df <-
      class_error_list %>%
      bind_rows()
    
    group_error <-
      class_error_df %>%
      left_join(school_error, by = "schoolid") %>%
      ungroup() %>%
      arrange(schoolid, classid) %>%
      mutate(zeta = class_e + school_e,
             rand_slope = class_z + school_z)
    
    group_error_final <- group_error
  }else{
    group_error_final <- 
      school_error %>% 
      mutate(zeta = school_e,
             rand_slope = school_z)
  }
  
  
  ####Create Response Surfaces####
  ##Linear Response Surface
  Xmat <- cbind(rep(1, nrow(X)), as.matrix(X))
  
  # generate betas 
  probA_raw <- rpois(n = 7, lambda = 2)
  probA <- probA_raw/sum(probA_raw)
  betaA <- sample(c(-1:5), ncol(Xmat), prob = probA, replace = TRUE)
  
  #z <- rbinom(nrow(X), 1, .5)
  # z_probs <- 
  #   X %>% 
  #   as.data.frame() %>% 
  #   ungroup() %>% 
  #   rowSums() %>% 
  #   scales::rescale(to = c(0, 1))
  # z <- rbinom(nrow(X), 1, z_probs)
  ya1.hat <- Xmat %*% betaA
  ya0.hat <- Xmat %*% betaA
  
  # add error
  group_error_full <- 
    groups %>% 
    dplyr::select(all_of(group_vars)) %>% 
    left_join(group_error_final, by = group_vars)
  
  if(random_slope == TRUE){
    ya1.hat <- ya1.hat + group_error_full$zeta + group_error_full$rand_slope
    ya0.hat <- ya0.hat + group_error_full$zeta
  }else{
    ya1.hat <- ya1.hat + group_error_full$zeta 
    ya0.hat <- ya0.hat + group_error_full$zeta
  }
  
  # add treatment effect 
  treatment_effect <- sd(ya0.hat)*tau
  ya1.hat <- ya1.hat + treatment_effect
  # ya1.hat <- ya1.hat + sd(ya1.hat) * 0.5
  
  # add indv. error 
  ya1 <- ya1.hat + rnorm(length(ya1.hat), 0, 1)
  ya0 <- ya0.hat + rnorm(length(ya0.hat), 0, 1)
  
  ya <- ifelse(z == 1, ya1, ya0)
  datA <- 
    data.frame(ya, ya1, ya0, z, X, groups) %>% 
    rename(y = ya, 
           y1 = ya1, 
           y0 = ya0) %>% 
    mutate(across(c(group_vars), ~as.factor(.)))
  
  ##Non-Linear Response Surface
  Xquad <- 
    X %>% 
    cbind(X^2 %>% 
            as.data.frame() %>% 
            dplyr::select(-c(minority, sex)) %>% 
            rename_with(~paste0(.x, "_quad"), everything()))
  Xquad_mat <- cbind(rep(1, nrow(Xquad)), as.matrix(Xquad))
  
  # generate betas 
  probB_raw <- rpois(n = 9, lambda = 2)
  probB <- probB_raw/sum(probB_raw)
  betaB <- sample(c(-2:6), ncol(Xquad_mat), prob = probB, replace = TRUE)
  
  
  #z <- rbinom(nrow(X), 1, .5)
  # z_probs <- 
  #   X %>% 
  #   as.data.frame() %>% 
  #   ungroup() %>% 
  #   rowSums() %>% 
  #   scales::rescale(to = c(0, 1))
  # z <- rbinom(nrow(X), 1, z_probs)
  yb1.hat <- Xquad_mat %*% betaB
  yb0.hat <- Xquad_mat %*% betaB
  
  # add error
  if(random_slope == TRUE){
    yb1.hat <- yb1.hat + group_error_full$zeta + group_error_full$rand_slope
    yb0.hat <- yb0.hat + group_error_full$zeta
  }else{
    yb1.hat <- yb1.hat + group_error_full$zeta 
    yb0.hat <- yb0.hat + group_error_full$zeta
  }
  
  
  # add treatment effect 
  treatment_effect <- sd(yb0.hat)*tau
  yb1.hat <- yb1.hat + treatment_effect
  # ya1.hat <- ya1.hat + sd(ya1.hat) * 0.5
  
  # add indv. error 
  yb1 <- yb1.hat + rnorm(length(yb1.hat))
  yb0 <- yb0.hat + rnorm(length(yb0.hat))
  
  yb <- ifelse(z == 1, yb1, yb0)
  datB <- 
    data.frame(yb, yb1, yb0, z, X, groups) %>% 
    rename(y = yb, 
           y1 = yb1, 
           y0 = yb0) %>% 
    mutate(across(c(group_vars), ~as.factor(.)))
  
  ##Non-Linear Response Surface w/ Interactions
  ytmp = rnorm(nrow(X))
  mod_bal <- glm(ytmp ~ (minority + sex + mathkind + mathgain + ses + mathknow + housepov + mathprep + yearstea)^2 + 
                   I(mathkind^2) + I(mathgain^2) + I(ses^2) + I(mathknow^2) + I(housepov^2) + I(mathprep^2) + I(yearstea^2), data = X, x = TRUE)
  coefs <- mod_bal$coef[-1]
  XX <- mod_bal$x[,-1]
  XX <- XX[,!is.na(coefs)]
  XXmat <- cbind(rep(1, nrow(X)), XX)
  
  # generate betas for main effects 
  probC1_raw <- rpois(n = 9, lambda = 2)
  probC1 <- probC1_raw/sum(probC1_raw)
  betaC1 <- sample(c(-2:6), ncol(Xquad_mat), prob = probC1, replace = TRUE)
  
  # generate betas for interaction effects
  betaC2 <- sample(c(-0.5, 0, 0.5, 1), ncol(XXmat) - ncol(Xquad_mat), prob = c(0.1, 0.7, 0.1, 0.1), replace = TRUE)
  betaC <- c(betaC1, betaC2)
  
  yc1.hat <- XXmat %*% betaC
  yc0.hat <- XXmat %*% betaC
  
  # add error
  if(random_slope == TRUE){
    yc1.hat <- yc1.hat + group_error_full$zeta + group_error_full$rand_slope
    yc0.hat <- yc0.hat + group_error_full$zeta
  }else{
    yc1.hat <- yc1.hat + group_error_full$zeta
    yc0.hat <- yc0.hat + group_error_full$zeta
  }
  
  
  # add treatment effect 
  treatment_effect <- sd(yc0.hat)*tau
  yc1.hat <- yc1.hat + treatment_effect
  # ya1.hat <- ya1.hat + sd(ya1.hat) * 0.5
  
  # add indv. error 
  yc1 <- yc1.hat + rnorm(length(yc1.hat))
  yc0 <- yc0.hat + rnorm(length(yc0.hat))
  
  yc <- ifelse(z == 1, yc1, yc0)
  datC <- 
    data.frame(yc, yc1, yc0, z, X, groups) %>% 
    rename(y = yc, 
           y1 = yc1, 
           y0 = yc0) %>% 
    mutate(across(c(group_vars), ~as.factor(.)))
  
  
  ####Calculate estimands####
  ##SATT
  datA_satt <- mean(datA$y1[z == 1] - datA$y0[z == 1])
  datB_satt <- mean(datB$y1[z == 1] - datB$y0[z == 1])
  datC_satt <- mean(datC$y1[z == 1] - datC$y0[z == 1])  
  
  ## ICATT
  datA_icatt <- datA$y1[z == 1] - datA$y0[z == 1]
  datB_icatt <- datB$y1[z == 1] - datB$y0[z == 1]
  datC_icatt <- datC$y1[z == 1] - datC$y0[z == 1]
  
  
  ## GSATT
  datA_gsatt <- with(datA, tapply(y1[z == 1], schoolid[z == 1], mean) - tapply(y0[z == 1], schoolid[z == 1], mean))
  datA_gsatt <- data.frame(gsatt = datA_gsatt, schoolid = names(datA_gsatt))
  datB_gsatt <- with(datB, tapply(y1[z == 1], schoolid[z == 1], mean) - tapply(y0[z == 1], schoolid[z == 1], mean))
  datB_gsatt <- data.frame(gsatt = datB_gsatt, schoolid = names(datB_gsatt))
  datC_gsatt <- with(datC, tapply(y1[z == 1], schoolid[z == 1], mean) - tapply(y0[z == 1], schoolid[z == 1], mean))
  datC_gsatt <- data.frame(gsatt = datC_gsatt, schoolid = names(datC_gsatt))
  
  
  satt_results <- data.frame(response_group = c("A", "B", "C"),
                             SATT = c(datA_satt, datB_satt, datC_satt))
  
  icatt_results <- data.frame(response_group = c(rep("A", length(datA_icatt)), 
                                                 rep("B", length(datB_icatt)),
                                                 rep("C", length(datC_icatt))),
                              ICATT = c(datA_icatt, datB_icatt, datC_icatt))
  
  gsatt_results <- data.frame(response_group = c(rep("A", nrow(datA_gsatt)), 
                                                 rep("B", nrow(datB_gsatt)),
                                                 rep("C", nrow(datC_gsatt))),
                              GSATT = rbind(datA_gsatt, datB_gsatt, datC_gsatt))
  
  
  datA <- datA %>% dplyr::select(-y1, -y0)
  datB <- datB %>% dplyr::select(-y1, -y0)
  datC <- datC %>% dplyr::select(-y1, -y0)
  
  out <- list("True SATT" = satt_results,
                          "True ICATT" = icatt_results, 
                          "True GSATT" = gsatt_results, 
                          "Linear Response Surface" = datA,
                          "Non-Linear Response Surface" = datB,
                          "Non-Linear Response Surface w/ Interactions" = datC)
  
  return(out)
}


