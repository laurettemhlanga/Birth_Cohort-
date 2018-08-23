packages <- c("deSolve", "ggplot2", "plyr", "glm2", "splines", "reshape2", "inctools")
lapply(packages, require, character.only = TRUE)


Birthrate <- function(t, birthrate = 0){
  birth = birthrate 
  return(birth)
}



Incidence_var <- function(t, const_inc=0.02){ 
  inc <- const_inc
  return (inc)
}



Backgrnd_Mortality_var <- function(t, Mort_peak = 0.01){
  Mort =  Mort_peak 
  return (Mort)
}



Exit_Recency_rate <-  function(t, exit_rate = 2){
  exit <- exit_rate
  return(exit)
} 




Excess_Mortality_var <- function(t, EMort_peak = 0.1
                                 ){
  EMort <- EMort_peak
  return (EMort)
}







base_model_recency <- function(time, y.vec, parms){
  #simulates an HIV epidemic taking the recency and non receny form into consideration.
  with(c(as.list(y.vec)),{
    dSdt <- Birthrate(time)-(Incidence_var(time) + Backgrnd_Mortality_var(time)) * Susceptible
    dRdt <- Incidence_var(time) * Susceptible - (Backgrnd_Mortality_var(time) + Excess_Mortality_var(time) +Exit_Recency_rate(time)) * Recent
    dNrdt <- Exit_Recency_rate(time) * Recent -(Backgrnd_Mortality_var(time) +Excess_Mortality_var(time))* Non_Recent
   
    list(c(dSdt, dRdt, dNrdt))
  })
}



base_model <- function(time, y.vec, parms){
  #simulates an HIV epidemic taking the recency and non receny form into consideration.
  with(c(as.list(y.vec)),{
    dSdt <- Birthrate(time)-(Incidence_var(time) + Backgrnd_Mortality_var(time)) * Susceptible
    dIdt <- Incidence_var(time) * Susceptible - (Backgrnd_Mortality_var(time) + Excess_Mortality_var(time)) * Infected
    list(c(dSdt, dIdt))
  })
}

#ran the simulations till equilibrium and the estimated incidence at age 20 and age 40 to compare the OISI 

Sim_Birth_Cohort <- function( Mintime = 15, Maxtime = 50, DT = 1,
                              #Pop_initial = c(Susceptible = 1,  Recent = 0, Non_Recent = 0, Infected = 0),
                              #Pop_initial = c(Susceptible = 1, Infected = 0),
                              Pop_initial = c(Susceptible = 1, Recent = 0, Non_Recent = 0),
                              model_function = NULL,
                              Recency = TRUE, Spline = TRUE
){
  
  #simulates a population with and SI model / SRN (susceptible, Receny and Non Receny) model
  #Outputs lists of spline functions or the population dataframe 
  #equilibrium  Recent = 0.00474, Non_Recent = 0.0086169
  
  if(Recency == TRUE){
    Population <- data.frame(lsoda(
      y = Pop_initial,              
      times = seq(Mintime, Maxtime, by = DT),             
      func = if(is.null(model_function)){base_model_recency},
      parms = NULL))
    Population$Prevalence = (Population$Recent + Population$Non_Recent)/(Population$Susceptible + Population$Recent + Population$Non_Recent)
    Population$Prev_Recency = Population$Recent /(Population$Recent + Population$Non_Recent)
    #Population = Population[,-5]
    Splines = as.list(c( Prev_Recency_Function = splinefun(x = Population$time, y = Population$Prev_Recency ),
                         Susceptible_Function = splinefun(x = Population$time, y = Population$Susceptible),
                         Prevalence_Function = splinefun(x = Population$time, y = Population$Prevalence),
                         Non_Recency_Function = splinefun(x = Population$time, y = Population$Non_Recent)))
  }else{
    Population <- data.frame(lsoda(
      y = Pop_initial,              
      times = seq(Mintime, Maxtime, by = DT),             
      func = if(is.null(model_function)){base_model},
      parms=NULL))
    Population$Prevalence = (Population$Infected)/(Population$Susceptible + Population$Infected)
    #Population = Population[,-3:-4]
    Splines = as.list(c(Susceptible_Function = splinefun(x = Population$time, y = Population$Susceptible),
                        Infected_Function = splinefun(x = Population$time, y = Population$Infected),
                        Prevalence_Function = splinefun(x = Population$time, y = Population$Prevalence) ))
  }
  if (Spline == TRUE){
    return(Splines)
  }else{
    return(Population)
  }
}

# ADD extraction of  the spline functions created.




Mincidence_PE <- function(prev, slope, excess_mortality){
  # calculates incidence only Mahiane estimator.
  return( 1 / (1 - prev) * slope + excess_mortality * prev)
}


Mincidence_SE <- function( Incidence, prev, 
                           slope, excess_mortality, 
                           prev_error, slope_error, excess_mortality_error = 0){
  #delta method to yeild the standard error of the Mahiane estimator
  return(sqrt(
    ((1-prev)^(-2) * slope + excess_mortality)^2 * (prev_error)^2  
    + ((1-prev)^(-2) * slope_error^2)
    + (prev*excess_mortality_error)^2
  )
  )
}






Midpoint_inc_from_prevalence <- function(N_iterations = 1, T1, T2, excess_mortality_estimate,  
                                         sample_size_1 = 1000, sample_size_2 = 1000, pop_prevalence,
                                         return_boots = FALSE, alpha = 0.05, DE1 = 1, DE2 = 1,
                                         sigma_prev_1 = 0, sigma_prev_2 = 0, Delta_T = 0
){
  # browser()
  
  
  if (is.vector(pop_prevalence) == TRUE){
    
    
    midpoint_prevalence = (pop_prevalence[1] + pop_prevalence[2]) / 2
    prevalence_derivative = (pop_prevalence[2] - pop_prevalence[1] ) / (T2 - T1)
    
    midpoint_prevalence_se = sqrt((0.5 ^ 2 *  ((sigma_prev_1) ^ 2) / DE1) + (0.5 ^ 2 *  ((sigma_prev_2) ^ 2) / DE2))
    prevalence_derivative_se = sqrt(( (1 / (T2 - T1)) ^ 2 * ((sigma_prev_1) ^ 2) / DE1) + (( -1 / (T2 - T1)) ^ 2* ((sigma_prev_2) ^ 2) / DE2))
    
    
    incidence = ((1 / (1 - midpoint_prevalence)) * prevalence_derivative) + midpoint_prevalence * excess_mortality_estimate
    incidence_se = Mincidence_SE(prev = midpoint_prevalence, slope = prevalence_derivative, prev_error = midpoint_prevalence_se, 
                                 slope_error = prevalence_derivative_se, excess_mortality_estimate)
    
    results = data.frame(Incidence = incidence, Incidence_se = incidence_se)
    
  }else{
    t=(T1+T2)/2
    mid_time <- data.frame(times = c(t))
    
    timepoints <- c(T1,T2)
    
    
    HIV_Pos_Counts_Boot_1 <- rbinom(n = N_iterations, size = sample_size_1, prob = pop_prevalence(T1))
    HIV_Pos_Counts_Boot_2 <- rbinom(n = N_iterations, size = sample_size_2, prob = pop_prevalence(T2))
    
    HIV_Pos_Counts_Exp_1 <- round(pop_prevalence(timepoints[1]) * sample_size_1)
    HIV_Pos_Counts_Exp_2 <-  round(pop_prevalence(timepoints[2]) * sample_size_2)
    
    
    #
    #
    # process Counts_Exp
    
    if(N_iterations == 1){
      sim_surv_time_prev <- data.frame(times = timepoints, 
                                       sample_sizes = c(sample_size_1 ,sample_size_2), 
                                       HIV_Pos_Counts = c(HIV_Pos_Counts_Exp_1, HIV_Pos_Counts_Exp_2))
      
      fit_glm <- glm2(formula = cbind(HIV_Pos_Counts,sample_sizes-HIV_Pos_Counts) ~ 1 + I(times), 
                      data = sim_surv_time_prev, 
                      family = binomial(link = "identity"))
      
      Slope_Exp_Naive <- (pop_prevalence(timepoints[2]) - pop_prevalence(timepoints[1]))/(T2-T1)
      Slope_Exp_PE <- summary(fit_glm)$coefficients[2,1] 
      Slope_Exp_SE <- summary(fit_glm)$coefficients[2,2]
      
      predictions <- predict(fit_glm, newdata = mid_time, se.fit = TRUE)
      Pred_Midpoint_Prev_PE <- predictions$fit
      Pred_Midpoint_Prev_SE <- predictions$se.fit
      
      
      Incidence_PE <- Mincidence_PE(prev = Pred_Midpoint_Prev_PE,slope = Slope_Exp_PE,excess_mortality = excess_mortality_estimate)
      
      Incidence_SE <- Mincidence_SE( prev = Pred_Midpoint_Prev_PE, 
                                     slope = Slope_Exp_PE, excess_mortality= excess_mortality_estimate, 
                                     prev_error = Pred_Midpoint_Prev_SE, slope_error = Slope_Exp_SE)
      
      CI_Lower <- Incidence_PE - 1.96 * Incidence_SE
      CI_Upper <- Incidence_PE + 1.96 * Incidence_SE
      
      results <- data.frame(Method='DM', 
                            Slope_PE = Slope_Exp_PE, Slope_SE = Slope_Exp_SE, 
                            Mid_Prev_PE = Pred_Midpoint_Prev_PE, Mid_Prev_SE = Pred_Midpoint_Prev_SE,
                            Incidence_PE = Incidence_PE, Incidence_SE = Incidence_SE,
                            CI_Upper = CI_Upper, CI_Lower = CI_Lower)
      
      
    }else{
      
      Slope_PE <- as.vector(rep(NA, N_iterations))
      Slope_SE <- as.vector(rep(NA, N_iterations))
      
      Pred_Prevalences_PE <- as.vector(rep(NA, N_iterations))
      Pred_Prevalences_SE <- as.vector(rep(NA, N_iterations))
      
      Incidence <- as.vector(rep(NA, N_iterations))
      Incidence_SE<- as.vector(rep(NA, N_iterations))
      CI_Lower <-  as.vector(rep(NA, N_iterations))
      CI_Upper <-  as.vector(rep(NA, N_iterations))
      
      samplesizes <- c(sample_size_1,sample_size_2)
      
      for(i in 1:N_iterations) {
        
        sim_surv_time_prev <- data.frame(timepoints, samplesizes, c(HIV_Pos_Counts_Boot_1[i],HIV_Pos_Counts_Boot_2[i]))
        
        fit_glm <- glm2(cbind(c(HIV_Pos_Counts_Boot_1[i],HIV_Pos_Counts_Boot_2[i]), 
                              c(sample_size_1-HIV_Pos_Counts_Boot_1[i], sample_size_2-HIV_Pos_Counts_Boot_2[i])) ~ 1 + 
                          I(timepoints), data = sim_surv_time_prev , family = binomial(link = "identity"))
        
        Slope_PE[i] <- summary(fit_glm)$coefficients[2,1] 
        Slope_SE[i] <- summary(fit_glm)$coefficients[2,2]
        
        predictions <- predict(fit_glm, newdata = data.frame(timepoints = c(t)), se.fit = TRUE)
        Pred_Prevalences_PE[i] <- predictions$fit
        Pred_Prevalences_SE[i] <- predictions$se.fit
        
        Incidence[i] <- Mincidence_PE(prev = Pred_Prevalences_PE[i], slope = Slope_PE[i],
                                      excess_mortality = excess_mortality_estimate)
        
        Incidence_SE[i] <- Mincidence_SE(prev = Pred_Prevalences_PE[i], 
                                         slope = Slope_PE[i], excess_mortality= excess_mortality_estimate, 
                                         prev_error = Pred_Prevalences_SE[i], slope_error = Slope_SE[i])
        
        CI_Lower[i] <- Incidence[i] - 1.96 * Incidence_SE[i]
        CI_Upper[i] <- Incidence[i] + 1.96 * Incidence_SE[i]
      }  
      
      Iterate_results <- data.frame(Method='Iterate', 
                                    Slope_PE = Slope_PE, Slope_SE = Slope_SE, 
                                    Mid_Prev_PE = Pred_Prevalences_PE, 
                                    Mid_Prev_SE = Pred_Prevalences_SE, 
                                    Incidence_PE = Incidence, Incidence_SE = Incidence_SE,
                                    CI_Upper = CI_Upper, CI_Lower = CI_Lower)
      
      
      
      
      if (return_boots == TRUE){
        
        results <- Iterate_results
        #results <- (Iterate_results)
      }
      
      # SUmmarise bootstrap iterations
      Boot_Inc_mean <- mean(Iterate_results$Incidence_PE)
      #Boot_Prev_mean <- mean(Iterate_results$Mid_Prev_PE)
      Boot_SE <-   sd(Iterate_results$Incidence_PE)
      Boot_CI <-   quantile(x = Iterate_results$Incidence_PE, probs = c(alpha/2, 1-alpha/2),type=2, names=FALSE) 
      
      Boot_summary  <- data.frame(Method='Boot', 
                                  Slope_PE = NA, Slope_SE = NA, 
                                  Mid_Prev_PE = NA, Mid_Prev_SE = NA,
                                  Incidence_PE = Boot_Inc_mean, Incidence_SE = Boot_SE,
                                  CI_Upper = Boot_CI[2], CI_Lower = Boot_CI[1])
      
      results <- rbind(results,Boot_summary)
    }}
  return(results)
  
}







