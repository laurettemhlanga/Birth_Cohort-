#Analysis of the dementia study in German

incidence_dem <- function(age){
  incidence <-  exp(-12.8 + 0.11 * age)
  return(incidence)
}

incidence_dem(57.5)

mortality_dem_0 <- function(age, tim = 0){
  mortality  <-  exp(-9.0 + (0.085 * age) - (tim * log(1.01)))
  return(mortality)
} 


mortality_dem_1 <- function(age, tim = 0, r = 2.63){
  mortality  <-  mortality_dem_0(age, tim) * r
  return(mortality)
} 


mortality_dem_excess <- function(age, tim = 0, r = 2.63){
  mortality  <-  mortality_dem_1(age, tim, r) - mortality_dem_0(age, tim) 
  return(mortality)
} 



base_model_dem <- function(time, y.vec, parms){
  #simulates an HIV epidemic taking the recency and non receny form into consideration.
  with(c(as.list(y.vec)),{
    dSdt <- -(incidence_dem(time) + mortality_dem_0 (time)) * Susceptible
    dIdt <- incidence_dem(time) * Susceptible - (mortality_dem_0( time) + mortality_dem_excess(time)) * Infected
    list(c(dSdt, dIdt))
  })
}





#Calling the functions from revise_birth_cohort for the simulation for a population with   
population <- Sim_Birth_Cohort( Mintime = 50, Maxtime = 65, DT = 1, Pop_initial = c(Susceptible = 1, 
                                                                      Infected = 0), 
                               model_function = NULL, model = base_model_dem,
                               Recency = FALSE, Spline = TRUE)


#Extracting the spline functions for each attribute ofinterest
#when Recency false 
prop_Prevalence <- population$Prevalence_Function
Prop_NOn_Recency <- population$Prev_Non_Receny_Function



#calculating incidence 
Inci = Midpoint_inc_from_prevalence(N_iterations=1, T1 = 55, T2 = 60, 
                                    sample_size_1 = 10000, sample_size_2  = 10000,
                                    pop_prevalence = prop_Prevalence, 
                                    excess_mortality_estimate = mortality_dem_excess(57.5), 
                                    return_boots = TRUE)





RMSE_STDE_BIAS_Cal_1 <- function( N_iterations, min_dT, Max_dT,
                                  TT, sample_size_1 = 10000, Incidence, sample_size_2 = 10000, 
                                  pop_prevalence = prop_Prevalence, excess_mortality = mortality_dem_excess(57.5),
                                  return_boots = FALSE,
                                  Delta_t = seq(2, 12,2)){
  survey_intervals = seq(min_dT, Max_dT, 1) 
  N_intervals <- length(survey_intervals)
  std_errors <- as.vector(rep(NA, N_intervals))
  bias <- as.vector(rep(NA, N_intervals))
  RMSE <- as.vector(rep(NA, N_intervals))
  meansInc = as.vector(rep(NA, N_intervals))
  
  counter <- 1
  for (delta in survey_intervals){
    T1 <- TT - delta
    T2 = TT + delta 
    rel<- Midpoint_inc_from_prevalence(N_iterations = N_iterations, T1 = T1, T2 = T2, 
                                       sample_size_1 = sample_size_1, sample_size_2 = sample_size_2 , 
                                       pop_prevalence = prop_Prevalence, excess_mortality = excess_mortality,
                                       return_boots = FALSE)
    
    meansInc[counter] <- rel$Incidence_PE[1]
    std_errors[counter] <-rel$Incidence_SE[1]
    bias[counter] <- Incidence(TT)- meansInc[counter]
    RMSE[counter] <- sqrt((std_errors[counter])^2+(bias[counter])^2)
    
    if (N_iterations > 1){
      meansInc <- rel$Incidence_PE[2]
      std_errors[counter] <-rel$Incidence_SE[2]
      bias[counter] <- Incidence(TT)- meansInc
      RMSE[counter] <- sqrt((std_errors[counter])^2+(bias[counter])^2)
    }
    counter <- counter + 1
  }
  return(data.frame(Delta_t = Delta_t, Bias= abs((bias/Incidence(TT))*100), RSE = (std_errors/Incidence(TT))*100, RMSE = ((RMSE/Incidence(TT))*100)))
}



DATA0 <- RMSE_STDE_BIAS_Cal_1(N_iterations = 1,
                              TT = 57.5, 
                              Incidence = incidence_dem, 
                              min_dT = 1, 
                              Max_dT = 6, 
                              excess_mortality = mortality_dem_excess(57.5))





dat0 <- melt(DATA0, id.vars = "Delta_t")

ggplot(dat0) + geom_smooth(aes(x=Delta_t, y=value, colour= variable), size  = 1) +
  scale_colour_manual(values=c("red","green","blue"))+
  labs(x = "Inter-Survey Interval (Years)", y = "Percent", color  = "")+
  theme_bw(base_size = 18, base_family = "") 





