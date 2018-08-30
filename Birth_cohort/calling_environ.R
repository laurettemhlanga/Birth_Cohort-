

#set the working directory to where project "Birth_cohort" is  saved  and source the script "Revised_birthcohort.R"


source('file:///H:/Birth_cohort/Birth_cohort/Revised_birthcohort.R', echo = FALSE)


#birth cohort 
y = Sim_Birth_Cohort(Spline = FALSE)

#prevalence plot 
dd <- ggplot(y) + geom_line(aes(x=time, y = Prevalence), size  = 1)+
  labs(x = "Time", y = "Prevalence")+
  theme(panel.grid.major = element_blank(),text = element_text(size=18),
        panel.grid.minor = element_blank(),legend.text=element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



#Extracting the spline functions for each attribute ofinterest
#when Recency is true
Prop_Susceptible <- Sim_Birth_Cohort()$Susceptible_Function
Prop_Prev_Recency <-  Sim_Birth_Cohort()$Prev_Recency_Function
Prop_NOn_Recency <- Sim_Birth_Cohort()$Prev_Non_Receny_Function  
Prop_Prevalence <- Sim_Birth_Cohort()$Prevalence_Function
#when Recency false 
Prop_Prevalence <- Sim_Birth_Cohort()$Prevalence_Function
Prop_NOn_Recency <- Sim_Birth_Cohort()$Prev_Non_Receny_Function    


# comparison for Recency estimator VS Mahiane estimator.



#estimating Incidence recency
sigma_Prop_Prevalence <-  sqrt(Prop_Prevalence(20)*(1-Prop_Prevalence(20))/10000)
sigma_Prev_Prop_Recency <-  sqrt(Prop_Prev_Recency(20)*(1-Prop_Prev_Recency(20))/(Prop_Prevalence(20)*10000))
Estimates <- incprops(PrevH = Prop_Prevalence(20), RSE_PrevH = sigma_Prop_Prevalence,
                      PrevR = Prop_Prev_Recency(20), RSE_PrevR =  sigma_Prev_Prop_Recency,
                      BS_Count = 100000, Boot = FALSE, MDRI = 170.4431, RSE_MDRI = 0.05,
                      FRR = 0.001190825, RSE_FRR = 0, BigT = 730)



 
#incidence estimation using the change prevalence approach 
Inci = Midpoint_inc_from_prevalence(N_iterations=10000, T1 = 17, T2 = 23, 
                                    sample_size_1 = 10000, sample_size_2  = 10000,
                                    pop_prevalence = Prop_Prevalence, 
                                    excess_mortality_estimate = Excess_Mortality_var(20), 
                                    return_boots = TRUE)






ggplot(Inci[c(1:10000),], aes(Incidence_PE)) +
  geom_histogram(bins = 120)+
  geom_vline(xintercept= Incidence_var(20), linetype="dashed", color = "blue", size = 1)+
 labs(x = "Incidence estimates", y = " Estimate (True Inc= 0.02)")+
  theme_bw(base_size = 18, base_family = "")+
  xlim(0,0.023)

  



#############################################################################################################
#############################################################################################################
#############################################################################################################
#investigation of the RMSE as a function of the time interval between surveys.


RMSE_STDE_BIAS_Cal_1 <- function( N_iterations, min_dT, Max_dT,
                                  TT, sample_size_1 = 10000, Incidence, sample_size_2 = 10000, 
                                  pop_prevalence = Prop_Prevalence, excess_mortality = 0.1,
                                  return_boots = FALSE,
                                  Delta_t = seq(2, 12,2)){
  survey_intervals = seq(min_dT, Max_dT, 1) 
  N_intervals <- length(survey_intervals)
  #TT = 23
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
                                       pop_prevalence = Prop_Prevalence, excess_mortality = excess_mortality,
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
  return(data.frame(Delta_t = Delta_t, Bias= (bias/Incidence(TT))*100, RSE = (std_errors/Incidence(TT))*100, RMSE = ((RMSE/Incidence(TT))*100)))
}



DATA0 <- RMSE_STDE_BIAS_Cal_1(N_iterations = 1 ,TT = 23, Incidence = Incidence_var, min_dT = 1, Max_dT = 6, excess_mortality = 0.1)





dat0 <- melt(DATA0, id.vars = "Delta_t")

ggplot(dat0) + geom_smooth(aes(x=Delta_t, y=value, colour= variable), size  = 1) +
  scale_colour_manual(values=c("red","green","blue"))+
  geom_hline(mapping=aes(yintercept=0.0511*100, colour = "RMSE_K"), linetype="longdash", color = "green",size = 1)+
  geom_hline(yintercept= 1.52, linetype="dashed", color = "red",size = 1)+
  geom_hline(yintercept= 5.331276, linetype="dashed", color = "blue", size = 1)+
  labs(x = "Inter-Survey Interval (Years)", y = "Percent", color  = "")+
  theme_bw(base_size = 18, base_family = "") 





