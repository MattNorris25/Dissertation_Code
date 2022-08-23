#Loading packages
library("survival")
library("survminer")
library("gcookbook")
library("ggplot2")
library("viridis")
library("gridExtra")
library("cowplot")

#Simulating trial data
Trial_data = function(sample_size_control, sample_size_intervention, recruitment_length, lambda_control_start, lambda_control_change, lambda_intervention_start, lambda_intervention_change){
  
  #Total sample size
  N = sample_size_control + sample_size_intervention
  
  #Control arm
  
  #Date of recruitment 
  recruitment_date_control = round(runif(sample_size_control, 0, recruitment_length), 0)
  recruitment_date_control = recruitment_date_control[order(recruitment_date_control, decreasing = FALSE)]
  
  #Time to event 
  time_to_event_control = round(rexp(sample_size_control, rate = (lambda_control_change*recruitment_date_control + lambda_control_start)), 0)
  event_date_control = recruitment_date_control + time_to_event_control
  
  #Intervention arm
  
  #Date of recruitment 
  recruitment_date_intervention = round(runif(sample_size_intervention, 0, recruitment_length), 0)
  recruitment_date_intervention = recruitment_date_intervention[order(recruitment_date_intervention, decreasing = FALSE)]
  
  #Time to event
  time_to_event_intervention = round(rexp(sample_size_intervention, rate = (lambda_intervention_change*recruitment_date_intervention + lambda_intervention_start)), 0)
  event_date_intervention = recruitment_date_intervention + time_to_event_intervention
  
  #Data frame creation
  df = data.frame(1:N, c(recruitment_date_intervention, recruitment_date_control), c(time_to_event_intervention, time_to_event_control), c(event_date_intervention, event_date_control), FALSE, FALSE, FALSE, c(rep("Intervention", sample_size_intervention), rep("Control", sample_size_control)))
  colnames(df) = c("index", "recruitment_date", "time_to_event", "event_date", "right_censored", "loss_to_follow_up", "uncensored", "arm")
  
  #End of trial date
  end_of_trial = round(median(c(event_date_intervention, event_date_control),0))
  
  #Right censoring
  for(i in 1:N){
    if(df[i, 4] > end_of_trial){
      df[i, 4] = end_of_trial
      df[i, 5] = TRUE
      df[i, 3] = end_of_trial - df[i, 2]
    }
  }
  
  #Loss to follow up
  unif = runif(N)
  for(i in 1:N){
    if(df[i, 5] == FALSE && unif[i] < 0.05){
      df[i, 6] = TRUE
      df[i, 4] = df[i, 2] + round(df[i, 3]*runif(1), 0)
      df[i, 3] = df[i, 4] - df[i, 2]
    }
  }
  
  #T/F on if patient is censored or not
  for(i in 1:N){
    if(df[i, 5] == FALSE && df[i, 6] == FALSE){
      df[i, 7] = TRUE
    }
  }
  
  #Removing patients recruited after the end of the trial
  keep = rep(TRUE, N)
  for(i in 1:N){
    if(df[i, 2] > end_of_trial){
      keep[i] = FALSE
    }
  }
  df = df[keep, ]
  
  return(df)
}

#Hazard ratio and trial success function
Hazard_Ratios_and_Trial_Success = function(number_of_simulations, prior_alpha_historical, prior_beta_historical, prior_alpha_intervention, prior_beta_intervention, historical_weight, sample_size_current, recruitment_length_current, lambda_control_start_hist, lambda_control_change, lambda_intervention_start_hist, lambda_intervention_change, sample_size_historical, recruitment_length_historical, start_date_difference, number_of_probability_simulations){
  Separate_Analysis_HR = Pooled_Analysis_HR = Regression_Model_HR = Power_Priors_HR = 0
  Separate_Analysis_Decision = Pooled_Analysis_Decision = Regression_Model_Decision = Power_Priors_Decision = 0
  
  for(i in 1:number_of_simulations){
    #Generating simulated data
    Current_trial = Trial_data(sample_size_control = sample_size_current/2, 
                               sample_size_intervention = sample_size_current/2, 
                               recruitment_length = recruitment_length_current, 
                               lambda_control_start = lambda_control_change*start_date_difference + lambda_control_start_hist, 
                               lambda_control_change = lambda_control_change, 
                               lambda_intervention_start = lambda_intervention_change*start_date_difference + lambda_intervention_start_hist, 
                               lambda_intervention_change = lambda_intervention_change)
    Historical = Trial_data(sample_size_control = sample_size_historical/2, 
                            sample_size_intervention = sample_size_historical/2, 
                            recruitment_length = recruitment_length_historical, 
                            lambda_control_start = lambda_control_start_hist, 
                            lambda_control_change = lambda_control_change, 
                            lambda_intervention_start = lambda_control_start_hist, 
                            lambda_intervention_change = lambda_control_change)
    
    #Using only intervention data from historical trial
    Historical = Historical[!(Historical$arm=="Control"), ]
    
    #Separate analysis HR and decision
    model = coxph(Surv(time_to_event, uncensored) ~ arm, data = Current_trial)
    Separate_Analysis_HR[i] = exp(summary(model)$coefficients[,1])
    Separate_Analysis_Decision[i] = summary(model)$coefficients[,5] < 0.05
    
    #Adjusting dates in current trial
    Current_trial$recruitment_date = Current_trial$recruitment_date + start_date_difference
    Current_trial$event_date = Current_trial$event_date + start_date_difference
    
    #Combining datasets
    Historical[, 8] = "Control"
    Current_trial_Historical_Combined = rbind(Current_trial, Historical)
    Current_trial_Historical_Combined[, 9] = c(rep(1, as.numeric(nrow(Current_trial))), rep(0, as.numeric(nrow(Historical))))
    colnames(Current_trial_Historical_Combined) = c("index", "recruitment_date", "time_to_event", "event_date", "right_censored", "loss_to_follow_up", "censored", "arm", "period")
    
    #Pooled analysis HR and decision
    model = coxph(Surv(time_to_event, censored) ~ arm, data = Current_trial_Historical_Combined)
    Pooled_Analysis_HR[i] = exp(summary(model)$coefficients[,1])
    Pooled_Analysis_Decision[i] = summary(model)$coefficients[,5] < 0.05
    
    #Regression model HR and decision
    model = coxph(Surv(time_to_event, censored) ~ arm + recruitment_date, data = Current_trial_Historical_Combined)
    Regression_Model_HR[i] = exp(summary(model)$coefficients[1,1])
    Regression_Model_Decision[i] = summary(model)$coefficients[1,5] < 0.05
    
    #Setting priors
    alpha_historical = prior_alpha_historical
    beta_historical = prior_beta_historical
    alpha_intervention = prior_alpha_intervention
    beta_intervention = prior_beta_intervention
    
    #Updating prior historical distribution
    alpha_historical = alpha_historical + historical_weight*sum(ifelse(Historical[, 7] == TRUE, 1, 0)) #Adding on number of non-censored patients in historical data
    beta_historical = beta_historical + historical_weight*(sum(Historical[,3])) #Adding on time to events for historical data 
    
    #Updating priors
    #Control Arm
    alpha_historical = alpha_historical + sum(Current_trial[(sum(Current_trial[, 8]  == "Intervention")+1):as.numeric(nrow(Current_trial)), 7]) #Adding on number of non-censored patients in control arm
    beta_historical = beta_historical + sum(Current_trial[(sum(Current_trial[, 8]  == "Intervention")+1):as.numeric(nrow(Current_trial)), 3]) #Adding on time to events for control arm
    
    #Intervention arm
    alpha_intervention = alpha_intervention + sum(Current_trial[1:(sum(Current_trial[, 8]  == "Intervention")), 7]) #Adding on number of non-censored patients in intervention arm
    beta_intervention = beta_intervention + sum(Current_trial[1:(sum(Current_trial[, 8]  == "Intervention")), 3]) #Adding on time to events for intervention arm
    
    #Estimating P(control < intervention given data)
    prob = sum(rgamma(number_of_probability_simulations, alpha_historical, beta_historical) > rgamma(number_of_probability_simulations, alpha_intervention, beta_intervention))/number_of_probability_simulations
    
    #Power priors HR and decision
    Power_Priors_HR[i] = (alpha_intervention*beta_historical)/(alpha_historical*beta_intervention)
    Power_Priors_Decision[i] = prob > 0.975 
    
  }
  #Metrics dataframe
  df = data.frame(Separate_Analysis_HR, ifelse(Separate_Analysis_Decision == 1, TRUE, FALSE), Pooled_Analysis_HR, ifelse(Pooled_Analysis_Decision == 1, TRUE, FALSE), Regression_Model_HR, ifelse(Regression_Model_Decision == 1, TRUE, FALSE), Power_Priors_HR, ifelse(Power_Priors_Decision == 1, TRUE, FALSE))
  colnames(df) = c("Separate Analysis HR", "Separate Analysis Decision", "Pooled Analysis HR", "Pooled Analysis Decision", "Regression Model HR", "Regression Model Decision", "Power Priors HR", "Power Priors Decision")
  return(df)
}

#Measures function
Measures = function(simulation_results, true_hr){
  #Measures
  Type_I_error = c(mean(simulation_results[, 2]), mean(simulation_results[, 4]), mean(simulation_results[, 6]), mean(simulation_results[, 8]))
  Bias = c(mean(simulation_results[, 1] - true_hr), mean(simulation_results[, 3] - true_hr), mean(simulation_results[, 5] - true_hr), mean(simulation_results[, 7] - true_hr))
  Variance = c(var(simulation_results[, 1]), var(simulation_results[, 3]), var(simulation_results[, 5]), var(simulation_results[, 7]))
  Mean_Squared_Error = Variance + Bias^2
  
  #Measures dataframe
  Measures = data.frame(Type_I_error, Bias, Variance, Mean_Squared_Error)
  rownames(Measures) = c("Separate Analysis", "Pooled Analysis", "Cox Proportional Hazards Model", "Power Priors")
  
  return(Measures)
}
