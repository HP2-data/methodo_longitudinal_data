#03/04/2025
#Longitudinal data, trajectories and time series: how to analyze them?
#An example of sleep data

library(dplyr) #for data manipulation

#--------------------Data set---------------------------------------------------
#' Simulation data for CPAP adherence (Normal positive data)
#'
#' @param nb_patient  integer, number of patients included in the study
#' @param nb_time_point integer, number of measuring points
#'
#' @returns Normal positive data set simulating CPAP adherence
sim_data <- function(nb_patient, nb_time_point){
  sim_df <- matrix(nrow = nb_patient, ncol = nb_time_point)
  #3 trajectories
  group <- list()
  group[[1]] <- seq(1,nb_patient/3)
  group[[2]] <- seq(nb_patient/3 + 1, (nb_patient/3)*2)
  group[[3]] <- seq(((nb_patient/3)*2)+1, nb_patient)
  
  #Varying depending on 3 time windows
  time_point <- nb_time_point/3
  time <- list()
  time[[1]] <- seq(1, time_point)
  time[[2]] <- seq(time_point + 1, time_point*2)
  time[[3]] <- seq((time_point*2) + 1, nb_time_point)
  
  j <- 1
  k <- 2
  l <- 3
  
  while(time_point < (nb_time_point + 1)){
    for(i in group[[j]]){
      seed <- i+31
      set.seed(seed)
      sim_df[i,time[[j]]] <- rnorm(length(time[[j]]), 2.2, 1)
      }
    for(i in group[[k]]){
      seed <- i+31
      set.seed(seed)
      sim_df[i,time[[j]]] <- rnorm(length(time[[j]]), 4.8, 0.6)
      }
    for(i in group[[l]]){
      seed <- i+30
      set.seed(seed)
      sim_df[i,time[[j]]] <- rnorm(length(time[[j]]), 6.3, 0.7)
    }
    
    j <- k
    k <- l
    l <- j -1
    time_point <- time_point + nb_time_point/3
  }
  
  sim_df <- sim_df %>%
    as.data.frame() %>%
    mutate(patient_id = seq(1,nb_patient,1))
  
  names(sim_df) <- c(paste0('T', seq(1, nb_time_point, 1)), 'patient_id')
  return(sim_df)
}


#' Simulation data for ESS score (categorical values from 0 to 24)
#'
#' @param nb_patient  integer, number of patients included in the study
#' @param nb_time_point integer, number of measuring points
#' @param score_max integer, the maximum number that the variable can take
#'
#' @returns Beta binomial data set simulating a categorical variable with score_max 
#' possible values
sim_data_discrete <- function(nb_patient, nb_time_point, score_max){
  sim_df <- matrix(nrow = nb_patient, ncol = nb_time_point)
  
  for(i in 1:nb_patient){
    seed <- i+31
    set.seed(seed)
    sim_df[i,] <- rbinom(nb_time_point, score_max, rbeta(nb_time_point, 9, 15))
  }
  sim_df <- sim_df %>%
    as.data.frame() %>%
    mutate(patient_id = seq(1,nb_patient,1))
  
  names(sim_df) <- c(paste0('T', seq(1, nb_time_point, 1)), 'patient_id')
  return(sim_df)
}

#--------------------LCA method-------------------------------------------------
#source: https://www.geeksforgeeks.org/latent-class-analysis-in-r/

#' Plot LCA model
#'
#' @param lca_model  LCA model fitted
#' @returns Barplots, for each time points, of the probability of belonging to each cluster
plot_lca <- function(lca_model) {
  probs <- lca_model$probs
  num_classes <- length(probs)
  par(mfrow = c(1, num_classes))
  
  for (i in 1:num_classes) {
    barplot(t(as.matrix(probs[[i]])), beside = TRUE, col = rainbow(ncol(probs[[i]])),
            main = paste("T", i), xlab = "CPAP adherence", ylab = "Probability")
  }
  legend("topright", legend = c('[0h,2h[', '[2h,4h[', '\u2265 4h'),
         fill = c('red', 'green', 'blue'))
}
