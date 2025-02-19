#10/10/2025
#Methodology for statistical method using longitudinal data
#Useful functions for sleep data example codes

library(dplyr) #for data manipulation

#--------------------Data set---------------------------------------------------
#' Simulation data for CPAP adherence (Normal positive data)
#'
#' @param nb_patient  integer, number of patients included in the study
#' @param nb_time_point integer, number of time points
#'
#' @returns Normal positive data set simulating CPAP adherence
sim_data <- function(nb_patient, nb_time_point){
  sim_df <- matrix(nrow = nb_patient, ncol = nb_time_point)
  
  for(i in 1:nb_patient){
    seed <- i+31
    set.seed(seed)
    sim_df[i,] <- rnorm(nb_time_point, 4, 1.5)
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
#' @param nb_time_point integer, number of time points
#' @param score_max integer, the maximum number that the variable can take
#'
#' @returns Binomial data set simulating a categorical variable with score_max 
#' possible values
sim_data_discrete <- function(nb_patient, nb_time_point, score_max){
  sim_df <- matrix(nrow = nb_patient, ncol = nb_time_point)
  
  for(i in 1:nb_patient){
    seed <- i+31
    set.seed(seed)
    sim_df[i,] <- sample(0:score_max, nb_time_point, replace = T)
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
#' @returns Barplots, for each time points, of the probability to be in a cluster
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