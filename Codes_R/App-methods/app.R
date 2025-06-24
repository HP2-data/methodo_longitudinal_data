#23/06/2025
#Longitudinal data, trajectories and time series: how to analyze them?
#An example of sleep data

#--------------------------Packages---------------------------------------------
library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(kml)
library(DT)
library(ggplot2)
library(rstatix)
library(Hmisc)
library(LMest)
#devtools::install_github("gitedric/trajeR")
library(trajeR)
library(lcmm)
library(lme4)
library(survival)
library(survminer)
library(forecast)
library(depmixS4)
select <- dplyr::select

source("../functions.R")

#-------------------------Data set/data preparation-----------------------------
load('../../Data/300Patients_90TimePoints/R_shiny_app_data/Sim_CPAP.RData')
load('../../Data/300Patients_90TimePoints/R_shiny_app_data/Sim_CPAP_cat.RData')
load('../../Data/300Patients_90TimePoints/R_shiny_app_data/Sim_ESS_cat.RData')
load('../../Data/300Patients_90TimePoints/R_shiny_app_data/Sim_ESS.RData')

Sim_CPAP <- read.csv("C:/Users/HP2/Desktop/Methodo_stat_donnees/Data/300Patients_90TimePoints/Sim_CPAP.csv",
                     sep=";")
clus_KML <- cld(Sim_CPAP, timeInData = 2:91, maxNA = 1, idAll = Sim_CPAP$patient_id)


#----------------------------R shiny application--------------------------------
#Prepare the page setting using sidebar, accordion results, etc.
###ui
ui <- page_navbar(
  title = 'Statistical methods to analyse longitudinal data',
  
  sidebar = sidebar(
    sliderInput(inputId = 'TimePoint',
                label = 'ANOVA/Chi2 - Time points:',
                min = 1,
                max = 90,
                value = 1),
    radioButtons(inputId = 'TimePoint2',
              label = 'Chi2 - Time points:',
              c('T1' = 'T1', 'T2' = 'T2')),
    sliderInput(inputId = 'nbClusters',
              label = 'Kmeans/LTA/GBTM/GMM - Number of clusters:',
              min = 2,
              max = 5,
              value = 3),
    sliderInput(inputId = 'Degree',
                label = 'GBTM - Degree of the curve:',
                min = 1,
                max = 3,
                value = 2)),
  
  nav_panel(
    title = 'ANOVA methods',
    accordion_panel("Table: Comparison of CPAP adherence at different time points",
                    DTOutput(outputId = 'ANOVA_table')),
    accordion_panel('Plot: Normal values assumption',
                    plotOutput(outputId = 'ANOVA_plot'))),
  nav_panel(
    title = 'Chi² method',
    "Comparison of CPAP adherence and ESS score at 2 time points",
    verbatimTextOutput(outputId = 'chi2_print')),
  nav_panel(
    title = 'Kmeans method',
    accordion_panel("Table: Patients' cluster membership",
                    DTOutput(outputId = 'clus_table')),
    accordion_panel('Plot: Number of patients per cluster ',
                    plotOutput(outputId = 'clus_plot'))),
  nav_panel(title = 'LTA method',
            accordion_panel('Table: LTA method parameters',
                            verbatimTextOutput(outputId = 'LTA_table')),
            accordion_panel('Plot: Transition probabilities',
            plotOutput(outputId = 'LTA_plot'))),
  nav_panel(title = 'GBTM method',
            accordion_panel("Table: Choice of method's parameters",
                            verbatimTextOutput(outputId = 'GBTM_table')),
            accordion_panel('Plot: Clusters trajectories',
                            plotOutput(outputId = 'GBTM_plot')),
            accordion_panel('Plot: Number of patients in each cluster',
                            plotOutput(outputId = 'GBTM_plot2'))),
  nav_panel(title = 'GMM method',
            'Table: GMM method parameters',
            verbatimTextOutput(outputId = 'GMM_table')),
  nav_panel(title = 'Mixed method',
            accordion_panel('Table: linear mixed model',
                            verbatimTextOutput(outputId = 'Mixed_table')),
            accordion_panel('Plot: Association between outcome and factors',
                            plotOutput(outputId = 'Mixed_plot'))),
  nav_panel(title = 'Survival method',
              accordion_panel('Plot: Survival probabilitiy to use CPAP over the study',
                              plotOutput((outputId = 'Survival_plot'))),
              accordion_panel("Plot: Survival probability to use CPAP over the study according to patients' sleepiness status",
                              plotOutput(outputId = 'Survival_plot2')),
              accordion_panel('Table: Comparison test for the use of CPAP over the study acording to groups of sleepiness patients',
                              verbatimTextOutput(outputId = 'Survival_table'))),
  nav_panel(title = 'ARIMA and Cross-correlation method',
            accordion_panel('Plot: ACF, PACF and time series for CPAP adherence',
                            plotOutput(outputId = 'ARIMA_plotCPAP')),
            accordion_panel('Plot: ACF, PACF and time series for ESS score',
                            plotOutput(outputId = 'ARIMA_plotESS')),
            accordion_panel('Table: ARIMA parameters for CPAP adherence',
                            verbatimTextOutput(outputId = 'ARIMA_tableCPAP')),
            accordion_panel('Table: ARIMA parameters for ESS score',
                            verbatimTextOutput(outputId = 'ARIMA_tableESS')),
            accordion_panel('Plot: QQplot to test normality distribution of CPAP adherence values',
                            plotOutput(outputId = 'ARIMA_plot2CPAP')),
            accordion_panel('Plot: ACF for ESS score values',
                            plotOutput(outputId = 'ARIMA_plot2ESS')),
            accordion_panel('Table: Ljung-Box test for CPAP adherence',
                            verbatimTextOutput(outputId = 'ARIMA_table2CPAP')),
            accordion_panel('Table: Ljung-Box test for ESS score',
                            verbatimTextOutput(outputId = 'ARIMA_table2ESS')),
            accordion_panel('Plot: Lags',
                            plotOutput(outputId = 'CCF_plot')),
            accordion_panel('Table: Linear regression association',
                            verbatimTextOutput(outputId = 'CCF_table'))),
  nav_panel(title = 'HMM method',
            accordion_panel('Table: HMM parameters',
                            verbatimTextOutput(outputId = 'HMM_table')),
            accordion_panel('Plot: State trajectory over time',
                            plotOutput(outputId = 'HMM_plot')),
            accordion_panel('Plot: Number of time points spent in each state',
                            plotOutput(outputId = 'HMM_plot2')))
)

#Functions to be executed to obtain the results of each statistical method
###server
server <- function(input, output){
  
  ##ANOVA method
    #Data set
  ANOVA_df <- Sim_CPAP %>%
    pivot_longer(cols = T1:T90, names_to = 'time', values_to = 'Adherence') %>%
    mutate_at(vars(patient_id, time), as.factor)
  
  output$ANOVA_table <- renderDT({
    #ANOVA test
    get_anova_table(anova_test(data = ANOVA_df, dv = Adherence, wid = patient_id,
                               within = time))
  })
  
  output$ANOVA_plot <- renderPlot({
    #Normal values assumption
    TimePoint <- paste0('T', input$TimePoint)
    ggpubr::ggqqplot(ANOVA_df[ANOVA_df$time == TimePoint,], 'Adherence',
                     title = paste0('Time point ', TimePoint))
  })
  
  ##Chi² method
  output$chi2_print <- renderPrint({
      #Comparison test between the ESS score (with 2 available questionnaires) and
      #CPAP adherence (the 1st two time points)
    TimePoint <- paste0('T', input$TimePoint)
    TimePoint2 <- input$TimePoint2
    
    if(TimePoint %nin% c('T1', 'T2')){
           print('Time point must be 1 or 2 !')}
       else{
      #Contingency table
    tableTimePoint1 <- t(table(as.matrix(select(Sim_ESS_cat, TimePoint)),
                             as.matrix(select(Sim_CPAP_cat, TimePoint))))
    
    tableTimePoint2 <- t(table(as.matrix(select(Sim_ESS_cat, TimePoint2)),
                               as.matrix(select(Sim_CPAP_cat, TimePoint2))))
    
    Chi2_test <- array(c(tableTimePoint1, tableTimePoint2),
          dim = c(3, 2, 2),
          dimnames = list(Adherence = c('[0h, 2h[', '[2h, 4h[', '\u2265 4h'),
                          ESS_score = c('No', 'Yes'),Time = c('T1', 'T2')))
    
    print(paste0('The pvalue of the Mantelhaen test is ',
          round(mantelhaen.test(Chi2_test, correct = F)$p.value, 3)))
       }
  }) 
  
  ##Kmeans method
  output$clus_table <- renderDT({
    #Number of clusters
    nbClusters <- input$nbClusters
    
    #K-means function
    kml(clus_KML, nbClusters = nbClusters, nbRedrawing = 15)
  })
  
  output$clus_plot <- renderPlot({
    #Number of clusters
    nbClusters <- input$nbClusters
    
    #Data set
    Sim_CPAP %>%
      mutate('CPAP adherence clusters' = getClusters(clus_KML, nbCluster = nbClusters)) %>%
      select(patient_id, 'CPAP adherence clusters') %>%
      rename(Patient = patient_id)
    
    #K-means function
    kml(clus_KML, nbClusters = nbClusters, nbRedrawing = 15)
    
    #Number of patients per cluster
    ggplot(Sim_CPAP, aes(x = Clusters_obs)) +
      geom_bar(stat = 'count') +
      theme_classic() +
      xlab('Clusters') + 
      ylab('Number of patients')
  })
  
  ##LTA method
    #Data set
  Sim_CPAP_LTA <- Sim_CPAP_cat %>%
    select(patient_id, T1:T90) %>%
    pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'CPAP_adherence') %>%
    mutate_at(vars(patient_id, Time, CPAP_adherence), as.factor) %>%
    separate(Time, sep = 1, into = c('T', 'Time')) %>%
    select(-'T') %>%
    mutate_all(as.numeric) %>%
    as.data.frame()
  
  output$LTA_table <- renderPrint({
    #Number of patients
    nbClusters <- input$nbClusters
    
    #LTA function
    model <- lmest(data = Sim_CPAP_LTA, index = c('patient_id', 'Time'),
                   start = 1, seed = 123, maxit = 5000, modBasic = 1, k = nbClusters)
    
    summary(model)
  })

  output$LTA_plot <- renderPlot({
    #Number of clusters
    nbClusters <- input$nbClusters
    
    #LTA function
    model <- lmest(data = Sim_CPAP_LTA, index = c('patient_id', 'Time'),
                   start = 1, seed = 123, maxit = 5000, modBasic = 1, k = nbClusters)
    
    #Transitions probability
    plot(model, what = "transitions")
  })
  
  ##GBTM method
    #Data set
  Sim_ESS_GBTM <- Sim_ESS %>%
    rename(ESS_baseline = T1) %>%
    mutate_at(vars(ESS_baseline, patient_id), as.factor)
  
  Sim_CPAP_GBTM <- Sim_CPAP %>%
    mutate_at(vars(patient_id), as.factor) %>%
    select(patient_id, T1:T5) %>%
    mutate(Time1 = 1, Time2 = 2, Time3 = 3, Time4 = 4, Time5 = 5) %>%
    left_join(select(Sim_ESS_GBTM, patient_id, ESS_baseline))
 
  output$GBTM_table <- renderPrint({
    #Number of clusters
    nbClusters <- input$nbClusters
    
    #Degree of trajectory curve
    Degree <- input$Degree
    
    #GBTM function
    GBTM_test <- trajeR(Sim_CPAP_GBTM[,2:6], Sim_CPAP_GBTM[,7:11], ng = nbClusters,
                        degre = rep(Degree, nbClusters), Model = 'CNORM')
    
    #Assumptions verification
    print(paste0('BIC: ',round(trajeRBIC(GBTM_test), 2),
                 ' and AIC: ', round(trajeRAIC(GBTM_test), 2)))
    
    adequacy(GBTM_test, Y = Sim_CPAP_GBTM[,2:6], A = Sim_CPAP_GBTM[,7:11])
  
  })
  
  output$GBTM_plot <- renderPlot({
    #Number of clusters
    nbClusters <- input$nbClusters
    
    #Degree of trajectory curve
    Degree <- input$Degree
    
    #GBTM function
    GBTM_test <- trajeR(Sim_CPAP_GBTM[,2:6], Sim_CPAP_GBTM[,7:11], ng = nbClusters,
                        degre = rep(Degree, nbClusters), Model = 'CNORM')
    
    #Trajectories
    plotrajeR(GBTM_test)
  
    })
  
    output$GBTM_plot2 <- renderPlot({
      #Number of clusters
      nbClusters <- input$nbClusters
      
      #Degree of trajectory curve
      Degree <- input$Degree
      
      #GBTM function
      GBTM_test <- trajeR(Sim_CPAP_GBTM[,2:6], Sim_CPAP_GBTM[,7:11], ng = nbClusters,
                          degre = rep(Degree, nbClusters), Model = 'CNORM')
      
      #Cluster to which the patient belongs
      Group_N <- GroupProb(GBTM_test, Y = Sim_CPAP_GBTM[,2:6], A = Sim_CPAP_GBTM[,7:11])

      Group_N <- Group_N %>%
        as.data.frame() %>%
        rowwise() %>%
        mutate(group = which.max(c_across(Gr1:paste0('Gr',nbClusters)))) %>%
        mutate_at(vars(group), as.factor)
    
      #Numer of patients per cluster
      ggplot(Group_N, aes(x = group)) +
        geom_bar(stat = 'count') +
        theme_classic() +
        xlab('Clusters') + 
        ylab('Number of patients')
    })
    
    ##GMM method
      #Data set
    Sim_CPAP_GMM <- Sim_CPAP %>%
      select(patient_id, T1:T5) %>%
      pivot_longer(cols = c(T1:T5), names_to = 'Time', values_to = 'CPAP_adherence') %>%
      mutate_at(vars(patient_id), as.numeric)
    
    output$GMM_table <- renderPrint({
      #Number of clusters
      nbClusters <- input$nbClusters
    
      #GMM function
      set.seed(123)
      GMM_test1 <- hlme(CPAP_adherence ~ Time, subject = 'patient_id',
                        random = ~ 1 + Time, ng = 1, data = Sim_CPAP_GMM)
      GMM_test <-  gridsearch(rep = 100, maxiter = 10, minit = GMM_test1,
                              hlme(CPAP_adherence ~ Time, subject = 'patient_id',
                                   random = ~ 1 + Time, ng = nbClusters, data = Sim_CPAP_GMM,
                                   mixture = ~ Time, nwg = T))
      
      #Parameters for choosing the best model
      summarytable(GMM_test1, GMM_test)
      
      #A posterior probability of belonging to each cluster
      postprob(GMM_test)
      
      #Model results
      summary(GMM_test)
    })
    
    ##Mixed model
      #Data set
    Sim_CPAP_mixed <- Sim_CPAP %>%
      pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'CPAP_adherence') %>%
      mutate_at(vars(patient_id, Time), as.factor) %>%
      separate(Time, sep = 1, into = c('T', 'Time')) %>%
      select(-'T') %>%
      mutate_at(vars(Time), ~ factor(.x, levels = sort(unique(as.numeric(.x))))) %>% #ordered values
      mutate_at(vars(Time), as.numeric) 
    
    Sim_ESS_mixed <- Sim_ESS_cat %>%
      select(T1, patient_id) %>%
      mutate_at(vars(patient_id, T1), as.factor)
    
    Sim_CPAP_mixed <- Sim_CPAP_mixed %>%
      left_join(Sim_ESS_mixed, by = 'patient_id') %>%
      rename(ESS_baseline = T1)
    
    #Linear mixed regression
    Mixed_test <- lmer(CPAP_adherence ~ Time + ESS_baseline + (1 + Time | patient_id),
                       data = Sim_CPAP_mixed)
    
    output$Mixed_table <- renderPrint({
      summary(Mixed_test)
    })
    
    output$Mixed_plot <- renderPlot({
      sjPlot::plot_model(Mixed_test)
    })
    
    ##Survival method
      #Data set
    Sim_CPAP_survival <- Sim_CPAP %>%
      pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'CPAP_adherence') %>%
      mutate_at(vars(patient_id), as.factor) %>%
      mutate_at(vars(Time), ~as.numeric(as.factor(.x))) %>%
      mutate_at(vars(CPAP_adherence), ~ifelse(.x == 0, 1, 0)) %>%
      group_by(patient_id) %>%
      summarise(Time = ifelse(max(CPAP_adherence) == 1, which.max(CPAP_adherence), max(Time)),
                CPAP_adherence = max(CPAP_adherence))  
    
    Sim_ESS_joint <- Sim_ESS_cat %>%
      select(patient_id, T1) %>%
      rename(Drowsy = T1) %>%
      mutate_at(vars(patient_id), as.factor) 
    
    Sim_CPAP_survival <- Sim_CPAP_survival %>%
      left_join(Sim_ESS_joint)
    
    output$Survival_plot <- renderPlot({
      #Survival function
      fit <- survfit(Surv(Time, CPAP_adherence) ~ 1, data = Sim_CPAP_survival)
      ggsurvplot(fit, risk.table = T)
      })

    output$Survival_plot2 <- renderPlot({
      #Survival function / according to sleepiness status of patients
      fit <- survfit(Surv(Time, CPAP_adherence) ~ Drowsy, data = Sim_CPAP_survival)
      ggsurvplot(fit, risk.table = T)
    
    })
    
    output$Survival_table <- renderPrint({
      #Comparison test of CPAP use between drowsy and non-drowsy patients
      survdiff(Surv(Time, CPAP_adherence) ~ Drowsy, data = Sim_CPAP_survival)
      })
    
    ##ARIMA and Cross-correlation method
      #Data set
    Sim_CPAP_ARIMA <- Sim_CPAP %>%
      select(patient_id, T1:T90) %>%
      filter(patient_id == 10) %>%
      pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'CPAP_adherence') %>%
      mutate_at(vars(patient_id), as.factor) %>%
      separate(Time, sep = 1, into = c('T', 'Time')) %>%
      select(-c('T', 'patient_id'))
    
    Sim_ESS_ARIMA <- sim_data_discrete(300, 90, 24) %>%
      select(patient_id, T1:T90) %>%
      filter(patient_id == 10) %>%
      pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'ESS_score') %>%
      mutate_at(vars(patient_id), as.factor) %>%
      separate(Time, sep = 1, into = c('T', 'Time')) %>%
      select(-c('T', 'patient_id'))
    
    #Time series creation
    ts_CPAP <- ts(Sim_CPAP_ARIMA$CPAP_adherence)
    ts_ESS <- ts(Sim_ESS_ARIMA$ESS_score)
    
    #ARIMA function 
    ARIMA_CPAP <- auto.arima(ts_CPAP)
    ARIMA_ESS <- auto.arima(ts_ESS)
    
    ESS_res <- ARIMA_ESS$residuals
    CPAP_res <- ARIMA_CPAP$residuals
    
    #Cross-correlation function
    CCF_CPAP_ESS <- ccf(CPAP_res, ESS_res)
    
    output$ARIMA_plotCPAP <- renderPlot({
      #Time series description
      fUnitRoots::urkpssTest(ts_CPAP, type = c("tau"), lags = c("short"),
                             use.lag = NULL, doplot = TRUE)
    })
    
    output$ARIMA_plotESS <- renderPlot({
      #Time series description
      fUnitRoots::urkpssTest(ts_ESS, type = c("tau"), lags = c("short"),
                             use.lag = NULL, doplot = TRUE)
    })
    
    output$ARIMA_tableCPAP <- renderPrint({
      #ARIMA parameters
      summary(ARIMA_CPAP)
    })
    
    output$ARIMA_tableESS <- renderPrint({ 
      #ARIMA parameters
      summary(ARIMA_ESS) 
    })
    
    output$ARIMA_plot2CPAP <- renderPlot({
      #Normal values assumption
      qqnorm(ARIMA_CPAP$residuals)
      qqline(ARIMA_CPAP$residuals)
    })
    
    output$ARIMA_plot2ESS <- renderPlot({
      #ACF function
      acf(ARIMA_ESS$residuals)
    })
    
    output$ARIMA_table2CPAP <-renderPrint({
      #Test Ljung-Box to check residuals patterns
      Box.test(ARIMA_CPAP$residuals, type = 'Ljung-Box')
    })
    
    output$ARIMA_table2ESS <-renderPrint({
      #Test Ljung-Box to check residuals patterns
      Box.test(ARIMA_ESS$residuals, type = 'Ljung-Box')
    })
    
    output$CCF_plot <- renderPlot({
      #ACF of the cross-correlation 
      plot(CCF_CPAP_ESS)
      })

    output$CCF_table <- renderPrint({
      #Linear regression using the results of cross-correlation method lags
      lag <- stats::lag
      final_data <- ts.intersect(ESS_res, CPAPlag14 = lag(CPAP_res, 14))
    
      reg_CCF <- lm(ESS_res ~ CPAPlag14, data = final_data)
      summary(reg_CCF)
      })
    
    
    ##Hiden Markov Model
      #Data set
    Sim_CPAP_HMM <- Sim_CPAP_cat %>%
      filter(patient_id %in% seq(100, 300)) %>%
      mutate_at(vars('T1':'T90'), ~ as.numeric(as.factor(.x))) %>%
      pivot_longer(cols = c('T1':'T90'), names_to = 'Time', values_to = 'CPAP_adherence')

      #To take into account patient_id: the lengths of individual, i.e. independent, time series.
    ntimes <- Sim_CPAP_HMM %>%
      group_by(patient_id) %>%
      summarise(N = length(CPAP_adherence))
    
    Sim_CPAP_HMM <- Sim_CPAP_HMM %>%
      arrange(patient_id)

      #HMM function
    HMM_test <- depmix(response = CPAP_adherence ~ 1, family = multinomial(), nstates = 2,
                       data = Sim_CPAP_HMM, ntimes = ntimes$N)
    
    set.seed(123)
    HMM_final <- fit(HMM_test)
  
    output$HMM_table <- renderPrint({
      #Statistical verification of final states results
      print(paste0('BIC: ', round(BIC(HMM_final), 2), ', AIC: ', round(AIC(HMM_final),2),
             ' and Loglikelihood: ', round(logLik(HMM_final), 2)))
      
      #Initial state probability
      init_prob <- matrix(getpars(HMM_final)[1:2], nrow = 1, byrow = T)
      colnames(init_prob) <- c('Non adherent', 'Adherent')
      print('Initial state probability')
      print(init_prob)
      
      #Transition matrix
      trans_prob <- matrix(getpars(HMM_final)[3:6], nrow = 2, byrow = T)
      colnames(trans_prob) <- c('Non adherent', 'Adherent')
      rownames(trans_prob) <- c('Non adherent', 'Adherent')
      print('Transition matrix')
      print(trans_prob)
      
      #Predict the hidden states up to 200 time points
      pred_states <- posterior(HMM_final, type = 'viterbi')[1:200,]
      print('Predict the hidden states up to 200 time points')
      print(pred_states)
    })
  
    #Predict the hidden states up to 200 time points
    pred_states <- posterior(HMM_final, type = 'viterbi')[1:200,]
    
    output$HMM_plot <- renderPlot({
      #Change of states during the study
      ggplot(pred_states, aes(y = state, x = seq(1, length(state)))) +
        geom_line() +
        ylab('State prediction') +
        xlab('Time points') +
        theme_classic() + 
        scale_y_discrete(breaks = c(1, 2), limits = factor(c(1, 2)))
    })
    
    output$HMM_plot2 <- renderPlot({
      #Number of time points spent in each state
      pred_states <- pred_states %>%
        mutate_at(vars(state), as.factor)
      
      ggplot(pred_states, aes(x = state)) +
        geom_bar(stat = 'count') + 
        theme_classic()
    })
    
}

#Run the app
###Application
shinyApp(ui = ui, server = server)
