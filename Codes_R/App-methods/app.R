library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(kml)
library(DT)
library(ggplot2)
library(rstatix)
library(Hmisc)


#Data set/data preparation
load('../../Data/300Patients_90TimePoints/Sim_CPAP.RData')
load('../../Data/300Patients_90TimePoints/ANOVA_df.RData')
load('../../Data/300Patients_90TimePoints/Sim_CPAP_cat.RData')
load('../../Data/300Patients_90TimePoints/Sim_ESS_cat.RData')

clus_KML <- cld(Sim_CPAP, timeInData = 2:91, maxNA = 1, idAll = Sim_CPAP$patient_id)


ui <- page_navbar(
  title = 'Statistical methods to analyse longitudinal data',
  
  sidebar = sidebar(
    sliderInput(inputId = 'nbClusters',
                label = 'Kmeans - Number of clusters:',
                min = 2,
                max = 5,
                value = 3),
    sliderInput(inputId = 'TimePoint',
                label = 'ANOVA/Chi2 - Time points:',
                min = 1,
                max = 90,
                value = 1),
    radioButtons(inputId = 'TimePoint2',
              label = 'Chi2 - Time points:',
              c('T1' = 'T1', 'T2' = 'T2'))),

  nav_panel(
    title = 'Kmeans method',
    accordion_panel("Table: Patients' cluster membership",
                    DTOutput(outputId = 'clus_table')),
    accordion_panel('Plot: Number of patients per cluster ',
                    plotOutput(outputId = 'clus_plot'))),
  nav_panel(
    title = 'ANOVA methods',
    accordion_panel("Table: Comparison of CPAP adherence at different time points",
                    DTOutput(outputId = 'ANOVA_table')),
    accordion_panel('Plot: Normal values assumption',
                    plotOutput(outputId = 'ANOVA_plot'))),
  nav_panel(
    title = 'Chi² methods',
    accordion_panel("Table: Comparison of CPAP adherence and ESS score at 2 time points",
                    verbatimTextOutput(outputId = 'chi2_print')))
)


server <- function(input, output){
  
  ##Kmeans method
  output$clus_table <- renderDT({
    nbClusters <- input$nbClusters
    
    kml(clus_KML, nbClusters = nbClusters, nbRedrawing = 15)
    
    Sim_CPAP %>%
      mutate('CPAP adherence clusters' = getClusters(clus_KML, nbCluster = nbClusters)) %>%
      select(patient_id, 'CPAP adherence clusters') %>%
      rename(Patient = patient_id)
    
  })
  
  output$clus_plot <- renderPlot({
    nbClusters <- input$nbClusters
    
    kml(clus_KML, nbClusters = nbClusters, nbRedrawing = 15)
    
    Sim_CPAP <- Sim_CPAP %>%
      mutate(Clusters_obs = getClusters(clus_KML, nbCluster = nbClusters)) %>%
      select(patient_id, Clusters_obs)
    
    ggplot(Sim_CPAP, aes(x = Clusters_obs)) +
      geom_bar(stat = 'count') +
      theme_classic() +
      xlab('Clusters') + 
      ylab('Number of patients')
  })
  
  ##ANOVA method
  output$ANOVA_table <- renderDT({
    get_anova_table(anova_test(data = ANOVA_df, dv = Adherence, wid = patient_id,
                               within = time))
  })
  
  output$ANOVA_plot <- renderPlot({
    TimePoint <- paste0('T', input$TimePoint)
    ggpubr::ggqqplot(ANOVA_df[ANOVA_df$time == TimePoint,], 'Adherence',
                     title = paste0('Time point ', TimePoint))
  })
  
  ##Chi² method
  output$chi2_print <- renderPrint({
    
    TimePoint <- paste0('T', input$TimePoint)
    TimePoint2 <- input$TimePoint2
    
    if(TimePoint %nin% c('T1', 'T2')){
           print('Time point must be 1 or 2 !')}
       else{
           
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
}

shinyApp(ui = ui, server = server)
