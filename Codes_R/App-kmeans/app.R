library(shiny)
library(bslib)
library(dplyr)
library(kml)
library(DT)
library(ggplot2)


#Data set/data preparation
load('../../Data/300Patients_90TimePoints/Sim_CPAP.RData')


clus_KML <- cld(Sim_CPAP, timeInData = 2:91, maxNA = 1, idAll = Sim_CPAP$patient_id)


ui <- page_sidebar(
  title = 'Clustering by K-means method',
  
  sidebar = sidebar(
    sliderInput(inputId = 'nbClusters',
                label = 'Number of clusters:',
                min = 2, 
                max = 5,
                value = 3)),
  accordion(
    open = c('Table'),
    accordion_panel("Table: patients' cluster membership",
                    DTOutput(outputId = 'clus_table')
    ),
    accordion_panel('Plot: Number of patients per cluster ',
                    plotOutput(outputId = 'clus_plot'))
    )
)


server <- function(input, output){
  
  output$clus_table <- renderDT({
    nbClusters <- input$nbClusters
    
    kml(clus_KML, nbClusters = nbClusters, nbRedrawing = 15)
    
    Sim_CPAP %>%
      mutate('CPAP adherence clusters' = getClusters(clus_KML, nbCluster = nbClusters)) %>%
      select(patient_id, 'CPAP adherence clusters') %>%
      rename(patient = patient_id)
    
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
}


shinyApp(ui = ui, server = server)

