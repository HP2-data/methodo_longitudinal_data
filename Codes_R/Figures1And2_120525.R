#12/05/2025
#Longitudinal data, trajectories and time series: how to analyze them?
#An example of sleep data

library(dplyr) #for data manipulation
library(tidyr) #to arrange data
library(TraMineR) #for chronogram
library(ggplot2) #for plots
library(ggalluvial) #for alluvial plots
library(survminer) #for survival plot
library(survival) #for survival plot
library(cowplot) #to add multiple figures in the same plot
select <- dplyr::select

source("functions.R") #all functions useful for this code

#--------------------Data set---------------------------------------------------
#In this section, functions are available in the 'functions.R' file
#For the simulation, we assume 90 time points (around 3 months) for 300 patients

Sim_CPAP_cat <- read.csv2("C:/Users/HP2/Desktop/Methodo_stat_donnees/Data/300Patients_90TimePoints/Sim_CPAP_cat.csv")
Sim_CPAP <- read.csv("C:/Users/HP2/Desktop/Methodo_stat_donnees/Data/300Patients_90TimePoints/Sim_CPAP.csv", sep=";")

#--------------------Figure 1---------------------------------------------------
#A: Chronogram
#5 time points

Sim_CPAP_chrono <- Sim_CPAP_cat %>%
  select('T1':'T5')

seqdplot(seqdef(Sim_CPAP_chrono), cex.plot=1.3, cex.legend = 1,
         cpal = c('#2171B5', '#6BAED6','#C6DBEF'))

#B: Alluvial plot
#3 time points
Sim_CPAP_alluvial <- Sim_CPAP_cat %>%
  select('T1':'T3') %>%
  mutate_at(vars('T1', 'T2', 'T3'), as.factor) %>%
  table() %>%
  as.data.frame() %>%
  mutate(`Evolution T1-T3` = case_when(`T1` == '[0h,2h[' & (`T3` == '[2h,4h[' | `T3` == '\u2265 4h') |
                                          `T1` == '[2h,4h[' & `T3` == '\u2265 4h' ~ 'Improvement',
                                       `T1` == `T3` ~ 'Stable',
                                       `T1` == '\u2265 4h' & (`T3` == '[2h,4h[' |`T3` == '[0h,2h[') |
                                         `T1` == '[2h,4h[' & `T3` ==  '[0h,2h[' ~ 'Decline')) %>%
  mutate_at(vars(`Evolution T1-T3`), as.factor) %>%
  mutate_at(vars(`Evolution T1-T3`), ~ factor(.x, levels = c('Decline', 'Stable', 'Improvement'))) %>%
  mutate_at(vars('T1', 'T2', 'T3'), ~factor(.x, levels = c('\u2265 4h', '[2h,4h[', '[0h,2h[')))

ggplot(Sim_CPAP_alluvial, aes(axis1 = `T1`, axis2 = `T2`, axis3 = `T3`, y = Freq)) +
  geom_alluvium(aes(fill = `Evolution T1-T3`)) + 
  geom_stratum() + 
  geom_text(stat = 'stratum', aes(label = after_stat(stratum)), size = 5) +
  theme_void() +
  scale_fill_manual(values = c('#C6DBEF', '#6BAED6', '#08519C')) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12))
  
#C: Mosaic plot
#2 time plots

Sim_CPAP_mosaic <- table(Sim_CPAP_cat$T1, Sim_CPAP_cat$T2)

mosaicplot(Sim_CPAP_mosaic, color = c('#08519C','#6BAED6', '#C6DBEF'),
           xlab = 'T1', ylab = 'T2', cex.axis = 1, main = '')


#--------------------Figure 2---------------------------------------------------
#A: Spaghetti and spline plots
#5 time points

Sim_CPAP_spaghetti <- Sim_CPAP %>%
  select('T1':'T5', 'patient_id') %>%
  pivot_longer(cols = c(T1, T2, T3, T4, T5), names_to = 'Time', values_to = 'CPAP adherence') %>%
  mutate_at(vars(Time), ~ as.numeric(as.factor(.x)))

ggplot(Sim_CPAP_spaghetti, aes(x = Time, y = `CPAP adherence`)) +
  geom_line(aes(group = patient_id, colour = patient_id),
            alpha = 0.1) +
  geom_point(aes(colour = patient_id)) +
  geom_smooth(method = 'loess') +
  theme_classic() + 
  theme(legend.position = 'none') +
  ylab('CPAP adherence (h/night)') +
  theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15))

#B: Boxplot
#3 time points: months and 5 patients

Sim_CPAP_boxplot <- Sim_CPAP %>%
  pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'CPAP adherence') %>%
  mutate_at(vars(Time), ~ as.numeric(as.factor(.x))) %>%
  group_by(patient_id) %>%
  mutate(Months = ifelse(Time %in% seq(1, 30), 1, ifelse(Time %in% seq(31, 60), 2, 3))) %>%
  filter(patient_id %in% seq(1, 5)) %>%
  mutate_at(vars(Months, patient_id), as.factor)

ggplot(Sim_CPAP_boxplot, aes(x = Months, y = `CPAP adherence`, color = patient_id)) +
  geom_boxplot() +
  theme_classic() +
  ylab('CPAP adherence (h/night)') +
  labs(colour = 'Patient') +
  scale_color_manual(values = c('#DEEBF7', '#9ECAE1', '#4292C6', '#08519C', '#08306B')) +
  theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),axis.text.y = element_text(size = 15),
        legend.title = element_text(size = 15),legend.text = element_text(size = 15))
  

#C: Survival plots
#All patients and all time points

Sim_CPAP_survival <- Sim_CPAP %>%
  pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'CPAP_adherence') %>%
  mutate_at(vars(patient_id), as.factor) %>%
  mutate_at(vars(Time), ~as.numeric(as.factor(.x))) %>%
  mutate_at(vars(CPAP_adherence), ~ifelse(.x >= 4, 1, 0)) %>%
  group_by(patient_id) %>%
  summarise(Time = ifelse(max(CPAP_adherence) == 1, which.max(CPAP_adherence), max(Time)),
            CPAP_adherence = ifelse(Time == 90, 0, CPAP_adherence)) 

fit1 <- survfit(Surv(Time, CPAP_adherence) ~ 1, data = Sim_CPAP_survival)
p1 <- ggsurvplot(fit1, data = Sim_CPAP_survival, risk.table = T, palette = c('#08519C'))  

#Stratification: ESS score
Sim_ESS_cat <- read.csv("C:/Users/HP2/Desktop/Methodo_stat_donnees/Data/300Patients_90TimePoints/Sim_ESS_cat.csv", sep=";")
Sim_ESS_joint <- Sim_ESS_cat %>%
  select(patient_id, T1) %>%
  rename(Drowsy = T1) %>%
  mutate_at(vars(patient_id), as.factor)

Sim_CPAP_survival <- Sim_CPAP_survival %>%
  mutate_at(vars(patient_id), as.factor) %>%
  left_join(Sim_ESS_joint)

fit2 <- survfit(Surv(Time, CPAP_adherence) ~ Drowsy, data = Sim_CPAP_survival)
p2 <- ggsurvplot(fit2, data = Sim_CPAP_survival, risk.table = T,
                 palette = c('#08519C', '#9ECAE1'))


splots <- list()
splots[[1]] <- p1
splots[[2]] <- p2
arrange_ggsurvplots(splots, print = T, ncol = 2, risk.table.height = 0.2)
