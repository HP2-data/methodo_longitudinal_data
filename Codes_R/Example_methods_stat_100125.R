#03/04/2025
#Longitudinal data, trajectories and time series: how to analyze them?
#An example of sleep data

library(dplyr) #for data manipulation
library(tidyr) #to arrange data
library(ggplot2) #for plots
select <- dplyr::select

source("functions.R") #all functions useful for this code

#--------------------Data set---------------------------------------------------
#In this section, functions are available in the 'functions.R' file
#For the simulation, we assume 90 time points (around 3 months) for 300 patients

#Simulation of a continuous variable describing 3 types of CPAP adherence (N(2.2, 1) 
#the 1st cluster (non-adherent); N(4.8, 0.6) for the 2nd cluster (adherent); N(6.3, 0.7)
#for the 3rd cluster (very adherent),with negative values replaced by 0)
Sim_CPAP <- sim_data(300, 90) %>%
  mutate_all(~ ifelse(round(.x) <= 0, 0, .x))

#Simulation of a categorical variable based on the CPAP adherence variable with
#a threshold at 2h and 4h
Sim_CPAP_cat <- Sim_CPAP %>%
  mutate_at(vars(T1:T90), ~ case_when(.x < 2 ~ '[0h,2h[',
                         .x >=2 & .x < 4 ~ '[2h,4h[',
                         .x >= 4 ~ '\u2265 4h')) %>%
  mutate_all(~ as.factor(.x))

#Simulation of a discrete variable describing ESS scores (values from 0 to 24)
Sim_ESS <- sim_data_discrete(300, 2, 24)

#Simulation of a categorical variable based on the ESS scores variable with
#a threshold at 10 (Yes = Drowsy patient)
Sim_ESS_cat <- Sim_ESS %>%
  mutate_at(vars(T1:T2), ~ ifelse(.x >= 10, 'Yes', 'No')) %>%
  mutate_at(vars(T1:T2, patient_id), as.factor)

#--------------------ANOVA method-----------------------------------------------
library(rstatix)
#ANOVA analysis: how can we study the probability of transition from one cluster to
#another between 2 consecutive points in time?
#Continuous variable: CPAP adherence
#All time points and all patients were included

#Data set/data preparation
ANOVA_df <- Sim_CPAP %>%
  pivot_longer(cols = T1:T90, names_to = 'time', values_to = 'Adherence') %>%
  mutate_at(vars(patient_id, time), as.factor)

#Normal values assumption
ggpubr::ggqqplot(ANOVA_df, 'Adherence', facet.by = 'time')

#ANOVA test
get_anova_table(anova_test(data = ANOVA_df, dv = Adherence, wid = patient_id,
                           within = time))
#p = 0.96 -> we can't say that the CPAP adherence is different at the different 
#time point

#--------------------chi² method-----------------------------------------------
#chi² Mantel-Haenszel: how can we study the probability of transition from one
#cluster to another between 2 consecutive points in time?
#2 nominal variables are conditionally independent in each stratum assuming that
#there is no 3-way interaction
#Selection of 2 time points (ESS score measured at 2 time points only) and
#calculation of contingency table for the chi² test
#All patients were included

#Contingency table for each time point
T1 <- table(Sim_ESS_cat$T1, Sim_CPAP_cat$T1)
T2 <- table(Sim_ESS_cat$T2, Sim_CPAP_cat$T2)

tableTimePoint1 <- t(table(as.matrix(select(Sim_ESS_cat, 'T1')),
                           as.matrix(select(Sim_CPAP_cat, 'T1'))))

tableTimePoint2 <- t(table(as.matrix(select(Sim_ESS_cat, 'T2')),
                           as.matrix(select(Sim_CPAP_cat, 'T2'))))

#Final table with all data information
Chi2_test <- array(c(tableTimePoint1, tableTimePoint2),
                   dim = c(3, 2, 2),
                   dimnames = list(Adherence = c('[0h, 2h[', '[2h, 4h[', '\u2265 4h'),
                                                    ESS_score = c('No', 'Yes'),
                                                    Time = c('T1', 'T2')))

#Mantel-Haenszel test
mantelhaen.test(Chi2_test, correct = F)
#p < 0.01 --> the OR does draw away from 1 
# --> there is a difference between groups

#--------------------K-means method---------------------------------------------
library(kml)
#KML: how can we analyze clusters trajectories to study and predict variations
#over time?
#Numerical variables: continuous CPAP adherence
#All patients and all time points were included

#Data set/data preparation
Kmeans <- Sim_CPAP %>%
  select(patient_id, T1:T90) %>%
  mutate_at(vars(patient_id), as.character)
  

#1st: transform the long format data into clusterLongData with the function cld()
clus_KML <- cld(Kmeans, timeInData = 2:91, maxNA = 1)

#2nd: run the kml function to create clusters
#KLM application
#We use 15 redrawings for each of the clusters and performed fo 2 to 6 clusters
kml(clus_KML, nbRedrawing = 15)
choice(clus_KML)

#Best model with 3 clusters: Calinski-Harabasz score higher for 3 clusters,
#N = 33% for each group, the cluster C decreased its adherence over the 90 time points,
#the cluster B decreased its adherence from the 31st time point to the 60th then
#increased up to the 90th time point, the cluster A increased from the 30th time
#point to the 60th time point then decreased over the last 30 time points

#--------------------LTA method-------------------------------------------------
#source: https://github.com/cran/LMest/blob/master/R/lmest.R
library(LMest)
#LTA: how can we compare the characterisitc of different clusters measured
#regularly over time?
#Discrete or categorical outcome, we chose categorical CPAP adherence
#All patients and all time points were included

#Data set/data preparation
#Data in long format
Sim_CPAP_LTA <- Sim_CPAP_cat %>%
  select(patient_id, T1:T90) %>%
  pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'CPAP_adherence') %>%
  mutate_at(vars(patient_id, Time, CPAP_adherence), as.factor) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-'T') %>%
  mutate_all(as.numeric) %>%
  as.data.frame()

#LTA application
#Tests for 2 to 4 clusters
lmest(data = Sim_CPAP_LTA, index = c('patient_id', 'Time'),
                  start = 1, seed = 123, maxit = 5000, modBasic = 1, k = 2:4)

#Best model for 3 clusters
LTA_final_model <- lmest(data = Sim_CPAP_LTA, index = c('patient_id', 'Time'),
                         start = 1, seed = 123, maxit = 5000, modBasic = 1, k = 3)

summary(LTA_final_model)

#Plot results
plot(LTA_final_model, what = "CondProb")
plot(LTA_final_model, what = "transitions")
plot(LTA_final_model, what = "marginal")

#--------------------GBTM method------------------------------------------------
#source: https://github.com/gitedric/trajeR/tree/master
devtools::install_github("gitedric/trajeR")
library(trajeR)
#GBTM: how can we analyze clusters trajectories to study and predict variations
#over time?
#In columns: patient id, CPAP adherence (numeric), time points (numeric)
#All patients but only 5 time points were included (faster calculation)

#Data set/data preparation
Sim_ESS_GBTM <- Sim_ESS %>%
  rename(ESS_baseline = T1) %>%
  mutate_at(vars(ESS_baseline), as.factor)

Sim_CPAP_GBTM <- Sim_CPAP %>%
  select(patient_id, T1:T5) %>%
  mutate(Time1 = 1, Time2 = 2, Time3 = 3, Time4 = 4, Time5 = 5) %>%
  left_join(select(Sim_ESS_GBTM, patient_id, ESS_baseline))

#GBTM application
#Test with 2 clusters and linear curve
GBTM_test12 <- trajeR(Sim_CPAP_GBTM[,2:6], Sim_CPAP_GBTM[,7:11], ng = 2,
                      degre = rep(1, 2), Model = 'CNORM')

trajeRBIC(GBTM_test12)
trajeRAIC(GBTM_test12)

#Test with 2 clusters and quadratic curve
GBTM_test22 <- trajeR(Sim_CPAP_GBTM[,2:6], Sim_CPAP_GBTM[,7:11], ng = 2,
                      degre = rep(2, 2), Model = 'CNORM')

trajeRBIC(GBTM_test22)
trajeRAIC(GBTM_test22)

#Linear better than quadratic according to likelihood but BIC diff < 10 so
#test with cubic curve

#Test with 2 clusters and cubic curve
GBTM_test32 <- trajeR(Sim_CPAP_GBTM[,2:6], Sim_CPAP_GBTM[,7:11], ng = 2,
                      degre = rep(3, 2), Model = 'CNORM')

trajeRBIC(GBTM_test32)
trajeRAIC(GBTM_test32)

#BIC difference > 10 so final model uses quadratic curve
#Then, test with quadratic curve but different number of clusters
#Test with 3 clusters and quadratic curve
GBTM_test23 <- trajeR(Sim_CPAP_GBTM[,2:6], Sim_CPAP_GBTM[,7:11], ng = 3,
                      degre = rep(2, 3), Model = 'CNORM')

trajeRBIC(GBTM_test23)
trajeRAIC(GBTM_test23)

adequacy(GBTM_test22, Y = Sim_CPAP_GBTM[,2:6], A = Sim_CPAP_GBTM[,7:11]) 
adequacy(GBTM_test23, Y = Sim_CPAP_GBTM[,2:6], A = Sim_CPAP_GBTM[,7:11]) 
#AvePP >= 0.7
#3 clusters better than 2 or 4 clusters for BIC criteria and likelihood too

#Plot trajectories
plotrajeR(GBTM_test23)

#Probability of belonging for each patient
Group_N <- GroupProb(GBTM_test23, Y = Sim_CPAP_GBTM[,2:6], A = Sim_CPAP_GBTM[,7:11])

#Number of patient in each group
Group_N <- Group_N %>%
  as.data.frame() %>%
  mutate(patient_id = seq(1, 300, 1)) %>%
  group_by(patient_id) %>%
  mutate(group = which.max(c(Gr1, Gr2, Gr3)))

questionr::freq(Group_N$group)

#--------------------GMM method-------------------------------------------------
#source: file:///C:/Users/HP2/Downloads/GMM%20in%20R_Dec2022_v2.pdf
library(lcmm)
#GMM: how can we described trajectories of longitudinal data with repeated measurement
#of follow-up?
#Continuous data for CPAP adherence
#All patients but only 5 time points were included (long runtimes with many time points)

#Data set/data preparation
#Data in long format
Sim_CPAP_GMM <- Sim_CPAP %>%
  select(patient_id, T1:T5) %>%
  pivot_longer(cols = c(T1:T5), names_to = 'Time', values_to = 'CPAP_adherence')

#GMM application
#Random effect for the intercept and slope
#test: 1, 2, 3, 4 clusters
set.seed(123)
GMM_test1 <- hlme(CPAP_adherence ~ Time, subject = 'patient_id',
                  random = ~ 1 + Time, ng = 1, data = Sim_CPAP_GMM)
GMM_test2 <- gridsearch(rep = 100, maxiter = 10, minit = GMM_test1,
                        hlme(CPAP_adherence ~ Time, subject = 'patient_id',
                             random = ~ 1 + Time, ng = 2, data = Sim_CPAP_GMM,
                             mixture= ~ Time, nwg = T))
GMM_test3 <- gridsearch(rep = 100, maxiter = 10, minit = GMM_test1,
                        hlme(CPAP_adherence ~ Time, subject = 'patient_id',
                             random = ~ 1 + Time, ng = 3, data = Sim_CPAP_GMM,
                             mixture= ~ Time, nwg = T))
GMM_test4 <- gridsearch(rep = 100, maxiter = 10, minit = GMM_test1,
                        hlme(CPAP_adherence ~ Time, subject = 'patient_id',
                             random = ~ 1 + Time, ng = 4, data = Sim_CPAP_GMM,
                             mixture= ~ Time, nwg = T))

#Parameters for choosing the best model
summarytable(GMM_test1, GMM_test2, GMM_test3, GMM_test4)
#Better BIC = 4 clusters, distribution of patients ok
summary(GMM_test4)

#Posterior probabilities
postprob(GMM_test4)

#Plot fixed effect in the longitudinal model
#Create a table with names of the variables, coefficients, CI_inf, CI_sup, P-values
#and cluster: "GMM_fixed"
ggplot(GMM_fixed, aes(colour = Cluster, y = beta, x = variable)) +
  geom_pointrange(aes(y = beta, ymin = IC_inf, ymax = IC_sup),
                  position = position_dodge(width = 0.25)) +
  ylab('h/night') +
  xlab('Time points') +
  theme_classic() +
  scale_color_manual(values = c('#C2D9ED', '#A4CCE3', '#539DCC', '#08306B')) +
  annotate('text', x = 1.07, y = -0.25, color = '#539DCC', label = '*') +
  annotate('text', x = 1.15, y = 0.47, color = '#08306B', label = '*') +
  annotate('text', x = 2.15, y = 0.52, color = '#08306B', label = '*') +
  annotate('text', x = 3.02, y = -0.64, color = '#A4CCE3', label = '*') +
  annotate('text', x = 3.15, y = 1.13, color = '#08306B', label = '*') +
  annotate('text', x = 4.02, y = -0.77, color = '#A4CCE3', label = '*') +
  annotate('text', x = 4.15, y = 1.13, color = '#08306B', label = '*')

#--------------------Mixed method-----------------------------------------------
library(lme4)
#Mixed: how can we estimate the relationship between the dependent variables and
#the fixed and random effects of the independents variables?
#Continuous or categorical variable / we show the example with continuous outcome:
#Time and ESS score as categorical variable
#All patients and all time points were included

#Data set/data preparation
#Data in long format
Sim_CPAP_mixed <- Sim_CPAP %>%
  pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'CPAP_adherence') %>%
  mutate_at(vars(patient_id, Time), as.factor) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-'T') %>%
  mutate_at(vars(Time), ~ factor(.x, levels = sort(unique(as.numeric(.x))))) %>% #ordered values
  mutate_at(vars(Time), as.numeric) #according to the research question (run or not this line)

#ESS at baseline (T1) = adjustment factor
Sim_ESS_mixed <- Sim_ESS_cat %>%
  select(T1, patient_id) %>%
  mutate_at(vars(patient_id, T1), as.factor)

Sim_CPAP_mixed <- Sim_CPAP_mixed %>%
  left_join(Sim_ESS_mixed, by = 'patient_id') %>%
  rename(ESS_baseline = T1)

#Mixed model application
#Random effect for the intercept, by patient id
Mixed_test <- lmer(CPAP_adherence ~ Time + ESS_baseline + (1 + Time | patient_id),
                   data = Sim_CPAP_mixed)
summary(Mixed_test)

#Plot model & results
sjPlot::plot_model(Mixed_test)
sjPlot::tab_model(Mixed_test)

#-------------------------Survival method------------------------------------------
library(survival) #For Surv() functions; fit survival models
library(survminer) #For ggsurvplot() function; plot the survival models
#survival: how can we determine the survival probability or cumulative risk of a
#population over a defined period of time.
#Continuous CPAP adherence and categorical ESS scores for Cox model
#All patients and all time points were included

#Data set for survival model
#The dead of the patient is the moment when its CPAP adherence is >= 4h
#Time = max between time when CPAP adherence == 1 and T90 (the end of the study)
#So the status is dead (CPAP_adherence = 1) when CPAP adherence >= 4h before the end of the study or 
#censored (CPAP_adherence = 0) when CPAP adherence < 4h over the study
Sim_CPAP_survival <- Sim_CPAP %>%
  pivot_longer(cols = c(T1:T90), names_to = 'Time', values_to = 'CPAP_adherence') %>%
  mutate_at(vars(patient_id), as.factor) %>%
  mutate_at(vars(Time), ~as.numeric(as.factor(.x))) %>%
  mutate_at(vars(CPAP_adherence), ~ifelse(.x >= 4, 1, 0)) %>%
  group_by(patient_id) %>%
  summarise(Time = ifelse(max(CPAP_adherence) == 1, which.max(CPAP_adherence), max(Time)),
            CPAP_adherence = ifelse(Time == 90, 0, CPAP_adherence)) 

#Survival model
fit <- survfit(Surv(Time, CPAP_adherence) ~ 1, data = Sim_CPAP_survival)
#Plot results
ggsurvplot(fit, risk.table = T)

#Data set for survival plot according to drowsy status
#ESS score: Drowsy patient = (ESS >= 10) = 0 /  Non-drowsy patient = (ESS < 10) = 1
#Time = max between time when ESS score == 1 and T90 (the end of the study)
Sim_ESS_joint <- Sim_ESS_cat %>%
  select(patient_id, T1) %>%
  rename(Drowsy = T1) %>%
  mutate_at(vars(patient_id), as.factor) 

Sim_CPAP_survival <- Sim_CPAP_survival %>%
  left_join(Sim_ESS_joint)

#Survival model
fit <- survfit(Surv(Time, CPAP_adherence) ~ Drowsy, data = Sim_CPAP_survival)
#Plot results
ggsurvplot(fit, risk.table = T)

#Comparison test
survdiff(Surv(Time, CPAP_adherence) ~ ESS_score, data = Sim_CPAP_survival) #p>0.05 => no difference

#--------------------ARIMA & CCF method-----------------------------------------
#source for CCF method: https://online.stat.psu.edu/stat510/lesson/8/8.2
library(forecast) #For auto.arima() function
#ARIMA & CCF: how can we analyze time series and evaluate the correlation between
#2 time series varying over time, coinciding or not over time intervals?
#Numerical outcome for time series
#1 time serie per patient: 1 for CPAP adherence and 1 for ESS score
#Then CCF between both time series
#All time points but 1 individual were included

#Data set/data preparation
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

##ARIMA application
ts_CPAP <- ts(Sim_CPAP_ARIMA$CPAP_adherence)
ts_ESS <- ts(Sim_ESS_ARIMA$ESS_score)

#Describe time serie
fUnitRoots::urkpssTest(ts_CPAP, type = c("tau"), lags = c("short"),
                       use.lag = NULL, doplot = TRUE)
tsstationary <- diff(ts_CPAP, differences = 1)
plot(tsstationary)
acf(ts_CPAP) #autocorrelation = correlation between the data and themselves
pacf(ts_CPAP) #partial autocorrelation = at lag k, correlation between all data 
#points that are exactly k steps apart- after accounting for their correlation 
#with the data between those k steps

#Automatic research for the best ARIMA model: auto.arima()
ARIMA_CPAP <- auto.arima(ts_CPAP)
summary(ARIMA_CPAP) #ARIMA(0, 1, 1)

ARIMA_ESS <- auto.arima(ts_ESS)
summary(ARIMA_ESS) #ARIMA(0, 0, 0)

#Test of validation
qqnorm(ARIMA_CPAP$residuals)
qqline(ARIMA_CPAP$residuals)
acf(ARIMA_CPAP$residuals)
Box.test(ARIMA_CPAP$residuals, type = 'Ljung-Box') #p>0.05 --> not significant
#=> no pattern in the residuals => ok

##CCF application
#CCF between ESS and CPAP variables modified (detrend) by the ARIMA function
ESS_res <- ARIMA_ESS$residuals
CPAP_res <- ARIMA_CPAP$residuals

CCF_CPAP_ESS <- ccf(CPAP_res, ESS_res)

#Validation
acf(ts.intersect(CPAP_res, ESS_res))
CCF_CPAP_ESS
plot(CCF_CPAP_ESS)
astsa::lag2.plot(CPAP_res, ESS_res, 15) #Weak correlation between ESS and CPAP  

#Can add it to regressions: ESS according to CPAP adherence parameter (at different lags) 
lag <- stats::lag
final_data <- ts.intersect(ESS_res, CPAPlag14 = lag(CPAP_res, 14))

reg_CCF <- lm(ESS_res ~ CPAPlag14, data = final_data)
summary(reg_CCF)
astsa::acf2(residuals(reg_CCF))

#Can add covariates to the ARIMA method --> ARIMAX

#-------------------------Hidden Markov method------------------------------------------
#source: https://www.geeksforgeeks.org/hidden-markov-model-in-r/
library(depmixS4)
#Hidden Markov: how can we assess changes in individual characteristics when these
#are not directly observable?
#1 known categorical variable: CPAP adherence with 3 states
#The last 200 patients and all time points were included

#Data set/data preparation
Sim_CPAP_HMM <- Sim_CPAP_cat %>%
  filter(patient_id %in% seq(100, 300)) %>%
  mutate_at(vars('T1':'T90'), ~ as.numeric(as.factor(.x))) %>%
  pivot_longer(cols = c('T1':'T90'), names_to = 'Time', values_to = 'CPAP_adherence')
  
#To take into account patient_id: the lengths of individual, i.e. independent, time series.
#If not specified, the responses are assumed to from a single time series.
ntimes <- Sim_CPAP_HMM %>%
  group_by(patient_id) %>%
  summarise(N = length(CPAP_adherence))

Sim_CPAP_HMM <- Sim_CPAP_HMM %>%
  arrange(patient_id) #To coincide with the ntimes table

#Hidden Markov application
#Number of hidden states: 2 -> Adherent and Non adherent; and EM algo
#multinomial() for categorical observations
#variable ~ 1 for independence of observation variable from all covariates 
HMM_test <- depmix(response = CPAP_adherence ~ 1, family = multinomial(), nstates = 2,
                   data = Sim_CPAP_HMM, ntimes = ntimes$N)

set.seed(123)
HMM_final <- fit(HMM_test)

#Test the model
BIC(HMM_final)
AIC(HMM_final)
logLik(HMM_final)

#Initial state probability
init_prob <- matrix(getpars(HMM_final)[1:2], nrow = 1, byrow = T)
colnames(init_prob) <- c('Non adherent', 'Adherent')
init_prob

#Transition matrix
trans_prob <- matrix(getpars(HMM_final)[3:6], nrow = 2, byrow = T)
colnames(trans_prob) <- c('Non adherent', 'Adherent')
rownames(trans_prob) <- c('Non adherent', 'Adherent')
trans_prob

#Predict the hidden states up to 200 time points
pred_states <- posterior(HMM_final, type = 'viterbi')[1:200,]
pred_states

ggplot(pred_states, aes(y = state, x = seq(1, length(state)))) +
  geom_line() +
  ylab('State prediction') +
  xlab('Time points') +
  theme_classic() + 
  scale_y_discrete(breaks = c(1, 2), limits = factor(c(1, 2)))

questionr::freq(pred_states$state)


#-------------------Other more complex methods----------------------------------
#Further examples of more complex methods

##-------------------------RI-CLPM method---------------------------------------
#source: https://rpubs.com/tpartridge/1204398
library(lavaan)
library(semPlot) #For the figure

#Ri-CLPM: how can we ass reciprocal causal effects or mediation effects taking
#into account time-invariant and trait-like stability?
#2 continuous variables: CPAP adherence and ESS score
#All patients but only 5 time points (2 time points too small so add 3 time
#points to the ESS score data set)

#Data set/data preparation
Sim_ESS_RICLPM <- sim_data_discrete(300, 5, 24) %>%
  rename(ESS1 = T1, ESS2 = T2, ESS3 = T3, ESS4 = T4, ESS5 = T5)

data_RICLPM <- Sim_CPAP %>%
  select(patient_id, T1, T15, T30, T55, T80) %>%
  rename(CPAP1 = T1, CPAP2 = T15, CPAP3 = T30, CPAP4 = T55, CPAP5 = T80) %>%
  left_join(Sim_ESS_RICLPM, by = 'patient_id')

#RI-CLPM application
#Function formula
riclpm <- '

# Define intercept factors

i_CPAP =~ 1*CPAP1+1*CPAP2+1*CPAP3+1*CPAP4+1*CPAP5
i_ESS =~ 1*ESS1+1*ESS2+1*ESS3+1*ESS4+1*ESS5

# Define single item latent variables

eta_CPAP1 =~ 1*CPAP1
eta_CPAP2 =~ 1*CPAP2 
eta_CPAP3 =~ 1*CPAP3
eta_CPAP4 =~ 1*CPAP4
eta_CPAP5 =~ 1*CPAP5
eta_ESS1 =~ 1*ESS1
eta_ESS2 =~ 1*ESS2
eta_ESS3 =~ 1*ESS3
eta_ESS4 =~ 1*ESS4
eta_ESS5 =~ 1*ESS5

# Autoregressive effects
eta_CPAP2 ~ a1*eta_CPAP1
eta_CPAP3 ~ a1*eta_CPAP2
eta_CPAP4 ~ a1*eta_CPAP3
eta_CPAP5 ~ a1*eta_CPAP4
eta_ESS2 ~ a2*eta_ESS1
eta_ESS3 ~ a2*eta_ESS2
eta_ESS4 ~ a3*eta_ESS3
eta_ESS5 ~ a3*eta_ESS4

# Crosslagged effects

eta_ESS2 ~ c1*eta_CPAP1
eta_ESS3 ~ c1*eta_CPAP2
eta_ESS4 ~ c1*eta_CPAP3
eta_ESS5 ~ c1*eta_CPAP4
eta_CPAP2 ~ c2*eta_ESS1
eta_CPAP3 ~ c2*eta_ESS2
eta_CPAP4 ~ c2*eta_ESS3
eta_CPAP5 ~ c2*eta_ESS4

# Some further constraints on the variance structure

# 1. Set error variances of the observed variables to zero

CPAP1 ~~ 0*CPAP1
CPAP2 ~~ 0*CPAP2
CPAP3 ~~ 0*CPAP3
CPAP4 ~~ 0*CPAP4
CPAP5 ~~ 0*CPAP5
ESS1 ~~ 0*ESS1
ESS2 ~~ 0*ESS2
ESS3 ~~ 0*ESS3
ESS4 ~~ 0*ESS4
ESS4 ~~ 0*ESS5

# 2. Let lavaan estimate the variance of the latent variables
eta_CPAP1 ~~ varCPAP1*eta_CPAP1
eta_CPAP2 ~~ varCPAP2*eta_CPAP2
eta_CPAP3 ~~ varCPAP3*eta_CPAP3
eta_CPAP4 ~~ varCPAP4*eta_CPAP4
eta_CPAP5 ~~ varCPAP5*eta_CPAP5
eta_ESS1 ~~ varESS1*eta_ESS1
eta_ESS2 ~~ varESS2*eta_ESS2
eta_ESS3 ~~ varESS3*eta_ESS3
eta_ESS4 ~~ varESS4*eta_ESS4
eta_ESS5 ~~ varESS5*eta_ESS5

# 3. We also want estimates of the intercept factor variances and an
#    estimate of their covariance
i_CPAP ~~ variCPAP*i_CPAP
i_ESS ~~ variESS*i_ESS
i_CPAP ~~ covi*i_ESS

# 4. We have to define that the covariance between the intercepts and
#    the latents of the first time point are zero

eta_CPAP1 ~~ 0*i_CPAP
eta_ESS1 ~~ 0*i_CPAP
eta_CPAP1 ~~ 0*i_ESS
eta_ESS1 ~~ 0*i_ESS

# 5. Finally, we estimate the covariance between the latents of x and y
#    of the first time point, the second time-point and so on. Note that
#    for the second to fourth time point the correlation is constrained to
#    the same value

eta_CPAP1 ~~ cov1*eta_ESS1
eta_CPAP2 ~~ e1*eta_ESS2
eta_CPAP3 ~~ e1*eta_ESS3
eta_CPAP4 ~~ e1*eta_ESS4
eta_CPAP5 ~~ e1*eta_ESS5

# The model also contains a mean structure and we have to define some
# constraints for this part of the model. The assumption is that we
# only want estimates of the mean of the intercept factors. All other means
# are defined to be zero:
CPAP1 ~ 0*1
CPAP2 ~ 0*1
CPAP3 ~ 0*1
CPAP4 ~ 0*1
CPAP5 ~ 0*1
ESS1 ~ 0*1
ESS2 ~ 0*1
ESS3 ~ 0*1
ESS4 ~ 0*1
ESS5 ~ 0*1
eta_CPAP1 ~ 1
eta_CPAP2 ~ 1
eta_CPAP3 ~ 1
eta_CPAP4 ~ 1
eta_CPAP5 ~1
eta_ESS1 ~ 1
eta_ESS2 ~ 1
eta_ESS3 ~ 1
eta_ESS4 ~ 1
eta_ESS5 ~ 1
i_CPAP ~ 1
i_ESS ~ 1

## Define correlations
cori := covi / (sqrt(variCPAP) * sqrt(variESS))
cor1 := cov1 / (sqrt(varCPAP1) * sqrt(varESS1))
cort2 := e1 / (sqrt(varCPAP2) * sqrt(varESS2))
cort3 := e1 / (sqrt(varCPAP3) * sqrt(varESS3))
cort4 := e1 / (sqrt(varCPAP4) * sqrt(varESS4))
cort5 := e1 / (sqrt(varCPAP5) * sqrt(varESS5))
'

#Fit the model
riclpm_fit <- sem(riclpm, estimator = "MLR", data = data_RICLPM, mimic = "Mplus",
                  missing = "FIML")

#Results
summary(riclpm_fit, fit.measures = TRUE, standardized = TRUE)

#Figure
semPaths(riclpm_fit, "std", layout = "tree2", edge.label.cex = 0.8, curvePivot = TRUE)


##-------------------------DTW method-------------------------------------------
#source: https://dtw.r-forge.r-project.org/
library(dtw)
#DTW: how can we look at the similarity between 2 time series?
#2 numerical variables for our example but it is possible to use numerical template
#and categorical query. Not possible with categorical template! Continuous CPAP
#adherence and discrete ESS scores
#All time points but only 1 individual were included

#Data set/data preparation
Sim_CPAP_DTW <- Sim_CPAP %>%
  filter(patient_id == 25) %>%
  as.numeric()

Sim_ESS_DTW <- Sim_ESS %>%
  filter(patient_id == 25) %>%
  as.numeric()

#DTW application
#Find the best match with the canonical recursion formula
align <- dtw(Sim_CPAP_DTW, Sim_ESS_DTW, keep = T)

#Alignment curve
plot(align, type = 'threeway', xlab = 'CPAP adherence over time',
     ylab = 'ESS score over time / reference index')

#Align and plot for the 90 measuring points
plot(dtw(Sim_CPAP_DTW[1:90], Sim_ESS_DTW[1:90], keep = T,
              step = rabinerJuangStepPattern(6,"c")),
     type = "twoway", offset = -2, ylab = 'CPAP adherence', xlab = 'Measuring points')


##--------------------LCA method------------------------------------------------
library(poLCA)
#LCA: how can we identify unmeasured clusters sharing common characteristics?
#Categorical variables for CPAP adherence
#All patients and all time points were included

#Data set/data preparation
LCA_data <- Sim_CPAP_cat %>%
  mutate_all(as.factor)

LCA_function <- as.formula(paste0('cbind(', paste0('T', seq(1,90), collapse = ','), ') ~ 1'))

#LCA application
#Choice of the number of clusters using BIC and AIC criterion
set.seed(123)
LCA_test <- poLCA(LCA_function, LCA_data, nclass = 2) #BIC = 35930.7 / AIC = 24593.6

set.seed(123)
poLCA(LCA_function, LCA_data, nclass = 4) #BIC = 24516.6 / AIC = 21838.8

set.seed(123)
poLCA(LCA_function, LCA_data, nclass = 5) #BIC = 25320.5 / AIC = 21972.3

set.seed(123)
LCA_test <- poLCA(LCA_function, LCA_data, nclass = 3) #BIC = 23544.7 / AIC = 21537.2
#number chosen = 3, bigger BIC and AIC for other number of clusters

#Class membership probabilities
LCA_test$P

#Item-response probabilities
LCA_test$probs

#Example visualization of item-response probabilities
#source: https://www.geeksforgeeks.org/latent-class-analysis-in-r/
plot_lca(LCA_test)


##-------------------------Joint method-----------------------------------------
library(nlme) #For lme() function; fit the mixed model
library(survival) #For Surv() and coxph() functions; fit survival models
library(survminer) #For ggsurvplot() function; plot the survival models
library(JMbayes2)   #For jm() function; fit joint model
library(sjPlot) #For plot_model() function; plot lme regression results
#joint: how can we describe the joint behavior of the evolution of a quantitative
#longitudinal marker and the time of occurrence of an event considering their joint
#density?
#Continuous CPAP adherence for mixed model and categorical ESS scores for Cox model
#All patients but only 7 time points were included

#Data set/data preparation
#1 covariate: sex of the patient
Sex <- sample(x = c('M', 'F'), replace = T, size = 300) %>%
  as.data.frame() %>%
  rename('Sex' = '.') %>%
  mutate(patient_id = seq(1,300,1)) %>%
  mutate_at(vars(patient_id, Sex), as.factor)

#Data set for mixed model
Sim_CPAP_joint <- Sim_CPAP %>%
  select(patient_id, T1:T5) %>%
  pivot_longer(cols = c(T1:T5), names_to = 'Time', values_to = 'CPAP_adherence') %>%
  mutate_at(vars(patient_id), as.factor) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-c('T')) %>%
  mutate_at(vars(Time), as.numeric) %>%
  group_by(patient_id) %>%
  left_join(Sex)

#Data set for Cox and Joint model
#ESS score: Drowsy patient = (ESS >= 10) = 0 /  Non-drowsy patient = (ESS < 10) = 1
#Time = max between time when ESS score == 1 and T7 (the end of the study)
Sim_ESS_joint <- sim_data_discrete(300, 5, 24) %>%
  select(patient_id, T1:T5) %>%
  pivot_longer(cols = c(T1:T5), names_to = 'Time', values_to = 'ESS_score') %>%
  mutate_at(vars(patient_id), as.factor) %>%
  group_by(patient_id) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-c('T')) %>%
  mutate_at(vars(Time), as.numeric) %>%
  left_join(Sim_CPAP_joint, by = c('patient_id', 'Time')) %>%
  mutate_at(vars(ESS_score), ~ ifelse(.x >= 10, 0, 1)) %>%
  group_by(patient_id) %>%
  mutate_at(vars(Time), ~ ifelse(max(ESS_score) == 1, which.max(ESS_score), max(Time))) %>%
  mutate_at(vars(ESS_score), ~ ifelse(Time == 5, 0, 1)) %>%
  distinct(patient_id, .keep_all = T) %>%
  left_join(Sex)

#Mixed model application
#CPAP adherence with random effect on intercept and slope, by patient 
control <- lmeControl(maxIter = 100, msMaxIter = 100, opt = 'optim')
reg_joint <- lme(CPAP_adherence ~  Time + Sex, data = Sim_CPAP_joint,
                 random = ~ 1 + Time| patient_id, control = control, method = 'ML')

summary(reg_joint)

#Validation
qqnorm(resid(reg_joint))
qqline(resid(reg_joint))

#Plot results
plot_model(reg_joint, vline.color = 'black', show.values = T,value.size = 3,
           value.offset = .3, title = 'CPAP adherence (h/night)', colors = 'black') +
  theme_classic()

#Cox model application
#Cox survival model: ESS score (0 == alive (drowsy patient) / 1 == dead (non-drowsy patient))
#Clustered by patient
Cox_joint <- coxph(Surv(Time, ESS_score) ~ Sex + cluster(patient_id),
                   data = Sim_ESS_joint, x = T, model = T)

summary(Cox_joint)

#Plot results
ggsurvplot(survfit(Cox_joint), data = Sim_ESS_joint, risk.table = T)
ggsurvplot(survfit(Surv(Time, ESS_score) ~ Sex + cluster(patient_id), data = Sim_ESS_joint),risk.table = T)

#Validation
cox.zph(Cox_joint) #hypothese de prop des risques relatifs (doit être > 0.05)

#Joint model application
set.seed(123)
fit_jm <- jm(Cox_joint, reg_joint, "Time", id_var = 'patient_id')

summary(fit_jm)

#Verification
ggtraceplot(fit_jm, 'alphas')
ggdensityplot(fit_jm, 'alphas')
ggdensityplot(fit_jm, 'betas')
ggtraceplot(fit_jm, 'betas')
#For CPAP adherence / alphas --> Not good 
#For the others / betas --> Not good 
#Because no convergence, superposition of the 3 chains

#Plot results for one patient (example patient 59)
Sim_ESS_joint <- Sim_ESS_joint %>%
  rename("Time_censored" = 'Time')

New_df <- Sim_CPAP_joint %>%
  left_join(select(Sim_ESS_joint, patient_id, Time_censored, ESS_score)) %>%
  filter(Time <= Time_censored) %>%
  filter(patient_id == 59)

#Mixed model: CPAP adherence
predLong1 <- predict(fit_jm, newdata = New_df, return_newdata = T)
plot(predLong1)
#Surv model: ESS score
predSurv <- predict(fit_jm, newdata = New_df, process ="event", return_newdata = T) 
plot(predSurv)
#Final plot with both models
plot(predLong1, predSurv)
