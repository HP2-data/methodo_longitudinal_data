#10/10/2025
#Methodology for statistical method using longitudinal data
#Example with sleep data

library(dplyr) #for data manipulation
library(tidyr) #to arrange data
library(ggplot2) #for plots

source("functions.R") #all functions useful for this code

#--------------------Data set---------------------------------------------------
#For the simulation, we assume 1000 time points for 50 patients, simulating
#CPAP adherence
Sim_CPAP <- sim_data(50, 1000) %>%
  mutate_all(~ ifelse(.x < 0, 0, .x))

#Simulation of a categorical variable based on the CPAP adherence variable
Sim_CPAP_cat <- Sim_CPAP %>%
  mutate_at(vars(T1:T1000), ~ case_when(.x < 2 ~ '[0h,2h[',
                         .x >=2 & .x < 4 ~ '[2h,4h[',
                         .x >= 4 ~ '\u2265 4h')) %>%
  mutate_all(~ as.factor(.x))

#Simulation of a categorical variable describing ESS scores (values from 0 to 24)
Sim_ESS <- sim_data_discrete(50, 1000, 24)


#--------------------ANOVA method-----------------------------------------------
library(rstatix)
#ANOVA analysis: does a continuous variable differ over time?
#Continuous variable: CPAP adherence
#All time points and all patients were included

#Data set/data preparation
ANOVA_df <- Sim_CPAP %>%
  pivot_longer(cols = T1:T1000, names_to = 'time', values_to = 'Adherence') %>%
  mutate_at(vars(patient_id, time), as.factor)

#ANOVA test
anova_test(data = ANOVA_df, dv = Adherence, wid = patient_id, within = time)
#p = 0.09 -> we can't say that the CPAP adherence is different at the different 
#time point

#--------------------chi² method-----------------------------------------------
#chi² Mantel-Haenszel: does a categorical variable differ over time?
#2 nominal variables are conditionally independent in each stratum assuming that
#there is no 3-way interaction
#Selection of 4 time points and calculation of contingency table for the chi² test
#All patients were included

#Data set/data preparation
ESS_cat <- Sim_ESS %>%
  mutate_at(vars(T1:T1000), ~ ifelse(.x >= 10, 'Yes', 'No')) %>%
  mutate_at(vars(T1:T1000, patient_id), as.factor) 

#Contingency table for each time point
T1 <- table(ESS_cat$T1, Sim_CPAP_cat$T1)
T2 <- table(ESS_cat$T2, Sim_CPAP_cat$T2)
T3 <- table(ESS_cat$T3, Sim_CPAP_cat$T3)
T4 <- table(ESS_cat$T4, Sim_CPAP_cat$T4)

#Final table with all data information
Chi2_test <- array(c(2, 11, 10, 1, 9, 17, 6, 6, 7, 1, 11, 19, 2, 10, 10, 4, 8, 16, 1, 10, 11, 3, 12, 13),
                   dim = c(3, 2, 4), dimnames = list(Adherence = c('[0h, 2h[', '[2h, 4h[', '\u2265 4h'),
                                                    ESS_score = c('No', 'Yes'),
                                                    Time = c('T1', 'T2', 'T3', 'T4')))

#Mantel-Haenszel test
mantelhaen.test(Chi2_test, correct = F)
#p = 0.16 --> the OR doesn't draw away from 1 
# --> the difference between groups is not really evident

#--------------------LCA method-------------------------------------------------
library(poLCA)
select <- dplyr::select
#LCA: how can we identify unmeasured clusters sharing common characteristics?
#Categorical variables for CPAP adherence
#All patients but only 5 time points were included

#Data set/data preparation
LCA_data <- Sim_CPAP_cat %>%
  select(c(T1:T5, patient_id))

LCA_function <- cbind(T1, T2, T3, T4, T5) ~ 1

#LCA application
#Choice of the number of clusters using BIC and AIC criterion
poLCA(LCA_function, LCA_data, nclass = 3) #BIC = 558.8 / AIC = 497.6
poLCA(LCA_function, LCA_data, nclass = 5) #BIC = 620.5 / AIC = 517.3
LCA_test <- poLCA(LCA_function, LCA_data, nclass = 2) #BIC = 533.5 / AIC = 493.4
#number chosen = 2

#Class membership probabilities
LCA_test$P

#Item-response probabilities
LCA_test$probs

#Example visualization of item-response probabilities
#source: https://www.geeksforgeeks.org/latent-class-analysis-in-r/
plot_lca(LCA_test)

#--------------------K-means method---------------------------------------------
library(kml)
#KML: how can we analyze clusters trajectories to study and predict variations
#over time?
#Numerical variables: continuous CPAP adherence
#All patients but only 5 time points were included

#Data set/data preparation
Kmeans <- Sim_CPAP %>%
  select(patient_id, T1:T5) %>%
  mutate_at(vars(patient_id), as.character)
  

#1st: transform the long format data into clusterLongData with the function cld()
clus_KML <- cld(Kmeans, timeInData = 2:6, maxNA = 1)

#2nd: run the kml function to create clusters
#KLM application
#We use 15 redrawings for each of the clusters and performed fo 2 to 6 clusters
kml(clus_KML, nbRedrawing = 15)
choice(clus_KML)
#Best model with 2 clusters: Calinski-Harabatz score higher for 2 clusters,
#N_A = 66% and N_B = 34%, A fairly stable around 4h and B increases after the 3rd
#time point from around 4h to around 5h30

#--------------------GMM method-------------------------------------------------
#source: file:///C:/Users/HP2/Downloads/GMM%20in%20R_Dec2022_v2.pdf
library(lcmm)
#GMM: how can we describe trajectories of longitudinal data with repeated measurement
#of follow-up?
#Continuous data for CPAP adherence
#All patients but only 5 time points were included

#Data set/data preparation
#Data in long format
Sim_CPAP_GMM <- Sim_CPAP %>%
  dplyr::select(patient_id, T1:T5) %>%
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
summarytable(GMM_test1, GMM_test2, GMM_test3, GMM_test4)

#Better BIC = 2 clusters, distribution of patients ok
summary(GMM_test2)
#Posterior probabilities
postprob(GMM_test2)

#Plot fixed effect in the longitudinal model
#Create a table with names of the variables, coefficients, CI_inf, CI_sup, P-values
#and cluster: "GMM_fixed"
ggplot(GMM_fixed, aes(colour = Cluster, y = beta, x = variable)) +
  geom_pointrange(aes(y = beta, ymin = IC_inf, ymax = IC_sup, shape = `P-value`),
                  position = position_dodge(width = 0.25)) +
  ylab('') +
  xlab('') +
  ggtitle('Fixed effect in the longitudinal model') +
  theme_classic()

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
  pivot_longer(cols = c(T1:T1000), names_to = 'Time', values_to = 'CPAP_adherence') %>%
  mutate_at(vars(patient_id, Time), as.factor) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-'T') %>%
  mutate_at(vars(Time), ~ factor(.x, levels = sort(unique(as.numeric(.x))))) #ordered values

Sim_ESS_mixed <- Sim_ESS %>%
  select(T1, patient_id) %>%
  mutate_at(vars(T1), ~ ifelse(.x >= 10, 'Yes', 'No')) %>%
  mutate_at(vars(patient_id, T1), as.factor)

Sim_CPAP_mixed <- Sim_CPAP_mixed %>%
  left_join(Sim_ESS_mixed, by = 'patient_id') %>%
  rename(ESS_baseline = T1)

#Mixed application
#Random effect for the intercept and the slope, by patient id
Mixed_test <- lmer(CPAP_adherence ~ Time + ESS_baseline + (1 | patient_id), data = Sim_CPAP_mixed)
summary(Mixed_test)

#Plot model & results
sjPlot::plot_model(Mixed_test)
sjPlot::tab_model(Mixed_test)

#--------------------GBTM method------------------------------------------------
#source: https://github.com/gitedric/trajeR/tree/master
devtools::install_github("gitedric/trajeR")
library(trajeR)
#GBTM: how can we analyze clusters trajectories to study and predict variations
#over time
#In columns: patient id, CPAP adherence (numeric), time points (numeric)
#All patients but only 5 time points were included

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

#Quadratic better than linear according to likelihood but BIC diff < 10 so test
#with cubic curve

#Test with 2 clusters and cubic curve
GBTM_test32 <- trajeR(Sim_CPAP_GBTM[,2:6], Sim_CPAP_GBTM[,7:11], ng = 2,
                      degre = rep(3, 2), Model = 'CNORM')

trajeRBIC(GBTM_test32)
trajeRAIC(GBTM_test32)
##Linear, quadratic better than cubic curve (BIC and likelihood)
#Consequently, quadratic curve for the final model

#Test with 3 clusters and quadratic curve
GBTM_test23 <- trajeR(Sim_CPAP_GBTM[,2:6], Sim_CPAP_GBTM[,7:11], ng = 3,
                      degre = rep(2, 3), Model = 'CNORM')

trajeRBIC(GBTM_test23)
trajeRAIC(GBTM_test23)

adequacy(GBTM_test22, Y = Sim_CPAP_GBTM[,2:6], A = Sim_CPAP_GBTM[,7:11]) 
adequacy(GBTM_test23, Y = Sim_CPAP_GBTM[,2:6], A = Sim_CPAP_GBTM[,7:11]) 
#AvePP >= 0.7 for both models but Prop = 0.1 for the 1st cluster of the 3-clusters model
#Prop = 0.1 for the 2nd cluster and Prop = 0.06 for the 3rd cluster of the 4-clusters model
#2 clusters better then 3 or 4 clusters for BIC criteria too

#Plot trajectories
plotrajeR(GBTM_test22)
GroupProb(GBTM_test22, Y = Sim_CPAP_GBTM[,2:6], A = Sim_CPAP_GBTM[,7:11])

#--------------------LTA method-------------------------------------------------
#source: https://github.com/cran/LMest/blob/master/R/lmest.R
library(LMest)
#LTA: how can we study the probability of transition from 1 cluster to another
#between 2 consecutive points in time?
#Discrete or categorical outcome, we chose categorical CPAP adherence
#All patients but only 500 time points were included

#Data set/data preparation
#Data in long format
Sim_CPAP_LTA <- Sim_CPAP_cat %>%
  select(patient_id, T1:T500) %>%
  pivot_longer(cols = c(T1:T500), names_to = 'Time', values_to = 'CPAP_adherence') %>%
  mutate_at(vars(patient_id, Time, CPAP_adherence), as.factor) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-'T') %>%
  mutate_all(as.numeric) %>%
  as.data.frame()

#LTA application
LTA_test <- lmest(data = Sim_CPAP_LTA, index = c('patient_id', 'Time'),
                         start = 1, seed = 123, maxit = 5000, modBasic = 1, k = 2:4)

LTA_test
#Best model for 2 clusters

LTA_final_model <- lmest(data = Sim_CPAP_LTA, index = c('patient_id', 'Time'),
                         start = 1, seed = 123, maxit = 5000, modBasic = 1, k = 2)

summary(LTA_final_model)

#Plot results
plot(LTA_final_model, what = "CondProb")
plot(LTA_final_model, what = "transitions")
plot(LTA_final_model, what = "marginal")

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
  select(patient_id, T1:T1000) %>%
  filter(patient_id == 10) %>%
  pivot_longer(cols = c(T1:T1000), names_to = 'Time', values_to = 'CPAP_adherence') %>%
  mutate_at(vars(patient_id), as.factor) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-c('T', 'patient_id'))

Sim_ESS_ARIMA <- Sim_ESS %>%
  select(patient_id, T1:T1000) %>%
  filter(patient_id == 10) %>%
  pivot_longer(cols = c(T1:T1000), names_to = 'Time', values_to = 'ESS_score') %>%
  mutate_at(vars(patient_id), as.factor) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-c('T', 'patient_id'))

##ARIMA application
ts_CPAP <- ts(Sim_CPAP_ARIMA$CPAP_adherence, frequency = 7)
ts_ESS <- ts(Sim_ESS_ARIMA$ESS_score, frequency = 7)

#Describe time serie
fUnitRoots::urkpssTest(ts_CPAP, type = c("tau"), lags = c("short"),
                       use.lag = NULL, doplot = TRUE)
tsstationary <- diff(ts_CPAP, differences=1)
plot(tsstationary)
acf(ts_CPAP) #autocorrelation = correlation between the data and themselves
pacf(ts_CPAP) #partial autocorrelation = at lag k, correlation between all data 
#points that are exactly k steps apart- after accounting for their correlation 
#with the data between those k steps

#Automatic research for the best ARIMA model: auto.arima()
ARIMA_CPAP <- auto.arima(ts_CPAP)
summary(ARIMA_CPAP) #ARIMA(0, 0, 0)

ARIMA_ESS <- auto.arima(ts_ESS)
summary(ARIMA_ESS) #ARIMA(1, 1, 0)(1, 0, 0)

#Test of validation
qqnorm(ARIMA_CPAP$residuals)
qqline(ARIMA_CPAP$residuals)
acf(ARIMA_CPAP$residuals)
Box.test(ARIMA_CPAP$residuals, type = 'Ljung-Box') #p>0.05 --> non significance
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
astsa::lag2.plot(CPAP_res, ESS_res, 10) #No correlation between ESS and CPAP with or without lag

#Can add it to regressions: ESS according to CPAP adherence parameter (at different lags) 
lag <- stats::lag
final_data <- ts.intersect(ESS_res, CPAPlag11 = lag(CPAP_res, -11),
                           CPAPlag1 = lag(CPAP_res, 1), 
                           CPAPlag14 = lag(CPAP_res, 14))

reg_CCF <- lm(ESS_res ~ CPAPlag11 + CPAPlag1 + CPAPlag14, data = final_data)
summary(reg_CCF)
astsa::acf2(residuals(reg_CCF))

#Can add covariates to the ARIMA method --> ARIMAX

#-------------------------DTW method--------------------------------------------
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

#Align and plot for the 1st 100 time points (to be more readable) but we can plot
#all time points
plot(dtw(Sim_CPAP_DTW[1:100], Sim_ESS_DTW[1:100], keep = T,
              step = rabinerJuangStepPattern(6,"c")),
     type = "twoway", offset = -2, main = 'Example for the first 100 time points',
     ylab = 'CPAP adherence', xlab = 'Time points')

#-------------------------Joint method------------------------------------------
library(nlme) #For lme() function; fit the mixed model
library(survival) #For Surv() and coxph() functions; fit survival models
library(survminer) #For ggsurvplot() function; plot the survival models
library(JMbayes2)   #For jm() function; fit joint model
library(sjPlot) #For plot_model() function; plot lme regression results
select <- dplyr::select
#joint:how can we describe the joint behavior of the evolution of a quantitative
#longitudinal marker and the time of occurrence of an event considering their joint
#density?
#Continuous CPAP adherence for mixed model and categorical ESS scores for Cox model
#All patients but only 7 time points were inlcuded

#Data set/data preparation
#1 covariate: sex of the patient
Sex <- sample(x = c('M', 'F'), replace = T, size = 50) %>%
  as.data.frame() %>%
  rename('Sex' = '.') %>%
  mutate(patient_id = seq(1,50,1)) %>%
  mutate_at(vars(patient_id, Sex), as.factor)

Sim_CPAP_joint <- Sim_CPAP %>%
  select(patient_id, T1:T7) %>%
  pivot_longer(cols = c(T1:T7), names_to = 'Time', values_to = 'CPAP_adherence') %>%
  mutate_at(vars(patient_id), as.factor) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-c('T')) %>%
  mutate_at(vars(Time), as.numeric) %>%
  group_by(patient_id) %>%
  left_join(Sex)
  

Sim_ESS_joint <- Sim_ESS %>%
  select(patient_id, T1:T7) %>%
  pivot_longer(cols = c(T1:T7), names_to = 'Time', values_to = 'ESS_score') %>%
  mutate_at(vars(patient_id), as.factor) %>%
  group_by(patient_id) %>%
  separate(Time, sep = 1, into = c('T', 'Time')) %>%
  select(-c('T')) %>%
  mutate_at(vars(Time), as.numeric) %>%
  left_join(Sim_CPAP_joint, by = c('patient_id', 'Time')) %>%
  mutate_at(vars(ESS_score), ~ ifelse(.x >= 10, 0, 1)) %>%
  group_by(patient_id) %>%
  mutate_at(vars(Time), ~ ifelse(max(ESS_score) == 1, which.max(ESS_score), max(Time))) %>%
  mutate_at(vars(ESS_score), ~ ifelse(Time == 7, 0, 1)) %>%
  distinct(patient_id, .keep_all = T) %>%
  left_join(Sex)

#Mixed model application
#1 mixed regression: CPAP adherence with random effect on patient 
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
#1 Cox survival model: ESS score (0 == ESS >=10 -> alive / 1 == ESS < 10 -> dead)
Cox_joint <- coxph(Surv(Time, ESS_score) ~ Sex + cluster(patient_id),
                             data = Sim_ESS_joint, x = T, model = T)

summary(Cox_joint)

#Plot results
ggsurvplot(survfit(Cox_joint), data = Sim_ESS_joint, risk.table = T)

#Validation
cox.zph(Cox_joint) #hypothese de prop des risques relatifs (doit être > 0.05)

#Joint model application
set.seed(123)
fit_jm <- jm(Cox_joint, reg_joint, "Time", id_var = 'patient_id')

summary(fit_jm)

#Verifications
ggtraceplot(fit_jm, 'alphas')
ggdensityplot(fit_jm, 'alphas')
ggdensityplot(fit_jm, 'betas')
ggtraceplot(fit_jm, 'betas')
#For CPAP adherence / alphas --> Not good but for the others / betas are pretty ok

#Plot results
Sim_ESS_joint <- Sim_ESS_joint %>%
  rename("Time_censored" = 'Time')

Ndf <- Sim_CPAP_joint %>%
  left_join(select(Sim_ESS_joint, patient_id, Time_censored, ESS_score)) %>%
  filter(Time <= Time_censored) %>%
  filter(patient_id == 49)

predLong1 <- predict(fit_jm, newdata = Ndf, return_newdata = T)
plot(predLong1)
predSurv <- predict(fit_jm, newdata = Ndf, process ="event", return_newdata = T) 
plot(predSurv)
plot(predLong1, predSurv)

#-------------------------Hidden Markov method------------------------------------------
#source: https://www.geeksforgeeks.org/hidden-markov-model-in-r/
library(depmixS4)
#Hidden Markov: how can we assess changes in individual characteristics when these
#are not directly observable?
#1 known categorical variable: CPAP adherence with 3 states
#Only 1 patient but all time points were inlcuded

#Data set/data preparation
Sim_CPAP_HMM <- Sim_CPAP_cat %>%
  filter(patient_id == 25) %>%
  as.numeric()

#Hidden Markov application
#Number of hidden states: 2 for example (Adherent and Non adherent) and EM algo
#multinomial() for categorical observations
#variable ~ 1 for independence of observation variable from all covariates 
HMM_test <- depmix(Sim_CPAP_HMM ~ 1, family = multinomial(), nstates = 2,
                             data = data.frame(Sim_CPAP_HMM))

set.seed(123)
HMM_final <- fit(HMM_test)

#Test the model
BIC(HMM_final)
AIC(HMM_final)
logLik(HMM_final)

#Transition matrix
trans_prob <- matrix(getpars(HMM_final)[1:4], nrow = 2, byrow = T)
colnames(trans_prob) <- c('Non adherent', 'Adherent')
rownames(trans_prob) <- c('Non adherent', 'Adherent')
trans_prob

#Emission probabilities
emission_prob <- matrix(getpars(HMM_final)[5:10], nrow = 2, byrow = T)
colnames(emission_prob) <- c('[0h;2h[', '[2h;4h[', '\u2265 4h')
rownames(emission_prob) <- c('Non adherent', 'Adherent')
emission_prob

#Predict the hidden states
pred_states <- posterior(HMM_final, type = 'viterbi')
pred_states

ggplot(pred_states, aes(y = state, x = seq(1, length(state)))) +
  geom_line() +
  ylab('State prediction') +
  xlab('Time points') +
  theme_classic() + 
  scale_y_discrete(breaks = c(1, 2), limits = factor(c(1, 2)))

questionr::freq(pred_states$state)
