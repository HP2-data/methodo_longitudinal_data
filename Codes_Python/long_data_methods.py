# Longitudinal data, trajectories and time series - An example of sleep data

# import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.formula.api import smf
from statsmodels.stats.anova import AnovaRM 
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.stattools import adfuller, ccf
from statsmodels.tsa.tsatools import lagmat
from statsmodels.tools.tools import add_constant
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.graphics.gofplots import qqplot
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.regression.linear_model import OLS
from pmdarima import auto_arima
from tslearn.clustering import TimeSeriesKMeans
from tslearn.utils import to_time_series_dataset
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import calinski_harabasz_score
from hmmlearn.hmm import CategoricalHMM
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# import functions from data_simulator file
from data_simulator import *
    
# ---------------------- Data Simulation ----------------------
# Simulate data using the data simulation function from data_simulator.py
# Simulate continuous CPAP data
Sim_CPAP = sim_data(300, 90)

# Replace negative or zero values with 0 (after rounding)
Sim_CPAP.iloc[:, :-1] = Sim_CPAP.iloc[:, :-1].applymap(lambda x: 0 if round(x) <= 0 else x)

# Convert continuous CPAP to categorical based on thresholds
def categorize_adherence(x):
    if x < 2:
        return '[0h,2h['
    elif 2 <= x < 4:
        return '[2h,4h['
    else:
        return '≥ 4h'

Sim_CPAP_cat = Sim_CPAP.copy()
Sim_CPAP_cat.iloc[:, :-1] = Sim_CPAP_cat.iloc[:, :-1].applymap(categorize_adherence)

# Convert to categorical type
Sim_CPAP_cat.iloc[:, :-1] = Sim_CPAP_cat.iloc[:, :-1].astype('category')

# Simulate ESS scores (discrete)
Sim_ESS = sim_data_discrete(300, 2, 24)

# Convert ESS scores to binary category
Sim_ESS_cat = Sim_ESS.copy()

for col in ['T1', 'T2']:
    Sim_ESS_cat[col] = Sim_ESS_cat[col].apply(lambda x: 'Yes' if x >= 10 else 'No')
    Sim_ESS_cat[col] = Sim_ESS_cat[col].astype('category')

# Optionally: make patient_id a categorical variable too
Sim_ESS_cat['patient_id'] = Sim_ESS_cat['patient_id'].astype('category')


# ---------------------- ANOVA method -----------------------
# Convert Sim_CPAP dataframe from wide to long-format
df_long = Sim_CPAP.melt(id_vars=['patient_id'], var_name='time', value_name='value')
anova_results = AnovaRM(data=df_long, depvar='value',
                  subject='patient_id', within=['time']).fit()
p_value = anova_results.anova_table['Pr > F'][0]
# print(f"ANOVA method, p-value: {p_value}")
# p = 0.96 -> we can't say that the CPAP adherence is significally different across time points

# ---------------------- Mantel-Haenszel method ----------------------

# chi² Mantel-Haenszel: how can we study the probability of transition from one cluster to another 
# between 2 consecutive points in time?
# 2 nominal variables are conditionally independent in each stratum assuming that there is no 3-way interaction
# Selection of 2 time points (ESS score measured at 2 time points only) and calculation of contingency table for the chi² test
# All patients were included

# The package statsmodels used here to perform the MH test supports only contingency tables of 2x2 format
# For this reason, the cpap adherence variable was mapped into 2 categories: High and Low

# Map values in  'T1', 'T2', 'T3', 'T4' columns of ess_cat and cpap_cat datasets
# Define the adherence map function for cpap_cat
def adherence_map(value):
    if value in ['[0h,2h[', '[2h,4h[']:  # Values indicating low cpap adherence
        return 'Low'
    else:  # Values indicating high cpap adherence
        return 'High'

# Map ESS categories (No = 0, Yes = 1) 
ess_map = {'No': 0, 'Yes': 1}

# Apply mapping
Sim_ESS_cat[['T1', 'T2', 'T3', 'T4']] = Sim_ESS_cat[['T1', 'T2', 'T3', 'T4']].apply(lambda col: col.map(ess_map))
Sim_CPAP_cat[['T1', 'T2', 'T3', 'T4']] = Sim_CPAP_cat[['T1', 'T2', 'T3', 'T4']].apply(lambda col: col.map(adherence_map))

# Create contingency tables for each time point (T1, T2, T3, T4)
tables = []
for time in ['T1', 'T2', 'T3', 'T4']:
    table = pd.crosstab(Sim_ESS_cat[time], Sim_CPAP_cat[time])
    # print(f"\nContingency table for {time} (Shape: {table.shape}):")
    # print(table)
    tables.append(table.values)

# Convert to 2x2xk format (transpose ensures correct shape)
contingency_tables = np.stack(tables, axis=2)
# print(contingency_tables.shape)  # Should be (2, 2, 4)

# Apply Mantel-Haenszel test
st = StratifiedTable(contingency_tables.astype(np.float64))

# Get the summary of test results
results = st.summary()

# Get specific results, such as p-value, pooled odds, confidence intervals
p_value = st.test_null_odds().pvalue
pooled_odds = st.logodds_pooled
conf_int = st.logodds_pooled_confint

# Print the results
# print("\nMantel-Haenszel Test Summary:")
# print(st.summary())
# print(f"p-value: {p_value}")
# print("\nPooled Log-Odds:", st.logodds_pooled)
# print("Pooled Log-Odds 95% CI:", st.logodds_pooled_confint())

# ---------------------- K-means ----------------------
# How can we analyze clusters trajectories to study and predict variations over time?
# Numerical variables: continuous CPAP adherence
# All patients and all time points were included

# Copy dataset to avoid modifying the original
df = Sim_CPAP.copy()

# Select time points columns (T1-T90)
time_points = [f'T{i}' for i in range(1, 91)]

# Standardize time series data
scaler = StandardScaler()
df_scaled = scaler.fit_transform(df[time_points])

# Convert to time-series format required for TimeSeriesKMeans
ts_data = to_time_series_dataset(df[time_points].values)

# Try clustering for k = 2 to 6
n_clusters_range = range(2, 7)
calinski_scores = []
models = {}

# Fit TimeSeriesKMeans for different cluster numbers
for n_clusters in n_clusters_range :
    model = TimeSeriesKMeans(n_clusters=n_clusters, metric="euclidean", random_state=42, n_init=15)
    labels = model.fit_predict(ts_data)
    score = calinski_harabasz_score(df[time_points], labels)
    calinski_scores.append(score)
    models[n_clusters] = model

# Identify the optimal number of clusters
optimal_clusters = n_clusters_range[np.argmax(calinski_scores)]
# print(f"Optimal number of clusters: {optimal_clusters}")

# Fit final TimeSeriesKMeans model with optimal clusters
best_model = models[optimal_clusters]
df['cluster'] = best_model.predict(ts_data)

# Reshape to long format for plotting
df_melted = df.melt(id_vars='cluster', value_vars=time_points, var_name='Time Point', value_name='Adherence')

# Plot adherence patterns per cluster
plt.figure(figsize=(12, 6))
sns.lineplot(data=df_melted, x='Time Point', y='Adherence',
             hue='cluster', estimator='mean', errorbar=None)

# Show only every 5th tick
ticks = list(range(0, len(time_points), 5))
labels = [time_points[i] for i in ticks]
plt.xticks(ticks, labels, rotation=45)

plt.title(f'Adherence Patterns for {best_model} Clusters (Time-Series K-Means)')
plt.ylabel('Mean Adherence (Standardized)')
plt.xlabel('Time Point')
plt.legend(title="Cluster")
plt.tight_layout()
plt.show()

# ---------------------- Hidden Markov Model ----------------------
# How can we assess changes in individual characteristics when these are not directly observable?
# 1 known categorical variable: CPAP adherence with 2 hidden states
# The last 200 patients and all time points were included

# Filter for patient IDs from 100 to 300
selected_ids = range(100, 301)
df = Sim_CPAP_cat[Sim_CPAP_cat['patient_id'].isin(selected_ids)].copy()

# Convert categorical adherence values to numerical with LabelEncoder
obs_cols = [col for col in df.columns if col.startswith('T')]
for col in obs_cols:
    le = LabelEncoder()
    df[col] = le.fit_transform(df[col])

# Reshape to long format for time series modeling
df_long = df.melt(id_vars='patient_id', value_vars=obs_cols,
                  var_name='Time', value_name='CPAP_adherence')
df_long = df_long.sort_values(['patient_id', 'Time'])

# Get sequence lengths (90 time points per patient)
ntimes = df_long.groupby('patient_id').size().values

# Prepare observed sequences
X = df_long['CPAP_adherence'].values.reshape(-1, 1)

# Fit HMM with 2 hidden states
model = CategoricalHMM(n_components=2, n_iter=100, random_state=42)
model.fit(X, lengths=ntimes)

# Predict hidden states
pred_states = model.predict(X)

# Add predictions to dataframe for plotting
df_long['State'] = pred_states

# Option 1 - Plot first 200 points (all patients)
plot_df = df_long.iloc[:200]
# Option 2 -  shown only 1 patient:
# plot_df = df_long[df_long['patient_id'] == df['patient_id'].sample(1).values[0]].iloc[:200] --> show only 1 patient

plt.figure(figsize=(12, 4))
sns.lineplot(data=plot_df, x=plot_df.index, y='State')
plt.title('Predicted Hidden States (First 200 Time Points)')
plt.xlabel('Time Points')
plt.ylabel('State Prediction')
plt.yticks([0, 1], ['Non-adherent', 'Adherent'])
plt.show()

# Print initial state probabilities and transition matrix
init_prob = model.startprob_
trans_prob = model.transmat_

print("Initial state probabilities:\n", init_prob)
print("Transition matrix:\n", trans_prob)


# ---------------------- Gaussian Mixture Model ----------------------
# GMM: How can we describe trajectories of longitudinal data with repeated measurements of follow-up?
# Continuous data for CPAP adherence
# All patients are included, but only 5 time points

# Select data (first 5 time points)
Sim_CPAP = Sim_CPAP[['patient_id'] + [f'T{i}' for i in range(1, 6)]].copy()

# Step 1: Data preparation
df_long = Sim_CPAP.melt(id_vars=["patient_id"], value_vars=[f'T{i}' for i in range(1, 6)],
                      var_name="Time", value_name="CPAP_adherence")
df_long['Time'] = df_long['Time'].str.extract('(\d+)').astype(int)

# Step 2: Fit GMMs for 1-4 clusters
X = df_long[['Time', 'CPAP_adherence']]
models = {}
bic_scores = []

for n in range(1, 5):
    gmm = GaussianMixture(n_components=n, random_state=42)
    gmm.fit(X)
    models[n] = gmm
    bic_scores.append(gmm.bic(X))

# Step 3: Choose best model by BIC
best_n = np.argmin(bic_scores) + 1
best_model = models[best_n]
df_long['Cluster'] = best_model.predict(X)

# Step 4: Plot CPAP adherence over time with clusters
plt.figure(figsize=(10, 6))
sns.lineplot(data=df_long, x='Time', y='CPAP_adherence', hue='Cluster',
             estimator='mean', errorbar=None, palette='tab10')
plt.title(f'CPAP Adherence Over Time (Best GMM: {best_n} Clusters)')
plt.xlabel('Time')
plt.ylabel('CPAP Adherence')
plt.legend(title='Cluster')
plt.tight_layout()
plt.show()

# Step 5: Plot posterior probabilities distribution
posterior_probabilities = best_model.predict_proba(X)
df_long['Posterior'] = posterior_probabilities.max(axis=1)

plt.figure(figsize=(10, 6))
sns.histplot(df_long['Posterior'], bins=30, kde=True)
plt.title('Distribution of Maximum Posterior Probabilities')
plt.xlabel('Posterior Probability')
plt.ylabel('Frequency')
plt.tight_layout()
plt.show()

# ---------------------- ARIMA & CCF method ----------------------
# How to analyze time series and evaluate the correlation between 2 time series varying over time, 
# coinciding or not over time intervals?
# Numerical outcome for time series
# 2 time series per patient: 1 for CPAP adherence and 1 for ESS score
# Then, CCF is performed between both time series
# All time points are included, but only 1 patient

# Simulate new ESS data
np.random.seed(42)
Sim_ESS_new = pd.DataFrame(
    np.random.randint(4, 19, size=(300, 90)),
    columns=[f'T{i}' for i in range(1, 91)]
)
Sim_ESS_new['patient_id'] = range(1, 301)

# Filter the data to keep 1 patient and transform to long format
cpap_ARIMA = Sim_CPAP[Sim_CPAP['patient_id'] == 10].melt(id_vars='patient_id', var_name='Time', value_name='CPAP_adherence')
cpap_ARIMA['Time'] = cpap_ARIMA['Time'].str.extract(r'T(\d+)')[0].astype(float)

ess_ARIMA = Sim_ESS_new[Sim_ESS_new['patient_id'] == 10].melt(id_vars='patient_id', var_name='Time', value_name='ESS_score')
ess_ARIMA['Time'] = ess_ARIMA['Time'].str.extract(r'T(\d+)')[0].astype(float)

# Transform to time series format
ts_Cpap = cpap_ARIMA.sort_values('Time')['CPAP_adherence'].reset_index(drop=True)
ts_Ess = ess_ARIMA.sort_values('Time')['ESS_score'].reset_index(drop=True)

# ADF test to check stationarity
print("ADF Test CPAP (p-value):", adfuller(ts_Cpap)[1])
print("ADF Test ESS (p-value):", adfuller(ts_Ess)[1])

# ACF & PACF plots
# ACF -> autocorrelation = correlation between the data and themselves
# PACF -> partial autocorrelation = at lag k, correlation between all data points that are exactly k steps apart
# after accounting for their correlation with data between k steps
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
plot_acf(ts_Cpap, ax=axes[0])
plot_pacf(ts_Cpap, ax=axes[1])
axes[0].set_title("ACF - CPAP")
axes[1].set_title("PACF - CPAP")
plt.tight_layout()
plt.show()

# Fit auto ARIMA and print results
model_cpap = auto_arima(ts_Cpap, seasonal=False, stepwise=True, suppress_warnings=True)
model_ess = auto_arima(ts_Ess, seasonal=False, stepwise=True, suppress_warnings=True)

print(model_cpap.summary())
print(model_ess.summary())

# Validation tests
# Residuals
res_cpap = pd.Series(model_cpap.resid())
res_ess = pd.Series(model_ess.resid())

# QQ plot
plt.figure(figsize=(6, 6))
qqplot(res_cpap, line='s', ax=plt.gca())
plt.title("QQ plot of CPAP ARIMA residuals")
plt.tight_layout()
plt.show()
# -> no pattern in the residuals => ok

# Ljung-Box test
ljungbox_results = acorr_ljungbox(res_cpap, lags=[20], return_df=True)
print("Ljung-Box test results:\n", ljungbox_results)
# p>0.05 -> result is not significant

# Cross-correlation function (CCF)
# CCF between ESS and CPAP variables modified (detrend) by the ARIMA function
ccf_vals = ccf(res_cpap, res_ess)
plt.figure(figsize=(10, 5))
plt.stem(range(len(ccf_vals)), ccf_vals)
plt.title("Cross-correlation function (CCF): CPAP vs ESS")
plt.xlabel("Lag")
plt.ylabel("Correlation")
plt.show()

# Lagged regression at lag 14
cpap_lags = lagmat(res_cpap, maxlag=14, trim='both')
final_df = pd.DataFrame({
    'ESS_resid': res_ess[-len(cpap_lags):].reset_index(drop=True),
    'CPAP_lag14': cpap_lags[:, 13]
})
X = add_constant(final_df['CPAP_lag14'])
y = final_df['ESS_resid']
reg = OLS(y, X).fit()
print(reg.summary())

# ACF of regression residuals
plot_acf(reg.resid, lags=20)
plt.title("ACF of Regression Residuals")
plt.show()

# Possible to add covariates to the ARIMA method --> ARIMAX 

# ---------------------- Mixed model ----------------------
# How to estimate the relationship between dependent variables and fixed and random effects of the independents variables?
# Continuous or categorical variable supported - we use the continuous outcome in this example:
# Time and ESS score as categorical variables
# All patients and all time points are included

# Reshape CPAP adherence data to long format
cpap_long = Sim_CPAP.melt(
    id_vars=["patient_id"],
    value_vars=[f"T{i}" for i in range(1, 91)],
    var_name="Time",
    value_name="CPAP_adherence"
)

# Extract numeric 'time' and convert to category (ordered)
cpap_long["Time"] = cpap_long["Time"].str.extract(r"(\d+)").astype(int)

# Ensure 'patient_id' is categorical
cpap_long["patient_id"] = cpap_long["patient_id"].astype("category")

# Prepare ESS_baseline (T1)
ess_baseline = Sim_ESS[["patient_id", "T1"]].rename(columns={"T1": "ESS_baseline"})
ess_baseline["ESS_baseline"] = ess_baseline["ESS_baseline"].astype("category")

# Merge CPAP data with ESS baseline
cpap_long = cpap_long.merge(ess_baseline, on="patient_id")

# Fit mixed effects model with random intercept and slope for 'Time'
model = smf.mixedlm(
    "CPAP_adherence ~ Time + ESS_baseline",
    data=cpap_long,
    groups="patient_id",
    re_formula="~Time"  # random intercept + slope for Time
)

result = model.fit()
print(result.summary())

# Plot: CPAP adherence over time by ESS_baseline
cpap_long["fitted"] = result.fittedvalues
plt.figure(figsize=(8,6))
sns.lineplot(data=cpap_long,x="Time", y="fitted", hue="ESS_baseline",estimator="mean", ci=None)
plt.title("CPAP Adherence Over Time by ESS Baseline")
plt.xlabel("Time (days)")
plt.ylabel("CPAP Adherence")
plt.show()

# ---------------------- Survival analysis ----------------------
# How to determine the survival probability or cumulative risk of a population over a defined period of time?
# Continuous CPAP adherence and categorical ESS scores for Cox model
# All patients and all time points were included

# The survival event in this case is the moment when the CPAP adherence reaches >= 4h
# Time = max between time when CPAP adherence = 1 and T90 (the end of the study)
# The survival event has occured (CPAP_adherence = 1) when CPAP adherence >= 4h before the end of the study or 
# the patient is censored (CPAP_adherence = 0) when CPAP adherence < 4h over the study

# Convert Sim_CPAP to long format
cpap_long = Sim_CPAP.melt(
    id_vars=["patient_id"],
    value_vars=[f"T{i}" for i in range(1, 91)],
    var_name="Time",
    value_name="CPAP_adherence")

# Convert Time to numeric
cpap_long["Time"] = cpap_long["Time"].str.extract(r"(\d+)").astype(int)

# Binary event: 1 if CPAP_adherence >= 4, else 0
cpap_long["event"] = (cpap_long["CPAP_adherence"] >= 4).astype(int)

# Group the survival data by patient
survival_df = (
    cpap_long
    .groupby("patient_id")
    .agg({
        "event": "max",
        "Time": lambda x: cpap_long.loc[x.index, "event"].idxmax() if x.any() else 90
    })
    .reset_index()
)

# Replace index of time with actual time values
survival_df["Time"] = survival_df["Time"].apply(lambda i: cpap_long.loc[i, "Time"] if i in cpap_long.index else 90)

# Set event=0 if patient never reached CPAP >= 4 by T90
survival_df["event"] = survival_df["Time"].apply(lambda t: 0 if t == 90 else 1)

# Join with Drowsy (ESS >= 10) from Sim_ESS_cat
ess_df = Sim_ESS_cat[["patient_id", "T1"]].rename(columns={"T1": "Drowsy"})
ess_df["Drowsy"] = ess_df["Drowsy"].map({"Yes": 1, "No": 0})

# Merge on survival_df
survival_df = survival_df.merge(ess_df, on="patient_id")

# Fit Kaplan-Meier model
kmf = KaplanMeierFitter()

# Plot KM curves for Drowsy and Non-Drowsy
plt.figure(figsize=(10, 6))

for label, group_df in survival_df.groupby("Drowsy"):
    kmf.fit(durations=group_df["Time"], event_observed=group_df["event"], label=f"Drowsy={label}")
    kmf.plot(ci_show=True)

plt.title("Survival: Time to CPAP Adherence ≥ 4h by Drowsiness (ESS)")
plt.xlabel("Time (days)")
plt.ylabel("Survival Probability")
plt.grid(True)
plt.tight_layout()
plt.show()

# Log-rank test for comparison
drowsy0 = survival_df[survival_df["Drowsy"] == 0]
drowsy1 = survival_df[survival_df["Drowsy"] == 1]

results = logrank_test(drowsy0["Time"], drowsy1["Time"], event_observed_A=drowsy0["event"], event_observed_B=drowsy1["event"])
results.print_summary()
# p>0.05 => no statistically significant difference in time to CPAP adherence between the drowsy and non-drowsy patient groups
