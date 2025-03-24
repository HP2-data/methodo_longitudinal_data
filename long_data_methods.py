#Longitudinal data, trajectories and time series: how to analyze them? An example of sleep data

# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import mixedlm, smf
from statsmodels.stats.anova import AnovaRM 
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.stattools import adfuller, ccf
from statsmodels.tsa.tsatools import lagmat
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.stats.diagnostic import acorr_ljungbox
from tslearn.clustering import TimeSeriesKMeans
from tslearn.utils import to_time_series_dataset
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import calinski_harabasz_score
from sklearn.cluster import KMeans
from hmmlearn.hmm import MultinomialHMM

# ---------------------- Load simulated data ----------------------
def find_csv_file(file_name, start_directory):
    for root, dirs, files in os.walk(start_directory):  # Start from the Downloads folder
        if file_name in files:
            return os.path.join(root, file_name)
    return None

# Replace this path with your actual Downloads directory path
downloads_folder = os.path.expanduser('~/Downloads')

# Searching for 'name.csv' in the Downloads folder
csv_file_path = find_csv_file('Sim_CPAP.csv', downloads_folder)


if csv_file_path:
    # If file is found, load it into a DataFrame
    df = pd.read_csv(csv_file_path)
    print("CSV file loaded successfully!")
    print(df.head())
else:
    print("CSV file not found.")


cpap = pd.read_csv('/Users/jandrici/Desktop/methods/Data/new/Sim_CPAP.csv', sep=';')
cpap_cat = pd.read_csv('/Users/jandrici/Desktop/methods/Data/new/Sim_CPAP_cat.csv', sep=';')
ess = pd.read_csv('/Users/jandrici/Desktop/methods/Data/new/Sim_ESS.csv', sep=';')
ess_cat = pd.read_csv('/Users/jandrici/Desktop/methods/Data/new/Sim_ESS_cat.csv', sep=';')

# ---------------------- ANOVA method -----------------------
# Convert cpap dataframe from wide to long-format
df_long = cpap.melt(id_vars=['patient_id'], var_name='time', value_name='value')
anova_results = AnovaRM(data=df_long, depvar='value',
                  subject='patient_id', within=['time']).fit()
p_value = anova_results.anova_table['Pr > F'][0]
# print(f"ANOVA method, p-value: {p_value}")

# ---------------------- Mantel-Haenszel method ----------------------
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
ess_cat[['T1', 'T2', 'T3', 'T4']] = ess_cat[['T1', 'T2', 'T3', 'T4']].apply(lambda col: col.map(ess_map))
cpap_cat[['T1', 'T2', 'T3', 'T4']] = cpap_cat[['T1', 'T2', 'T3', 'T4']].apply(lambda col: col.map(adherence_map))

# Create contingency tables for each time point (T1, T2, T3, T4)
tables = []
for time in ['T1', 'T2', 'T3', 'T4']:
    table = pd.crosstab(ess_cat[time], cpap_cat[time])
    # print(f"\nContingency table for {time} (Shape: {table.shape}):")
    # print(table)
    tables.append(table.values)

# Convert to 2x2xk format (transpose ensures correct shape)
contingency_tables = np.stack(tables, axis=2)
# print(contingency_tables.shape)  # Should be (2, 2, 4)

# Create StratifiedTable using statsmodels package
st = smstats.StratifiedTable(contingency_tables.astype(np.float64))

# Get the summary of Mantel-Haenszel test results
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
# Copy dataset to avoid modifying the original
df = cpap.copy()

# Select patient ID and time points
time_points = ['T1', 'T2', 'T3', 'T4', 'T5']
df = df[['patient_id'] + time_points]

# Standardize time series data
scaler = StandardScaler()
df[time_points] = scaler.fit_transform(df[time_points])

# Convert to time-series format required for TimeSeriesKMeans
ts_data = to_time_series_dataset(df[time_points].values)

# Define range for number of clusters (2 to 6)
n_clusters_range = range(2, 7)
calinski_scores = []
models = {}

# Fit TimeSeriesKMeans for different cluster numbers
for n_clusters in n_clusters_range:
    model = TimeSeriesKMeans(n_clusters=n_clusters, metric="euclidean", random_state=42, n_init=10)
    labels = model.fit_predict(ts_data)
    score = calinski_harabasz_score(df[time_points], labels)
    
    calinski_scores.append(score)
    models[n_clusters] = model

# Identify the best number of clusters
optimal_clusters = n_clusters_range[np.argmax(calinski_scores)]
# print(f"Optimal number of clusters: {optimal_clusters}")

# Fit final TimeSeriesKMeans model with optimal clusters
best_model = models[optimal_clusters]
df['cluster'] = best_model.predict(ts_data)

# Display cluster distribution
print('_____________________')
print('Cluster distribution:')
print(df['cluster'].value_counts(normalize=True))

# Reshape for visualization (long format)
df_melted = df.melt(id_vars=['patient_id', 'cluster'],
                     value_vars=time_points,
                     var_name='Time Point',
                     value_name='Adherence')

# Plot adherence patterns per cluster
plt.figure(figsize=(10, 6))
sns.lineplot(data=df_melted, x='Time Point', y='Adherence',
             hue='cluster', estimator=np.mean, errorbar=None, palette='viridis')
plt.title(f'Adherence Patterns for {optimal_clusters} Clusters (Time-Series K-Means)')
plt.ylabel('Mean Adherence (Standardized)')
plt.legend(title="Cluster")
plt.show()

# ---------------------- Hidden Markov Model ----------------------
# Filter the dataset with categorical adherence for patient with ID 25 (example)
patient_id = 25
cpap_cat_patient = cpap_cat[cpap_cat['patient_id'] == patient_id].drop(columns=['patient_id']).T

# Convert categorical data to numeric (mapping categories)
category_map = {'High': 0, 'Low': 1, '≥ 4h': 2, '[2h,4h[': 3, '[0h,2h[': 4}
cpap_cat_patient = cpap_cat_patient.applymap(category_map.get)

# Reshape data to a suitable format for HMM (one column per observation)
X = cpap_cat_patient.values.reshape(-1, 1)

# Create and fit the HMM model with 2 hidden states (adherent and non-adherent)
model = MultinomialHMM(n_components=2, n_iter=1000, random_state=42)
model.fit(X)

# Predict the hidden states
pred_states = model.predict(X)

# Evaluate the model: Log-Likelihood
log_likelihood = model.score(X)

# Calculate number of parameters (startprob, transmat, emissions)
n_params = model.n_components + model.n_components * model.n_components + model.n_components * len(np.unique(X))

# Calculate AIC and BIC
aic = 2 * n_params - 2 * log_likelihood
bic = np.log(len(X)) * n_params - 2 * log_likelihood

print(f"AIC: {aic}")
print(f"BIC: {bic}")
print(f"Log-Likelihood: {log_likelihood}")

# Extract the initial state probabilities
init_prob = model.startprob_
print("Initial State Probabilities:")
print(init_prob)

# Extract the transition probabilities (state-to-state transition matrix)
trans_prob = model.transmat_
print("Transition Probabilities:")
print(trans_prob)

# Plot the predicted hidden states over time
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(pred_states) + 1), pred_states, marker='o', linestyle='-', color='b')
plt.yticks([0, 1], ['Non adherent', 'Adherent'])
plt.xlabel('Time points')
plt.ylabel('State prediction')
plt.title('Hidden States Prediction Over Time')
plt.show()

# Frequency of the predicted states (how many times each state occurs)
state_freq = pd.Series(pred_states).value_counts()
print("State Frequencies:")
print(state_freq)

# ---------------------- Gaussian Mixture Model ----------------------
# GMM: How can we describe trajectories of longitudinal data with repeated measurements of follow-up?
# Continuous data for CPAP adherence
# All patients are included, but only 5 time points

# Data preparation (pivot long format for the CPAP adherence data)
cpap_long = cpap.melt(id_vars=["patient_id"], value_vars=["T1", "T2", "T3", "T4", "T5"], 
                                var_name="Time", value_name="CPAP_adherence")

# Convert 'Time' column to numeric
cpap_long['Time'] = cpap_long['Time'].str.extract('(\d+)').astype(int)

# Prepare data for Gaussian Mixture Model (X: patient_id, CPAP_adherence, Time)
X = cpap_long[['Time', 'CPAP_adherence']]

# Step 1: Fit GMM with different numbers of clusters (n = 1, 2, 3, 4)
models = {}
bic_scores = []

for n in range(1, 5):  # Trying 1 to 4 clusters
    gmm = GaussianMixture(n_components=n, random_state=42)
    gmm.fit(X)
    models[n] = gmm
    bic_scores.append(gmm.bic(X))

# Find the best model based on the lowest BIC score
best_n = np.argmin(bic_scores) + 1
best_model = models[best_n]

# Print the BIC scores and the best model
print("BIC scores for different number of clusters:")
for n, bic in zip(range(1, 5), bic_scores):
    print(f"n = {n}: BIC = {bic}")

print(f"Best model: {best_n} clusters (based on BIC)")

# Step 2: Get posterior probabilities (cluster membership)
posterior_probabilities = best_model.predict_proba(X)
cpap_long['Cluster'] = np.argmax(posterior_probabilities, axis=1)

# Step 3: Visualize the clusters using a plot of CPAP adherence over Time
plt.figure(figsize=(10, 6))
sns.scatterplot(data=cpap_long, x='Time', y='CPAP_adherence', hue='Cluster', palette='Set1')
plt.title('CPAP Adherence Over Time (Clustered)')
plt.xlabel('Time')
plt.xticks(np.arange(0,6, step=1))
plt.ylabel('CPAP Adherence')
plt.legend(title='Cluster')
plt.show()

# Step 4: Extract fixed effects and plot
cluster_means = best_model.means_
cluster_covariances = best_model.covariances_

# Create a DataFrame with cluster means and covariances for visualization
fixed_effects = pd.DataFrame({
    'Cluster': [f'Cluster {i+1}' for i in range(best_n)],
    'Mean': cluster_means[:, 0],  # Use the correct index for CPAP adherence (the first column)
    'Covariance': cluster_covariances[:, 0, 0],  # Covariance for CPAP adherence (first index for CPAP)
})

# Plot fixed effects
plt.figure(figsize=(10, 6))
sns.pointplot(x='Cluster', y='Mean', data=fixed_effects, markers='o', color='b')
plt.title('Fixed Effects in the Longitudinal GMM Model')
plt.ylabel('Mean CPAP Adherence')
plt.xlabel('Cluster')
plt.show()

# Plot the posterior probabilities (cluster membership probabilities) for each patient
plt.figure(figsize=(10, 6))
for i in range(best_n):
    sns.kdeplot(cpap_long[cpap_long['Cluster'] == i]['CPAP_adherence'], label=f'Cluster {i+1}')

plt.title('Distribution of CPAP Adherence by Cluster')
plt.xlabel('CPAP Adherence')
plt.ylabel('Density')
plt.legend(title='Cluster')
plt.show()

# ---------------------- ARIMA & CCF method ----------------------
# ARIMA & CCF: how can we analyze time series and evaluate the correlation between 2 time series varying over time, 
# coinciding or not over time intervals?
# Numerical outcome for time series
# 2 time series per patient: 1 for CPAP adherence and 1 for ESS score
# Then, CCF is performed between both time series
# All time points are included, but only 1 patient

# Data Preparation (Reshape to Long Format)
cpap_ARIMA = cpap[['patient_id'] + [f'T{i}' for i in range(1, 91)]].loc[cpap['patient_id'] == 10]
cpap_ARIMA = cpap_ARIMA.melt(id_vars=['patient_id'], var_name='Time', value_name='CPAP_adherence')
cpap_ARIMA['Time'] = cpap_ARIMA['Time'].str.extract('(\d+)').astype(int)

ess_ARIMA = ess[['patient_id'] + [f'T{i}' for i in range(1, 91)]].loc[ess['patient_id'] == 10]
ess_ARIMA = ess_ARIMA.melt(id_vars=['patient_id'], var_name='Time', value_name='ESS_score')
ess_ARIMA['Time'] = ess_ARIMA['Time'].str.extract('(\d+)').astype(int)

# Convert CPAP and ESS data to time series with frequency 7 (weekly data)
ts_Cpap = pd.Series(cpap_ARIMA['CPAP_adherence'].values, index=pd.date_range(start='1/1/2000', periods=len(cpap_ARIMA), freq='D'))
ts_Ess = pd.Series(ess_ARIMA['ESS_score'].values, index=pd.date_range(start='1/1/2000', periods=len(ess_ARIMA), freq='D'))

# Check Stationarity using ADF Test
adf_result_Cpap = adfuller(ts_Cpap)
adf_result_Ess = adfuller(ts_Ess)

print(f"ADF Test for CPAP: p-value = {adf_result_Cpap[1]}")
print(f"ADF Test for ESS: p-value = {adf_result_Ess[1]}")

# Plot the ACF and PACF for CPAP
plt.figure(figsize=(12, 6))
plt.subplot(121)
plot_acf(ts_Cpap, lags=20, ax=plt.gca())
plt.title('ACF of CPAP Adherence')

plt.subplot(122)
plot_pacf(ts_Cpap, lags=20, ax=plt.gca())
plt.title('PACF of CPAP Adherence')
plt.tight_layout()
plt.show()

# ARIMA application
# Fit ARIMA models for both CPAP and ESS
ARIMA_Cpap = ARIMA(ts_Cpap, order=(0, 0, 0), seasonal_order=(1, 0, 1, 7)).fit()
print(f"ARIMA CPAP Summary:\n{ARIMA_Cpap.summary()}")

ARIMA_Ess = ARIMA(ts_Ess, order=(0, 0, 0), seasonal_order=(0, 0, 0, 7)).fit()
print(f"ARIMA ESS Summary:\n{ARIMA_Ess.summary()}")

# Residuals from the ARIMA models
Cpap_residuals = ARIMA_Cpap.resid
Ess_residuals = ARIMA_Ess.resid

# QQ plot for CPAP residuals
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
sns.histplot(Cpap_residuals, kde=True)
plt.title('QQ Plot for ARIMA CPAP Residuals')
plt.show()

# Ljung-Box test for CPAP residuals
lb_test_Cpap = acorr_ljungbox(Cpap_residuals, lags=[20], return_df=True)
print(f"Ljung-Box Test for CPAP residuals:\n{lb_test_Cpap}")

# CCF (Cross-Correlation Function) between CPAP and ESS residuals
ccf_result = ccf(Cpap_residuals, Ess_residuals, adjusted=True)

# Plot the CCF
plt.figure(figsize=(10, 6))
plt.plot(ccf_result)
plt.title('CCF between CPAP and ESS Residuals')
plt.xlabel('Lag')
plt.ylabel('CCF')
plt.show()

# Lag plot for CPAP and ESS residuals with lag of 10
plt.figure(figsize=(10, 6))
lag = 10
plt.plot(Cpap_residuals[:-lag], Ess_residuals[lag:], 'o', label='Lag 10')
plt.title(f'Lag plot (Lag = {lag})')
plt.xlabel('CPAP Residuals')
plt.ylabel('ESS Residuals')
plt.legend()
plt.show()

# Regression of ESS residuals with CPAP residuals at different lags
# Lag CPAP residuals
Cpap_lagged = lagmat(Cpap_residuals, maxlag=14, trim='both')

# Create final dataset for regression
final_data = pd.DataFrame({
    'ESS_residuals': Ess_residuals[14:],  # Adjust length to match lagged data
    'Cpap_lag11': Cpap_lagged[:, 10],  # Lag 11 (index 10)
    'Cpap_lag1': Cpap_lagged[:, 0],  # Lag 1 (index 0)
    'Cpap_lag14': Cpap_lagged[:, 13]  # Lag 14 (corrected index)
})

# Fit linear regression model
X = sm.add_constant(final_data[['Cpap_lag11', 'Cpap_lag1', 'Cpap_lag14']])
y = final_data['ESS_residuals']
reg_model = sm.OLS(y, X).fit()
print(f"Regression Summary:\n{reg_model.summary()}")

# ACF of residuals of regression model
plot_acf(reg_model.resid, lags=20)
plt.title('ACF of Regression Residuals')
plt.show()

# ---------------------- Mixed model ----------------------
#Mixed model: how can we estimate the relationship between the dependent variables and 
# the fixed and random effects of the independents variables?
# Continuous or categorical variable supported - we use the continuous outcome in this example:
# Time and ESS score as categorical variables
# All patients and all time points are included

# Trnasform CPAP adherence dataset to long format
cpap_long = cpap.melt(id_vars=["patient_id"], value_vars=[f"T{i}" for i in range(1, 91)],
                       var_name="Time", value_name="CPAP_adherence")

# Convert 'Time' and 'patient_id' to categorical variables
cpap_long['Time'] = cpap_long['Time'].str.extract('(\d+)').astype(int).astype("category")
cpap_long['patient_id'] = cpap_long['patient_id'].astype("category")

# Get ESS baseline score (T1)
ess_baseline = ess[['patient_id', 'T1']].rename(columns={"T1": "ESS_baseline"})
ess_baseline['ESS_baseline'] = ess_baseline['ESS_baseline'].astype("category")

# Merge CPAP data with ESS baseline
cpap_long = cpap_long.merge(ess_baseline, on="patient_id")

# Fit the Mixed Effects Model (random intercept per patient)
mixed_model = smf.mixedlm("CPAP_adherence ~ Time + ESS_baseline", 
                          data=cpap_long, groups=cpap_long["patient_id"])
mixed_result = mixed_model.fit()

# Print model summary
print(mixed_result.summary())

# Convert 'Time' to numeric for plotting
cpap_long["Time"] = cpap_long["Time"].astype(int)

# Plot the results (optional)
sns.lmplot(data=cpap_long, x="Time", y="CPAP_adherence", hue="ESS_baseline", aspect=1.5)
plt.title("CPAP Adherence Over Time by ESS Baseline")
plt.show()
