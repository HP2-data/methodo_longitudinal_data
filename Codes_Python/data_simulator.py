### Data Simulation function

# Imports
import numpy as np
import pandas as pd

# Simulate CPAP adherence data with 3 trajectory groups and time windows
def sim_data(nb_patient, nb_time_point):
    sim_df = np.zeros((nb_patient, nb_time_point))
    
    group_size = nb_patient // 3
    groups = {
        0: range(0, group_size),
        1: range(group_size, 2 * group_size),
        2: range(2 * group_size, nb_patient)
    }

    time_point = nb_time_point // 3
    time_windows = {
        0: range(0, time_point),
        1: range(time_point, 2 * time_point),
        2: range(2 * time_point, nb_time_point)
    }

    means = [2.2, 4.8, 6.3]
    sds = [1.0, 0.6, 0.7]

    # Rotate groups through time windows
    for j in range(3):
        k = (j + 1) % 3
        l = (j + 2) % 3
        for i in groups[j]:
            np.random.seed(i + 31)
            sim_df[i, list(time_windows[j])] = np.random.normal(means[j], sds[j], len(time_windows[j]))
        for i in groups[k]:
            np.random.seed(i + 31)
            sim_df[i, list(time_windows[j])] = np.random.normal(means[k], sds[k], len(time_windows[j]))
        for i in groups[l]:
            np.random.seed(i + 30)
            sim_df[i, list(time_windows[j])] = np.random.normal(means[l], sds[l], len(time_windows[j]))

    df = pd.DataFrame(sim_df, columns=[f'T{i+1}' for i in range(nb_time_point)])
    df['patient_id'] = np.arange(1, nb_patient + 1)
    return df

# Simulate ESS score data (categorical values from 0 to 24)
def sim_data_discrete(nb_patient, nb_time_point, score_max):
    sim_df = np.zeros((nb_patient, nb_time_point), dtype=int)

    for i in range(nb_patient):
        np.random.seed(i + 31)
        probs = np.random.beta(9, 15, nb_time_point)
        sim_df[i] = np.random.binomial(score_max, probs)

    df = pd.DataFrame(sim_df, columns=[f'T{i+1}' for i in range(nb_time_point)])
    df['patient_id'] = np.arange(1, nb_patient + 1)
    return df

# Example usage
# nb_patient = 10
# nb_time_point = 6
# score_max = 24

# cpap_data = sim_data(nb_patient, nb_time_point)
# ess_data = sim_data_discrete(nb_patient, 6, score_max)

# print("CPAP Adherence Data:")
# print(cpap_data)
# print("\nESS Score Data:")
# print(ess_data)
