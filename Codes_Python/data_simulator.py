import numpy as np
import pandas as pd

# Simulating data for CPAP adherence (Normal positive data)
def sim_data(nb_patient, nb_time_point):
    np.random.seed(42) 
    
    # Create the matrix with random normal values (mean=4, std=1.5)
    sim_matrix = np.random.normal(4, 1.5, (nb_patient, nb_time_point))
    
    # Convert the matrix into a dataframe
    sim_cpap = pd.DataFrame(sim_matrix)
    
    # Add patient_id as the last column
    sim_cpap['patient_id'] = np.arange(1, nb_patient + 1)
    
    # Rename the columns
    time_columns = [f'T{i}' for i in range(1, nb_time_point + 1)]
    sim_cpap.columns = time_columns + ['patient_id']
    
    return sim_cpap

# Simulating data for ESS score (categorical values from 0 to score_max)
def sim_data_discrete(nb_patient, nb_time_point, score_max):
    np.random.seed(42) 
    
    # Create the matrix with random categorical values in range(0, score_max + 1)
    sim_matrix = np.random.randint(0, score_max + 1, size=(nb_patient, nb_time_point))
    
    # Convert the matrix into a DataFrame
    sim_ess = pd.DataFrame(sim_matrix)
    
    # Add patient_id as the last column
    sim_ess['patient_id'] = np.arange(1, nb_patient + 1)
    
    # Rename the columns
    time_columns = [f'T{i}' for i in range(1, nb_time_point + 1)]
    sim_ess.columns = time_columns + ['patient_id']
    
    return sim_ess

# Example usage
nb_patient = 10  # Number of patients
nb_time_point = 5  # Number of time points
score_max = 24  # Max ESS score

cpap_data = sim_data(nb_patient, nb_time_point)
ess_data = sim_data_discrete(nb_patient, nb_time_point, score_max)

# Display the results (optional)
print("CPAP Adherence Data:")
print(cpap_data)
print("\nESS Score Data:")
print(ess_data)
