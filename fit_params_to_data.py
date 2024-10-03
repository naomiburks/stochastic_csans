from stochastic_csans import CSANModel
from scipy import optimize
from constants import SIMPLE_D41_PARAMETERS
import csv


model = CSANModel(4, 8)
param_names = list(SIMPLE_D41_PARAMETERS.keys())
timepoints = [i * 4 for i in range(13)]


with open('data/D41_Ave.csv', newline='', ) as csvfile:
    data_reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    data = []
    for row in data_reader:
        data.append([float(num) for num in row])

def fit_to_data(data, guess, initial_condition):
    pass

def _sample_error(parameters, data):
    parameters = get_param_dict(parameters)
    print(parameters)
    error = 0
    for sim_index in range(8):
        initial_state = model.get_empty_state()
        initial_state.Ts[0] = 100
        if sim_index != 0:
            initial_state.Es[0] = 300
        if sim_index >= 2:
            initial_state.D = 4**(8 - sim_index) * 20
        
        sample = model.generate_simulation_data(parameters, initial_state, timepoints)["data"][0]
        real_data = [data[i][sim_index] for i in range(8)]
        normalized_sample = [sample[timepoint].T() / 1000 for timepoint in timepoints]

        for real_point, sample_point in zip(real_data, normalized_sample):
            error += (real_point - sample_point) ** 2
    print(error)
    return error

def get_param_dict(param_list):
    return {name: param for name, param in zip(param_names, param_list)}

def get_param_list(param_dict):
    return [param_dict[name] for name in param_names]

if __name__ == "__main__":
    initial_state = model.get_empty_state()
    
    initial_guess = get_param_list(SIMPLE_D41_PARAMETERS)
    bounds = [(0, guess * 2) for guess in initial_guess]
    opt_res = optimize.minimize(_sample_error, initial_guess, args=data, bounds=bounds)
    print(opt_res)
    

    
