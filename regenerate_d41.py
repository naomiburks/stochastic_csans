import numpy as np
from scipy.integrate import odeint 
from stochastic_csans import State, CSANModel
from constants import SIMPLE_D41_PARAMETERS as params
from plot import plot_all_aggregate_normalized
import csv

data = []

model = CSANModel(4, 8)

initial_state = model.get_empty_state()
initial_state.D = 0
initial_state.Es[0] = 0
initial_state.Ts[0] = 1000
timepoints = [i * 4 for i in range(13)]

result = model.generate_simulation_data(params, initial_state, timepoints, sample_count=3)["data"]

for sample in result:
    row = []
    for timepoint in timepoints:
        state = sample[timepoint]
        tumor_count = sum(state.Ts) + sum([sum(l) for l in state.ETs]) + sum([sum(l) for l in state.EDTs])
        row.append(tumor_count)
    data.append(row)


initial_state.Es[0] = 3000
result = model.generate_simulation_data(params, initial_state, timepoints, sample_count=3)["data"]
for sample in result:
    row = []
    for timepoint in timepoints:
        state = sample[timepoint]
        tumor_count = sum(state.Ts) + sum([sum(l) for l in state.ETs]) + sum([sum(l) for l in state.EDTs])
        row.append(tumor_count)
    data.append(row)

for i in range(6):
    print(i)
    initial_state.D = 4 ** i * 40
    print(initial_state.D)
    result = model.generate_simulation_data(params, initial_state, timepoints, sample_count=3)["data"]
    for sample in result:
        row = []
        for timepoint in timepoints:
            state = sample[timepoint]
            tumor_count = sum(state.Ts) + sum([sum(l) for l in state.ETs]) + sum([sum(l) for l in state.EDTs])
            row.append(tumor_count)
        data.append(row)


with open('data/D41Sim.csv', 'w', newline='') as csvfile:
    datawriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    datawriter.writerows(data)