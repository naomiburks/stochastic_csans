from models.stochastic_csans import CSANModel
from plot import make_average_free_csans_plot, make_average_TCell_plot, make_average_tumor_cell_plot, make_total_tumor_cell_plot, plot_all_aggregate_normalized
import constants




duration = 48
timesteps = 480
times_to_plot = [duration * i * 120 / timesteps for i in range(5)]

model = CSANModel(10, 10)
state = model.get_empty_state()
state.D = 10000
state.Es[0] = 100
state.Ts[0] = 100
parameters = constants.TEST_PARAMETERS

def plot_with_initial(initial_d, initial_e):
    print(f"Dosage: {initial_d} CSAN, {initial_e} T Cells")
    state.D = initial_d
    state.Es[0] = initial_e
    full_result = model.generate_simulation_data(parameters, state, [duration * (i / timesteps) for i in range(timesteps + 1)], sample_count=10)["data"]
    plot_all_aggregate_normalized(full_result)
#make_average_free_csans_plot(full_result)




plot_with_initial(0, 0)

plot_with_initial(0, 100)

for i in range(4):
    plot_with_initial(4 ** i, 100)
    
