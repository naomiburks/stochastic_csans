from stochastic_csans import CSANModel
import matplotlib.pyplot as plt
import constants




duration = 48
timesteps = 480
times_to_plot = [duration * i * 120 / timesteps for i in range(5)]


e_receptors = 8
t_receptors = 8
sample_count = 100
model = CSANModel(e_receptors, t_receptors)
state = model.get_empty_state()
state.D = 10000
state.Es[0] = 100
state.Ts[0] = 100
parameters = constants.TEST_PARAMETERS
xs = list(range(e_receptors + 1))
def plot_from_initial(initial_d, initial_e):
    print(f"Dosage: {initial_d} CSAN, {initial_e} T Cells")
    state.D = initial_d
    state.Es[0] = initial_e
    full_result = model.generate_simulation_data(parameters, state, [duration * (i / timesteps) for i in range(timesteps + 1)], sample_count=sample_count)["data"]
    fig, ax = plt.subplots()
    bar_count = len(times_to_plot)
    for j, time in enumerate(times_to_plot):
        e_sum = [0 for _ in range(e_receptors + 1)]
        e_squares = [0 for _ in range(e_receptors + 1)]
        for sample in full_result:
            e_counts = sample[time].Es
            e_cell_total = sum(e_counts)
            for i, count in enumerate(e_counts):
                e_sum[i] = e_sum[i] + count / e_cell_total
                e_squares[i] = e_squares[i] + (count / e_cell_total) ** 2
        e_ave = [sum / sample_count for sum in e_sum]
        e_ave_sq = [sq / sample_count for sq in e_squares]
        e_standard_dev = [(e_ave_sq[i] - e_ave[i] ** 2) ** 0.5 for i in range(e_receptors + 1)]
        """ax.errorbar(xs, e_ave, e_standard_dev, label = f"t = {time}", capsize=4)"""
        true_xs = [x - 0.32 + 0.64 * j / (bar_count - 1) for x in xs]
        ax.bar(true_xs, e_ave, width=0.16, label = f"t = {time}")
 
    ax.set_ylabel("Fraction of T Cells")
    ax.set_xlabel("Number of Receptors Bound")
    ax.set_xticks(list(range(e_receptors + 1)))
    ax.set_ybound(0, 1)
    ax.legend()
    plt.savefig(f"plots/receptor_distribution_tcell_{initial_d}_{initial_e}")

    fig, ax = plt.subplots()
    
    for j, time in enumerate(times_to_plot):
        extinction_count = 0
        t_sum = [0 for _ in range(t_receptors + 1)]
        t_squares = [0 for _ in range(t_receptors + 1)]
        for sample in full_result:
            t_counts = sample[time].Ts
            t_cell_total = sum(t_counts)
            if t_cell_total == 0:
                extinction_count += 1
                continue
            for i, count in enumerate(t_counts):
                t_sum[i] = t_sum[i] + count / t_cell_total
                t_squares[i] = t_squares[i] + (count / t_cell_total) ** 2
        
        
        if sample_count == extinction_count:
            continue
        t_ave = [sum / (sample_count - extinction_count) for sum in t_sum]
        t_ave_sq = [sq / (sample_count - extinction_count) for sq in t_squares]
        t_standard_dev = [(t_ave_sq[i] - t_ave[i] ** 2) ** 0.5 for i in range(t_receptors + 1)]
        """ax.errorbar(xs, t_ave, t_standard_dev, label = f"t = {time}", capsize=4)"""
        true_xs = [x - 0.32 + 0.64 * j / (bar_count - 1) for x in xs]
        ax.bar(true_xs, t_ave, width=0.16, label = f"t = {time}")
    ax.set_ylabel("Fraction of Tumor Cells")
    ax.set_xlabel("Number of Receptors Bound")
    ax.set_ybound(0, 1)
    ax.legend()
    plt.savefig(f"plots/tumor_distribution_tcell_{initial_d}_{initial_e}")



#make_average_free_csans_plot(full_result)



"""
plot_from_initial(0, 0)

plot_from_initial(0, 100)
"""
for i in range(5):
    plot_from_initial(4 ** (i + 2), 100)
    
