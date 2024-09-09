import matplotlib.pyplot as plt


def make_average_free_csans_plot(result):
    csans_count = {key : 0 for key in result[0].keys()}
    for sample_path in result:
        for key, val in sample_path.items():
            csans_count[key] = csans_count[key] + val.D

    sample_count = len(result)
    csans_average = {}
    for key in result[0].keys():
        csans_average[key] = csans_count[key] / sample_count
    
    fig, ax = plt.subplots()

    xs = list(csans_average.keys())
    ys = [csans_average[key] for key in xs]

    ax.plot(xs, ys)
    ax.set_xlabel("Time")
    ax.set_ylabel("Free CSANs")
    plt.show()



def make_average_TCell_plot(result, times_to_track):
    sample_count = len(result)
    times = list(result[0].keys())
    e_receptors = len(result[0][times[0]].Es) - 1
    t_receptors = len(result[0][times[0]].Ts) - 1
    xs = list(range(e_receptors + 1))
    fig, ax = plt.subplots()
    for time in times_to_track:
        total_TCells = [0 for _ in range(e_receptors + 1)]
        for sample in result:
            TCell_timepoint = sample[time].Es
            for i, TCell_count in enumerate(TCell_timepoint):
                total_TCells[i] = total_TCells[i] + TCell_count
        ys = [total / sample_count for total in total_TCells]
        ax.plot(xs, ys, label=f"t = {time}")
    

    ax.set_xlabel("Number of Bound Receptors")
    ax.set_ylabel("Number of T Cells")
    plt.legend()
    plt.show()

    



def make_average_tumor_cell_plot(result, times_to_track):
    sample_count = len(result)
    times = list(result[0].keys())
    e_receptors = len(result[0][times[0]].Es) - 1
    t_receptors = len(result[0][times[0]].Ts) - 1
    xs = list(range(t_receptors + 1))
    fig, ax = plt.subplots()
    for time in times_to_track:
        total_tumor_cells = [0 for _ in range(t_receptors + 1)]
        for sample in result:
            tumor_cell_timepoint = sample[time].Ts
            for i, tumor_cell_count in enumerate(tumor_cell_timepoint):
                total_tumor_cells[i] = total_tumor_cells[i] + tumor_cell_count
        ys = [total / sample_count for total in total_tumor_cells]
        ax.plot(xs, ys, label=f"t = {time}")
    

    ax.set_xlabel("Number of Bound Receptors")
    ax.set_ylabel("Number of Tumor Cells")
    plt.legend()
    plt.show()

def make_total_tumor_cell_plot(result):

    

    sample_count = len(result)
    times = list(result[0].keys())
    e_receptors = len(result[0][times[0]].Es) - 1
    t_receptors = len(result[0][times[0]].Ts) - 1
    xs = times
    fig, ax = plt.subplots()
    ys = []
    for time in times:
        total_tumor_cells = 0
        for sample in result:
            total_tumor_cells += sum(sample[time].Ts)    
        ys.append(total_tumor_cells / sample_count)
    ax.plot(xs, ys)
    

    ax.set_xlabel("Time")
    ax.set_ylabel("Number of Tumor Cells")
    initial_count = sum(result[0][times[0]].Ts)
    ax.set_ybound(0, initial_count * 2)
    plt.show()


def plot_all_aggregate_normalized(result):
    sample_count = len(result)
    times = list(result[0].keys())
    e_receptors = len(result[0][times[0]].Es) - 1
    t_receptors = len(result[0][times[0]].Ts) - 1
    fig, ax = plt.subplots()
    initial_state = result[0][0]
    d = {time: 0 for time in times}
    t = {time: 0 for time in times}
    e = {time: 0 for time in times}
    et = {time: 0 for time in times}
    edt = {time: 0 for time in times}
    d_total = {time: 0 for time in times}
    bound_d_per_e_total = {time: 0 for time in times}
    bound_d_per_t_total = {time: 0 for time in times}

    for sample in result:
        for time in times:
            state = sample[time]
            d[time] = d[time] + state.D
            t[time] = t[time] + sum(state.Ts)
            e[time] = e[time] + sum(state.Es)
            et[time] = et[time] + sum([sum(row) for row in state.ETs])
            edt[time] = edt[time] + sum([sum(row) for row in state.EDTs])
            d_total[time] = d_total[time] + state.total_D()
            bound_d_per_e_total[time] = bound_d_per_e_total[time] + state.bound_D_per_E()
            bound_d_per_t_total[time] = bound_d_per_t_total[time] + state.bound_D_per_T()

    dscalar = max(initial_state.D, 1)
    tscalar = max(initial_state.Ts[0], 1)
    escalar = max(initial_state.Es[0], 1)
    y_d = [d[time] / sample_count / dscalar for time in times]
    y_t = [t[time] / sample_count / tscalar for time in times]
    y_e = [e[time] / sample_count / escalar for time in times]
    y_et = [et[time] / sample_count / tscalar for time in times]
    y_edt = [edt[time] / sample_count / tscalar for time in times]
    y_d_total = [d_total[time] / sample_count / dscalar for time in times]
    y_bound_d_per_e = [bound_d_per_e_total[time] / sample_count / e_receptors for time in times]
    y_bound_d_per_t = [bound_d_per_t_total[time] / sample_count / t_receptors for time in times]

    x = times
    ax.plot(x, y_d, label="Normalized Free CSANs")
    ax.plot(x, y_t, label="Normalized Free Tumor Cells")
    ax.plot(x, y_e, label="Normalized Free T Cells")
    ax.plot(x, y_et, label="Normalized T Cell - Tumor Bindings")
    ax.plot(x, y_edt, label="Normalized T Cell - CSAN - Tumor Bindings")
    ax.plot(x, y_d_total, label = "Normalized CSANs (All)")
    ax.plot(x, y_bound_d_per_e, label = "Fraction of Free T Cell Receptors Bound to CSAN")
    ax.plot(x, y_bound_d_per_t, label = "Fraction of Free Tumor Receptors Bound to CSAN")
    ax.set_ybound(0, 1.5)
    ax.set_xlabel("Time")
    ax.set_ylabel("Normalized Entity Counts")
    plt.legend()
    plt.show()
    
    

