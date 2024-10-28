from deterministic_csans import CSANDetModel
import matplotlib.pyplot as plt
from copy import copy
from constants import ZEROS, EFFECTIVE_PARAMS as PARAMS

model = CSANDetModel(20, 20)
initial_state = model.get_initial_state(0, 30000, 10000)

#parameters = copy(ZEROS)
#parameters["b_E"] = 0.1
#parameters["b_T"] = 0.1

parameters = PARAMS

timepoint_count = 100
max_time = 48

timepoints = [max_time / (timepoint_count - 1) * i for i in range(timepoint_count)]

initial_doses = [0] + list([640000 / (4 ** i) for i in range(6)])
print(initial_doses)

fig, ax = plt.subplots()

for dose in initial_doses:
    initial_state = model.get_initial_state(dose, 30000, 10000)
    res = model.run(parameters, initial_state, timepoints)
    d_res = [float(model._d(timepoint)) for timepoint in res]
    e_res = [float(sum(model._e(timepoint))) for timepoint in res]
    t_res = [float(sum(model._t(timepoint))) for timepoint in res]
    print(d_res)
    print(e_res)
    print(t_res)

    #ax.plot(timepoints, d_res, label = "CSAN count")
    #ax.plot(timepoints, e_res, label = "Effector cell Count")
    ax.plot(timepoints, t_res, label = f"DOSE = {dose}")

ax.legend()
ax.set_xlim(0, 50)
ax.set_ylim(0, 50000)
plt.show()