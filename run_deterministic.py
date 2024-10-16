from deterministic_csans import CSANDetModel
from constants import NO_TRIMER_PARAMS



m = CSANDetModel(1, 1)
params = NO_TRIMER_PARAMS
times = [i * 4 for i in range(13)]
initial_state = m.get_initial(10, 10, 10)

res = m.run(params, initial_state, times)
print(res)
for data in res:
    m.show(data)
