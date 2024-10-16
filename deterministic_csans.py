from model import Model
from itertools import product
from scipy.special import binom
from scipy.integrate import odeint
from constants import ZEROS
from copy import copy

class CSANDetModel(Model):
    
    def __init__(self, e_receptors, t_receptors):
        self.e_receptors = e_receptors
        self.t_receptors = t_receptors
        self.dimension = 3 + e_receptors + t_receptors + 2 * e_receptors * t_receptors

    def run(self, parameters, initial_state, timepoints):
        f = lambda state, t: self.derivative(state, parameters)
        result = odeint(f, initial_state, timepoints)
        return result


    def _get_d_index(self):
        return 0
    
    def _get_e_index(self, bound_receptors):
        return bound_receptors + 1
    
    def _get_t_index(self, bound_receptors):
        return self.e_receptors + 2 + bound_receptors

    def _get_et_index(self, e_bound, t_bound):
        return self.e_receptors + self.t_receptors + 3 + e_bound * self.t_receptors + t_bound
    
    def _get_edt_index(self, e_bound, t_bound):
        return self.e_receptors + self.t_receptors + 3 + self.e_receptors * self.t_receptors + e_bound * self.t_receptors + t_bound

    def d(self, state):
        return state[self._get_d_index()]
    
    def t(self, state):
        total = 0
        for i in range(self.t_receptors + 1):
            total += state[self._get_t_index(i)]
        return total 
    
    def e(self, state):
        total = 0
        for i in range(self.e_receptors + 1):
            total += state[self._get_e_index(i)]
        return total

    def get_derivative(self, state, parameters):
        derivative = [0] * self.dimension

        # get CSANs derivative
        csans_change = 0
        
        # bindings to effectors
        for i in range(self.e_receptors + 1):
            index = self._get_e_index(i)
            csans_change -= parameters["lambda_E"] * state[index] * (self.e_receptors - i) * state[self._get_d_index()]
        
        # bindings to tumors
        for i in range(self.t_receptors + 1):
            index = self._get_t_index(i)
            csans_change -= parameters["lambda_T"] * state[index] * (self.t_receptors - i) * state[self._get_d_index()]
        
        # tumor cells dying
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            index = self._get_et_index(i, j)
            csans_change += parameters["d_ET"] * state[index] * j
            index = self._get_edt_index(i, j)
            csans_change += parameters["d_EDT"] * state[index] * j
        
        # effector cells dying
        for i in range(self.e_receptors + 1):
            index = self._get_e_index(i)
            csans_change += parameters["d_E"] * state[index] * i

        derivative[0] = csans_change


        # get Es derivatives
        for i in range(self.e_receptors + 1):
            e_change = 0
            
            # incoming births
            for j in range(self.e_receptors + 1):
                index = self._get_e_index(j)
                e_change += 2 * parameters["b_E"] * binom(j, i) / (2 ** j) * state[index]
            
            # outgoing births
            index = self._get_e_index(i)
            e_change -= parameters["b_E"] * state[index]
            
            # deaths
            e_change -= parameters["d_E"] * state[index]

            # incoming dimer kills
            for j in range(self.t_receptors):
                if i == self.e_receptors:
                    continue
                index = self._get_et_index(i, j)
                e_change += state[index] * parameters["d_ET"]
                if i == 0:
                    continue
                index = self._get_edt_index(i - 1, j)
                e_change += state[index] * parameters["d_EDT"]

            # incoming csan bindings
            if i != 0:
                index = self._get_e_index(i - 1)
                e_change += state[index] * (self.e_receptors - i + 1) * parameters["lambda_E"]
                

            # outgoing csan bindings
            if i != self.e_receptors:
                index = self._get_e_index(i)
                e_change += state[index] * (self.e_receptors - i) * parameters["lambda_E"]

            # incoming csan absorbtions
            pass
            # outgoing csan absorbtions
            pass
            # dimer/trimer formation
            pass
            
            # put change in derivative
            index = self._get_e_index(i)
            derivative[index] = e_change


        # get Ts derivatives
        for i in range(self.t_receptors + 1):
            t_change = 0

            # incoming births
            for j in range(self.t_receptors + 1):
                index = self._get_t_index(j)
                t_change += 2 * parameters["b_T"] * binom(j, i) / (2 ** j) * state[index]
            
            # outgoing births
            index = self._get_t_index(i)
            t_change -= parameters["b_T"] * state[index]

            # deaths
            t_change -= parameters["d_T"] * state[index]

            # incoming csan bindings
            if i != 0:
                index = self._get_t_index(i - 1)
                t_change += state[index] * parameters["lambda_T"] * (self.t_receptors - i + 1) * state[self._get_d_index()]

            # outgoing csan bindings
            if i != self.t_receptors:
                index = self._get_t_index(i)
                t_change -= state[index] * parameters["lambda_T"] * (self.t_receptors - i) * state[self._get_d_index()]

            # incoming csan absorbtions
            pass

            # outgoing csan absorbtions
            pass

            # dimer/trimer formation
            pass

            index = self._get_t_index(i)
            derivative[index] = t_change

        # get ETs derivatives
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            # ET formation
            pass
            # ET killing
            pass

        # get EDTs derivatives 
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            # EDT formation
            pass

        return derivative
    
    def get_initial_state(self, d, e, t):
        state = [0] * (self.e_receptors + self.t_receptors + 2 * self.e_receptors * self.t_receptors + 3)
        state[self._get_d_index()] = d
        state[self._get_e_index(0)] = e
        state[self._get_t_index(0)] = t
        return state

    def run(self, parameters, initial_state, timepoints):
        f = lambda x, t: self.get_derivative(x, parameters)
        result = odeint(f, initial_state, timepoints)
        return result

if __name__ == "__main__":
    model = CSANDetModel(5, 5)
    initial_state = model.get_initial_state(10, 10, 10)
    parameters = copy(ZEROS)
    parameters["b_E"] = 1
    res = model.run(parameters, initial_state, [i / 10 for i in range(30)])
    print(res)