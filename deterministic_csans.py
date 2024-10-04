from model import Model
from itertools import product
from scipy.special import binom

class CSANDetModel(Model):
    
    def __init__(self, e_receptors, t_receptors):
        self.e_receptors = e_receptors
        self.t_receptors = t_receptors
        self.dimension = 3 + e_receptors + t_receptors + 2 * e_receptors * t_receptors

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

    def get_derivative(self, state, parameters):
        derivative = [0] * self.dimension

        # get CSANs derivative
        csans_change = 0
        
        # bindings to effectors
        for i in range(self.e_receptors + 1):
            index = self._get_e_index(i)
            csans_change -= parameters["lambda_E"] * state[index] * (self.e_receptors - i)
        
        # bindings to tumors
        for i in range(self.t_receptors + 1):
            index = self._get_t_index(i)
            csans_change -= parameters["lambda_T"] * state[index] * (self.t_receptors - i)
        
        # tumor cells dying
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            index = self._get_et_index(i, j)
            csans_change += parameters["d_ET"] * state[index] * j
            index = self._get_edt_index(i, j)
            csans_change += parameters["d_EDT"] * state[index] * j
        
        derivative[0] = csans_change


        # get Es derivatives
        for i in range(self.e_receptors + 1):
            e_change = 0
            
            # incoming births
            for j in range(self.e_receptors + 1):
                index = self._get_e_index(j)
                e_change += parameters["b_E"] * binom(j, i) / (2 ** j) * state[index]
            index = self._get_e_index(i)
            
            # outgoing births
            e_change -= parameters["b_E"] * state[index]
            
            # deaths
            e_change -= parameters["d_E"] * state[index]

            # incoming dimer/trimer kills
            for j in range(self.t_receptors):
                index = self._get_et_index(i, j)
                e_change += state[index] * parameters["d_ET"]
                if j == 0:
                    continue
                index = self._get_edt_index(i, j - 1)
                e_change += state[index] * parameters["d_ET"]

            # incoming csan bindings
            


            # outgoing csan bindings

            # incoming csan absorbtions

            # outgoing csan absorbtions

            # dimer/trimer formation
            pass
        
        # get Ts derivatives
        for i in range(self.t_receptors + 1):
            # births
            # deaths
            # csan bindings
            # absorbtions
            # dimer/trimer formation
            pass
        
        # get ETs derivatives
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            pass

        # get EDTs derivatives 
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            pass
