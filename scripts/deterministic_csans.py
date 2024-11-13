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
        assert isinstance(bound_receptors, int)
        assert bound_receptors >= 0
        assert bound_receptors <= self.e_receptors

        return bound_receptors + 1
    
    def _get_t_index(self, bound_receptors):
        assert isinstance(bound_receptors, int)
        assert bound_receptors >= 0
        assert bound_receptors <= self.t_receptors
        
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

    def _single_success_prob(self, e_bound, t_bound, parameters, e_csan, t_csan):
        if e_csan and t_csan:
            return 0
        if e_csan: 
            return e_bound / self.e_receptors * (1 - t_bound / self.t_receptors) * parameters["p_ED|T"] 
        if t_csan:
            return (1 - e_bound / self.e_receptors) * t_bound / self.t_receptors * parameters["p_E|DT"]
        return (1 - e_bound / self.e_receptors) * (1 - t_bound / self.t_receptors) * parameters["p_E|T"]
        
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

        # ET killings
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            index = self._get_et_index(i, j)
            csans_change += parameters["d_ET"] * state[index] * j

        ## EDT killings
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            index = self._get_edt_index(i, j)
            csans_change += parameters["d_EDT"] * state[index] * j

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
                    
            for j in range(self.t_receptors):
                if i == 0:
                    continue
                index = self._get_edt_index(i - 1, j)
                e_change += state[index] * parameters["d_EDT"]

            # incoming csan bindings
            if i != 0:
                index = self._get_e_index(i - 1)
                e_change += state[self._get_d_index()] * state[index] * (self.e_receptors - i + 1) * parameters["lambda_E"]
                

            # outgoing csan bindings
            if i != self.e_receptors:
                index = self._get_e_index(i)
                e_change -= state[self._get_d_index()] * state[index] * (self.e_receptors - i) * parameters["lambda_E"]

            # incoming csan absorbtions
            pass
            # outgoing csan absorbtions
            pass
            # dimer/trimer formation
            for j in range(self.t_receptors + 1):
                t_index = self._get_t_index(j)
                e_index = self._get_e_index(i)
                
                e_bound_pro = i / self.e_receptors
                t_bound_pro = j / self.t_receptors
                success_prob = (1 - e_bound_pro) * (1 - t_bound_pro) * parameters["p_E|T"] + \
                    (e_bound_pro) * (1 - t_bound_pro) * parameters["p_ED|T"] + \
                    (1 - e_bound_pro) * (t_bound_pro) * parameters["p_E|DT"]
                total_prob = 1 - (1 - success_prob) ** parameters["M"]
                
                e_change -= state[t_index] * state[e_index] * parameters["lambda_ET"] * total_prob 

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
            for j in range(self.e_receptors + 1):
                t_index = self._get_t_index(i)
                e_index = self._get_e_index(j)
                
                e_bound_pro = i / self.e_receptors
                t_bound_pro = j / self.t_receptors
                success_prob = (1 - e_bound_pro) * (1 - t_bound_pro) * parameters["p_E|T"] + \
                    (e_bound_pro) * (1 - t_bound_pro) * parameters["p_ED|T"] + \
                    (1 - e_bound_pro) * (t_bound_pro) * parameters["p_E|DT"]
                total_prob = 1 - (1 - success_prob) ** parameters["M"]
                
                t_change -= state[t_index] * state[e_index] * parameters["lambda_ET"] * total_prob 

            index = self._get_t_index(i)
            derivative[index] = t_change

        # get ETs derivatives
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            et_change = 0
            
            # E|T formation
            e_index = self._get_e_index(i)
            t_index = self._get_t_index(j)
            single_ET = self._single_success_prob(i, j, parameters, False, False)
            single_ED_T = self._single_success_prob(i, j, parameters, True, False)
            single_E_DT = self._single_success_prob(i, j, parameters, False, True)
            single = single_ET + single_E_DT + single_ED_T
            if single_ET != 0:
                et_change += parameters["lambda_ET"] * (1 - (1 - single) ** parameters["M"]) * state[e_index] * state[t_index] * single_ET / single

            # ET killing
            index = self._get_et_index(i, j)
            et_change -= state[index] * parameters["d_ET"]
            

            index = self._get_et_index(i, j)
            derivative[index] = et_change


        # get EDTs derivatives 
        for i, j in product(range(self.e_receptors), range(self.t_receptors)):
            edt_change = 0
            # ED|T formation
            e_index = self._get_e_index(i + 1)
            t_index = self._get_t_index(j)
            single_ET = self._single_success_prob(i + 1, j, parameters, False, False)
            single_ED_T = self._single_success_prob(i + 1, j, parameters, True, False)
            single_E_DT = self._single_success_prob(i + 1, j, parameters, False, True)
            single = single_ET + single_E_DT + single_ED_T
            if single_ED_T != 0:
                edt_change += parameters["lambda_ET"] * (1 - (1 - single) ** parameters["M"]) * state[e_index] * state[t_index] * single_ED_T / single
            
            # E|DT formation
            e_index = self._get_e_index(i)
            t_index = self._get_t_index(j + 1)
            single_ET = self._single_success_prob(i, j + 1, parameters, False, False)
            single_ED_T = self._single_success_prob(i, j + 1, parameters, True, False)
            single_E_DT = self._single_success_prob(i, j + 1, parameters, False, True)
            single = single_ET + single_E_DT + single_ED_T
            if single_E_DT != 0:
                edt_change += parameters["lambda_ET"] * (1 - (1 - single) ** parameters["M"]) * state[e_index] * state[t_index] * single_E_DT / single

            index = self._get_edt_index(i, j)
            derivative[index] = edt_change
            
            # EDT killing
            index = self._get_edt_index(i, j)
            edt_change -= state[index] * parameters["d_EDT"]
            

            index = self._get_edt_index(i, j)
            derivative[index] = edt_change


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

    def _d(self, state):
        return state[0]
    def _e(self, state):
        return [state[self._get_e_index(i)] for i in range(self.e_receptors + 1)]
    def _t(self, state):
        return [state[self._get_t_index(i)] for i in range(self.t_receptors + 1)]
    def _et(self, state):
        return [[state[self._get_et_index(i, j)] for i in range(self.e_receptors + 1)] for j in range(self.t_receptors + 1)]
    def _edt(self, state):
        return [[state[self._get_edt_index(i, j)] for i in range(self.e_receptors + 1)] for j in range(self.t_receptors + 1)]
