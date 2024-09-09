from event import TimeIndependentEvent as TIE, EventModel
import numpy as np

class State:
    def __init__(self, D, Es, Ts, EDTs, ETs):
        self.D = D
        self.Es = Es
        self.Ts = Ts
        self.EDTs = EDTs
        self.ETs = ETs

    def __repr__(self):
        return f"D: {self.D}. Es: {self.Es}. Ts: {self.Ts}. EDTs: {self.EDTs}. ETs: {self.ETs}"

    def verify_coherence(self):
        if self.D < 0:
            return False
        for E in self.Es:
            if E < 0:
                return False
        for T in self.Ts:
            if T < 0:
                return False
        for EDT_by_E in self.EDTs:
            for EDT in EDT_by_E:
                if EDT < 0:
                    return False
        for ET_by_E in self.ETs:
            for ET in ET_by_E:
                if ET < 0:
                    return False
        return True

    def E(self):
        return sum(self.Es)
    
    def T(self):
        return sum(self.Es)
    
    def EDT(self):
        return sum([sum(row) for row in self.EDTs])
    
    def ET(self):
        return sum([sum(row) for row in self.ETs])

    def total_D(self):
        total = self.D
        for i, E in enumerate(self.Es):
            total += i * E
        for i, T in enumerate(self.Ts):
            total += i * T
        for i, row in enumerate(self.EDTs):
            for j, EDT in enumerate(row):
                total += (i + j + 1) * EDT
        for i, row in enumerate(self.ETs):
            for j, ET in enumerate(row):
                total += (i + j) * ET
        return total

    def bound_D_per_E(self):
        total_bound_D = 0
        total_E = 0
        for i, E in enumerate(self.Es):
            total_bound_D += i * E
            total_E += E
        if total_bound_D == 0:
            return total_bound_D
        return total_bound_D / total_E
    
    def bound_D_per_T(self):
        total_bound_D = 0
        total_T = 0
        for i, T in enumerate(self.Ts):
            total_bound_D += i * T
            total_T += T
        if total_bound_D == 0:
            return total_bound_D
        return total_bound_D / total_T


class CSANModel(EventModel):
    def __init__(self, e_receptors, t_receptors):
        self.e_receptors = e_receptors
        self.t_receptors = t_receptors

        events = []
        for bound_CSAN_count in range(e_receptors + 1):
            events.append(TCellDeath(bound_CSAN_count))
            events.append(TCellBinding(bound_CSAN_count, e_receptors))
            events.append(TCellInternalization(bound_CSAN_count, e_receptors))
        for bound_CSAN_count in range(t_receptors + 1):
            events.append(TumorBirth(bound_CSAN_count))
            events.append(TumorBinding(bound_CSAN_count, t_receptors))
            events.append(TumorInternalization(bound_CSAN_count, t_receptors))

        for TCell_bound_count in range(e_receptors + 1):
            for tumor_bound_count in range(t_receptors + 1):
                events.append(ED_TFormation(TCell_bound_count, tumor_bound_count, e_receptors, t_receptors))
                events.append(E_DTFormation(TCell_bound_count, tumor_bound_count, e_receptors, t_receptors))
                events.append(E_TFormation(TCell_bound_count, tumor_bound_count, e_receptors, t_receptors))
                
        for TCell_bound_count in range(e_receptors):
            for tumor_bound_count in range(t_receptors):
                events.append(TrimerDeath(TCell_bound_count, tumor_bound_count, e_receptors, t_receptors)) 
                events.append(DimerDeath(TCell_bound_count, tumor_bound_count, e_receptors, t_receptors))
        self.events = events
    
    
    def get_empty_state(self):
        D = 0
        Es = [0 for _ in range(self.e_receptors + 1)]
        Ts = [0 for _ in range(self.t_receptors + 1)]
        EDTs = [[0 for _ in range(self.t_receptors)] for __ in range(self.e_receptors)]
        ETs = [[0 for _ in range(self.t_receptors)] for __ in range(self.e_receptors)]
        
        return State(D, Es, Ts, EDTs, ETs)
    
class TCellDeath(TIE):
    def __init__(self, bound_CSAN_count):
        self.bound_CSAN_count = bound_CSAN_count

    def get_max_rate(self, state, parameters):
        return state.Es[self.bound_CSAN_count] * parameters["d_E"]

    def implement(self, state):
        state.Es[self.bound_CSAN_count] = state.Es[self.bound_CSAN_count] - 1
        state.D = state.D + self.bound_CSAN_count

class TumorBirth(TIE):
    def __init__(self, bound_CSAN_count):
        self.bound_CSAN_count = bound_CSAN_count
    
    def get_max_rate(self, state, parameters):
        return state.Ts[self.bound_CSAN_count] * parameters["b_T"]
        

    def implement(self, state):
        state.Ts[self.bound_CSAN_count] = state.Ts[self.bound_CSAN_count] - 1
        child_csans = np.random.binomial(self.bound_CSAN_count, 0.5)
        child_csans2 = self.bound_CSAN_count - child_csans
        state.Ts[child_csans] = state.Ts[child_csans] + 1
        state.Ts[child_csans2] = state.Ts[child_csans2] + 1
        return state
 
class TCellBinding(TIE):
    def __init__(self, bound_CSAN_count, receptor_count):
        self.bound_CSAN_count = bound_CSAN_count
        self.receptor_count = receptor_count
        self.free_receptor_count = self.receptor_count - self.bound_CSAN_count
    
    def get_max_rate(self, state, parameters):
        return state.D * self.free_receptor_count * state.Es[self.bound_CSAN_count] * parameters["lambda_E"]
        
    
    def implement(self, state):
        state.D = state.D - 1
        state.Es[self.bound_CSAN_count] = state.Es[self.bound_CSAN_count] - 1
        state.Es[self.bound_CSAN_count + 1] = state.Es[self.bound_CSAN_count + 1] + 1
        return state

class TumorBinding(TIE):
    def __init__(self, bound_CSAN_count, receptor_count):
        self.bound_CSAN_count = bound_CSAN_count
        self.receptor_count = receptor_count
        self.free_receptor_count = self.receptor_count - self.bound_CSAN_count
    
    def get_max_rate(self, state, parameters):
        return state.D * self.free_receptor_count * state.Ts[self.bound_CSAN_count] * parameters["lambda_T"]
         
    def implement(self, state):
        state.D = state.D - 1
        state.Ts[self.bound_CSAN_count] = state.Ts[self.bound_CSAN_count] - 1
        state.Ts[self.bound_CSAN_count + 1] = state.Ts[self.bound_CSAN_count + 1] + 1
        return state

class TCellInternalization(TIE):
    def __init__(self, bound_CSAN_count, receptor_count):
        self.bound_CSAN_count = bound_CSAN_count
        self.receptor_count = receptor_count
        self.free_receptor_count = self.receptor_count - self.bound_CSAN_count
    
    def get_max_rate(self, state, parameters):
        return self.bound_CSAN_count * state.Es[self.bound_CSAN_count] * parameters["mu_E"]
        
    
    def implement(self, state):
        state.Es[self.bound_CSAN_count] = state.Es[self.bound_CSAN_count] - 1
        state.Es[self.bound_CSAN_count - 1] = state.Es[self.bound_CSAN_count - 1] + 1
        return state

class TumorInternalization(TIE):
    def __init__(self, bound_CSAN_count, receptor_count):
        self.bound_CSAN_count = bound_CSAN_count
        self.receptor_count = receptor_count
        self.free_receptor_count = self.receptor_count - self.bound_CSAN_count
    
    def get_max_rate(self, state, parameters):
        return self.bound_CSAN_count * state.Ts[self.bound_CSAN_count] * parameters["mu_T"]
        
    
    def implement(self, state):
        state.Ts[self.bound_CSAN_count] = state.Ts[self.bound_CSAN_count] - 1
        state.Ts[self.bound_CSAN_count - 1] = state.Ts[self.bound_CSAN_count - 1] + 1
        return state


class ED_TFormation(TIE):
    def __init__(self, TCell_bound_count, tumor_bound_count, TCell_receptors, tumor_receptors):
        self.TCell_bound_count = TCell_bound_count
        self.tumor_bound_count = tumor_bound_count
        self.TCell_receptors = TCell_receptors
        self.tumor_receptors = tumor_receptors

    def get_max_rate(self, state, model_parameters):
        p_e_t = model_parameters["p_E|T"]
        p_ed_t = model_parameters["p_ED|T"]
        p_e_dt = model_parameters["p_E|DT"]
        lambda_ET = model_parameters["lambda_ET"]
        m = model_parameters["M"]
        c_e_t = (self.TCell_receptors - self.TCell_bound_count) * (self.tumor_receptors - self.tumor_bound_count) * p_e_t / (self.TCell_receptors * self.tumor_receptors)
        c_e_dt = (self.TCell_receptors - self.TCell_bound_count) * self.tumor_bound_count * p_e_dt / (self.TCell_receptors * self.tumor_receptors)
        c_ed_t = self.TCell_bound_count * (self.tumor_receptors - self.tumor_bound_count) * p_ed_t / (self.TCell_receptors * self.tumor_receptors)

        c_f = 1 - c_e_t - c_e_dt - c_ed_t
        if c_f == 1:
            return 0

        return lambda_ET * (1 - c_f ** m) * state.Es[self.TCell_bound_count] * state.Ts[self.tumor_bound_count] * c_ed_t / (c_e_t + c_e_dt + c_ed_t)

    def implement(self, state):
        state.Es[self.TCell_bound_count] = state.Es[self.TCell_bound_count] - 1
        state.Ts[self.tumor_bound_count] = state.Ts[self.tumor_bound_count] - 1
        state.EDTs[self.TCell_bound_count - 1][self.tumor_bound_count] = state.EDTs[self.TCell_bound_count - 1][self.tumor_bound_count] + 1

class E_DTFormation(TIE):
    def __init__(self, TCell_bound_count, tumor_bound_count, TCell_receptors, tumor_receptors):
        self.TCell_bound_count = TCell_bound_count
        self.tumor_bound_count = tumor_bound_count
        self.TCell_receptors = TCell_receptors
        self.tumor_receptors = tumor_receptors

    def get_max_rate(self, state, model_parameters):
        p_e_t = model_parameters["p_E|T"]
        p_ed_t = model_parameters["p_ED|T"]
        p_e_dt = model_parameters["p_E|DT"]
        lambda_ET = model_parameters["lambda_ET"]
        m = model_parameters["M"]
        c_e_t = (self.TCell_receptors - self.TCell_bound_count) * (self.tumor_receptors - self.tumor_bound_count) * p_e_t / (self.TCell_receptors * self.tumor_receptors)
        c_e_dt = (self.TCell_receptors - self.TCell_bound_count) * self.tumor_bound_count * p_e_dt / (self.TCell_receptors * self.tumor_receptors)
        c_ed_t = self.TCell_bound_count * (self.tumor_receptors - self.tumor_bound_count) * p_ed_t / (self.TCell_receptors * self.tumor_receptors)

        c_f = 1 - c_e_t - c_e_dt - c_ed_t

        
        if c_f == 1:
            return 0
        return lambda_ET * (1 - c_f ** m) * state.Es[self.TCell_bound_count] * state.Ts[self.tumor_bound_count] * c_e_dt / (c_e_t + c_e_dt + c_ed_t)

    def implement(self, state):
        state.Es[self.TCell_bound_count] = state.Es[self.TCell_bound_count] - 1
        state.Ts[self.tumor_bound_count] = state.Ts[self.tumor_bound_count] - 1
        state.EDTs[self.TCell_bound_count][self.tumor_bound_count - 1] = state.EDTs[self.TCell_bound_count][self.tumor_bound_count - 1] + 1

class E_TFormation(TIE):
    def __init__(self, TCell_bound_count, tumor_bound_count, TCell_receptors, tumor_receptors):
        self.TCell_bound_count = TCell_bound_count
        self.tumor_bound_count = tumor_bound_count
        self.TCell_receptors = TCell_receptors
        self.tumor_receptors = tumor_receptors

    def get_max_rate(self, state, model_parameters):
        p_e_t = model_parameters["p_E|T"]
        p_ed_t = model_parameters["p_ED|T"]
        p_e_dt = model_parameters["p_E|DT"]
        lambda_ET = model_parameters["lambda_ET"]
        m = model_parameters["M"]
        c_e_t = (self.TCell_receptors - self.TCell_bound_count) * (self.tumor_receptors - self.tumor_bound_count) * p_e_t / (self.TCell_receptors * self.tumor_receptors)
        c_e_dt = (self.TCell_receptors - self.TCell_bound_count) * self.tumor_bound_count * p_e_dt / (self.TCell_receptors * self.tumor_receptors)
        c_ed_t = self.TCell_bound_count * (self.tumor_receptors - self.tumor_bound_count) * p_ed_t / (self.TCell_receptors * self.tumor_receptors)

        c_f = 1 - c_e_t - c_e_dt - c_ed_t

        
        if c_f == 1:
            return 0
        return lambda_ET * (1 - c_f ** m) * state.Es[self.TCell_bound_count] * state.Ts[self.tumor_bound_count] * c_e_t / (c_e_t + c_e_dt + c_ed_t)

    def implement(self, state):
        state.Es[self.TCell_bound_count] = state.Es[self.TCell_bound_count] - 1
        state.Ts[self.tumor_bound_count] = state.Ts[self.tumor_bound_count] - 1
        state.ETs[self.TCell_bound_count][self.tumor_bound_count] = state.ETs[self.TCell_bound_count][self.tumor_bound_count] + 1


class TrimerDeath(TIE):
    def __init__(self, TCell_bound_count, tumor_bound_count, TCell_receptors, tumor_receptors):
        self.TCell_bound_count = TCell_bound_count
        self.tumor_bound_count = tumor_bound_count
        self.TCell_receptors = TCell_receptors
        self.tumor_receptors = tumor_receptors

    def get_max_rate(self, state, model_parameters):
        d_EDT = model_parameters["d_EDT"]
        return d_EDT * state.EDTs[self.TCell_bound_count][self.tumor_bound_count]
    
    def implement(self, state):
        state.EDTs[self.TCell_bound_count][self.tumor_bound_count] = state.EDTs[self.TCell_bound_count][self.tumor_bound_count] - 1
        state.D = state.D + self.tumor_bound_count
        state.Es[self.TCell_bound_count + 1] = state.Es[self.TCell_bound_count + 1] + 1 

class DimerDeath(TIE):
    def __init__(self, TCell_bound_count, tumor_bound_count, TCell_receptors, tumor_receptors):
        self.TCell_bound_count = TCell_bound_count
        self.tumor_bound_count = tumor_bound_count
        self.TCell_receptors = TCell_receptors
        self.tumor_receptors = tumor_receptors
    
    def get_max_rate(self, state, model_parameters):
        d_ET = model_parameters["d_ET"]
        return d_ET * state.ETs[self.TCell_bound_count][self.tumor_bound_count]
    
    def implement(self, state):
        state.ETs[self.TCell_bound_count][self.tumor_bound_count] = state.ETs[self.TCell_bound_count][self.tumor_bound_count] - 1
        state.D = state.D + self.tumor_bound_count
        state.Es[self.TCell_bound_count] = state.Es[self.TCell_bound_count] + 1 
