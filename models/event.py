from copy import deepcopy
from random import random
from abc import abstractmethod

import numpy as np
from scipy.linalg import solve
from models.model import Model

class Event:
    """
    Abstract class used by Event Models. 

    To instantiate an event, you must override:
     - rate function | (state, time, parameters) -> nonnegative number
     - max rate function | (state, parameters) -> nonnegative number
     - implement function | state -> state
    """

    @abstractmethod
    def get_rate(self, state, time, model_parameters):
        """Returns the instantaneous rate at which the event occurs at the current state and time."""

    @abstractmethod
    def get_max_rate(self, state, model_parameters):
        """Returns the maximum possible rate (across time) at which the event could occur."""

    @abstractmethod
    def implement(self, state,**kwargs):
        """Returns the new state after event is implemented.
        Mutates the state to become the new state."""

class TimeIndependentEvent(Event):
    """Class for events whose rate does not depend on time."""
    def get_rate(self, state, time, model_parameters):
        return self.get_max_rate(state, model_parameters)

class ConstantEvent(TimeIndependentEvent):
    """Constant events occur only when the model is in a certain state. 
    They transfor the state into a new state."""
    def __init__(self, initial_state, end_state, get_rate_from_parameters):
        self.initial_state = initial_state
        self.end_state = end_state
        self.get_rate_from_parameters = get_rate_from_parameters

    def get_max_rate(self, state, model_parameters):
        if state == self.initial_state:
            return self.get_rate_from_parameters(model_parameters)
        return 0

    def implement(self, state):
        return self.end_state

class EventModel(Model):
    """
    Abstract class to handle event-driven time-independent models.
    Parameters not included in instantiation.

    To understand a model you should be comfortable with some concepts: 
        - A state space. This is the type of data that the model will manipulate.
        - A set of possible events. These events occur stochastically and modify the state.
        Event types can be any subclasses of RateBoundedEvent.
        - Parameters. A single model may behave in different ways under different parameters.
        This is generally accomplished via having the parameters affect the event rates. 


    In order to instantiate a model, we must specify its events. 
    The state space is implicitly determined. Attempting to run models on states 
    outside its state space will generally lead to errors. 
    The parameters must be provided at runtime.

    In order to run a model, we must specify:
        - parameters
        - initial state
        - duration
    """

    name = "Event-Driven Model"

    def __init__(self, events: list[Event]):
        """"""
        self.events = events

    def run(self, parameters: dict, initial_state, duration: float, max_num_steps=None, verbose=False):
        """
        Returns the result of running the model. Does not mutate any of the arguments.  
        """
        current_state = deepcopy(initial_state)
        current_time = 0
        num_steps = 0
        while True:
            num_steps += 1
            if max_num_steps is not None and num_steps > max_num_steps:
                raise RuntimeError(
                    "Maximum number of steps for single simulation exceeded")
            rates = []
            for event in self.events:
                rates.append(event.get_max_rate(current_state, parameters))
            total_rate = sum(rates)
            if total_rate == 0:
                break
            else:
                random_number = random()
                waiting_time = - np.log(random_number) / total_rate
            current_time += waiting_time
            
            if current_time > duration:
                break
            event_index = random() * total_rate
            found_event = None
            for event, rate in zip(self.events, rates):
                if rate >= event_index:
                    found_event = event
                    found_max_rate = rate
                    found_rate = event.get_rate(current_state, current_time, parameters)
                    break
                event_index -= rate
            if found_event is None:
                raise RuntimeError("Event was not able to be found!")
            if found_rate / found_max_rate > random():
                found_event.implement(current_state)
        
        if verbose:
            print(current_state)
        return current_state

    


class ConstantEventModel(EventModel):
    def __init__(self, events: list[ConstantEvent]):
        self.state_space = []
        for event in events:
            for state in [event.initial_state, event.end_state]:
                if state not in self.state_space:
                    self.state_space.append(state)

        super().__init__(events)

    def get_stable_distribution(self, parameters):
        """
        Returns the left eigenvector of the transition matrix whose eigenvalue is 0 with appropriate norm.
        Accomplishes this by solving the equation Ax=b, where b = [1, 0, .., 0] and A is the matrix
        that takes the transpose of the transition matrix and replaces the first row with [1, 0, ..., 0].

        Will return a Singular Matrix Error if the space of these eigenvectors is 2+ dimensional,
        which corresponds to having multiple stable states (only possible for reducible markov processes). 
        """

        matrix_to_solve = self._get_transition_matrix(parameters).T
        row = np.zeros(len(self.state_space))
        row[0] = 1
        matrix_to_solve[0] = row
        unscaled_stable_vector = solve(matrix_to_solve, row)
        scale = sum(list(unscaled_stable_vector))
        stable_vector = list(unscaled_stable_vector / scale)

        stable_probabilities = {}
        for state, probability in zip(self.state_space, stable_vector):
            stable_probabilities[state] = probability
        return stable_probabilities

    def _get_transition_matrix(self, parameters):
        state_count = len(self.state_space)
        transition_matrix = np.zeros([state_count, state_count])
        for event in self.events:
            initial_index = self._get_list_index(
                self.state_space, event.initial_state)
            end_index = self._get_list_index(self.state_space, event.end_state)
            rate = event.get_max_rate(event.initial_state, parameters)
            transition_matrix[initial_index,
                              end_index] = transition_matrix[initial_index, end_index] + rate
            transition_matrix[initial_index,
                              initial_index] = transition_matrix[initial_index, initial_index] - rate
        return transition_matrix

    @staticmethod
    def _get_list_index(items, item):
        for i, entry in enumerate(items):
            if entry == item:
                return i





