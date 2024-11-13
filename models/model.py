"""
Models are implemented in a parameter-agnostic way. 
Parameters instead are supplied at the time of running.
"""
# pylint:disable=arguments-differ
from abc import abstractmethod

class Model:
    """Abstract class. Contains a state space and function to run for a duration."""
    name = "Abstract Model"

    @abstractmethod
    def run(self, parameters, initial_state, duration: float, **kwargs):
        """Returns the result of running the model on initial_state for a duration with given parameters."""

    def generate_simulation_data(self, parameters: dict, initial_state, timepoints: list, sample_count: int = 1,
                                 **kwargs):
        """
        Returns result of run between timepoints starting from initial_state sample_count times.
        
        Data is output in json-style:
        {
            "model": [model name],
            "parameters": [parameter dictionary],
            "data": [timepoint data dictionary],
        }
        """


        simulation_result = {
            "parameters": parameters,
            "model": self.name,
            "data": [],
            "timepoints": timepoints
        }

        percent_completed = -10


        for i in range(sample_count):
            if sample_count > 100 and i * 100 / sample_count >= percent_completed + 10:
                percent_completed += 10
                print(f"{i}/{sample_count} completed")

            timepoint_data = {0: initial_state}
            last_time = 0
            current_state = initial_state
            for time in timepoints:
                duration = time - last_time
                current_state = self.run(
                    parameters, current_state, duration, **kwargs)
                timepoint_data[time] = current_state
                last_time = time
            simulation_result["data"].append(timepoint_data)
        return simulation_result

