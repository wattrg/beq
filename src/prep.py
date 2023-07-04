#! /usr/bin/python3

import json
import os
import sys

from enum import Enum
import numpy as np

kB = 1.380649e-23

def maxwellian_distribution(flow_state, vel, gas, dv, n_vel_increments):
    """
    Compute the maxwellian velocity distribution for a given flow state

    Parameters
    ----------
    flow_state: FlowState
        The macroscopic flow quantities
    vel: float
        The velocity to evaluate the distribution function at
    mass: float
        The mass of a single gas particle
    dv: float
        The velocity increment
    n_vel_increments: int
        The number of velocity increments

    Returns
    -------
    list<float>
        The distribution function
    """
    mass = gas.mass
    mass_on_two_pi_kB_T = mass / (2 * np.pi * kB * flow_state.T)
    exponent = -mass_on_two_pi_kB_T * (vel - flow_state.v)**2
    return mass_on_two_pi_kB_T**(3.2) * np.exp(exponent)

class FlowState:
    __slots__ = ["rho", "T", "v"]

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            self.__setattr__(key, value)

class GasModel:
    _values = ["mass"]
    __slots__ = _values

    def __init__(self, species=None, **kwargs):
        if species:
            self._load_species(species)

        for key, value in kwargs.items():
            self.__setattr__(key, value)

    def _load_species(self, species):
        beq = os.environ["BEQ"]
        with open(f"{beq}/resources/species/{species}.json", "r") as species_data:
            data = json.load(species_data)

        for key, value in data.items():
            self.__setattr__(key, value)

    def to_dict(self):
        json_values = {}
        for value in self._values:
            json_values[value] = getattr(self, value)
        return json_values


class _JsonData:
    def __init__(self, **kwargs):
        beq = os.environ["BEQ"]
        with open(f"{beq}/resources/defaults/{self._default}", "r") as default_json:
            defaults = json.load(default_json)

        for key in self._json_values:
            self.__setattr__(key, defaults[key])

        for key, value in kwargs.items():
            self.__setattr__(key, value)

    def to_dict(self):
        json_values = {}
        for value in self._json_values:
            json_values[value] = getattr(self, value)
        return json_values

    def set_values(self, dict):
        for key, value in dict.items():
            self.__setattr__(key, value)



class Config(_JsonData):
    _json_values = [
        "title"
    ]

    _python_values = [
        "initial_condition",
        "solver",
        "domain",
        "equation",
        "gas_model",
    ]

    __slots__ = _json_values + _python_values

    _default = "config.json"

    def __init__(self, config=None, **kwargs):
        # start by filling out defaults
        super().__init__(**kwargs)
        
        # overwrite values with values from the config dictionary, if provided
        if config:
            for key, value in config.items():
                if key == "solver":
                    self.solver = make_solver(value)
                elif key == "domain":
                    self.domain = Domain(**value)
                elif key == "gas_model":
                    self.gas_model = GasModel(**value)
                elif key == "equation":
                    self.equation = make_equation(value)
                else:
                    self.__setattr__(key, value)

    def _write_json_config(self):
        if not os.path.exists("config"):
            os.mkdir("config")

        json_values = self.to_dict()
        json_values["solver"] = self.solver.to_dict()
        json_values["domain"] = self.domain.to_dict()
        json_values["equation"] = self.equation.to_dict()
        if hasattr(self, "gas_model"):
            json_values["gas_model"] = self.gas_model.to_dict()

        with open("config/config.json", "w") as file:
            json.dump(json_values, file, indent=4)

    def _write_direct_initial_condition(self):
        dx = self.domain.length / self.domain.number_cells
        with open("solution/phi_0.beq", "w") as ic:
            for i in range(self.domain.number_cells):
                x = (i + 0.5) * dx
                value = self.initial_condition(x)
                ic.write(f"{value}\n")

    def _write_equilibrium_initial_condition(self):
        dx = self.domain.length / self.domain.number_cells
        dv = self.equation.dv
        n_vel_inc = self.equation.n_vel_increments
        max_v = (n_vel_inc + 0.5) * dv
        vels = np.arange(0.5*dv, max_v, dv)
        vels = np.append(vels, -vels)
        with open("solution/phi_0.beq", "w") as ic:
            for vel in vels:
                for i in range(self.domain.number_cells):
                    x = (i + 0.5) * dx
                    flow_state = self.initial_condition(x)
                    dist = maxwellian_distribution(
                        flow_state, vel, self.gas_model, self.equation.dv, self.equation.n_vel_increments
                    )
                    ic.write(f"{dist}\n")

    def _write_initial_condition(self):
        if not os.path.exists("solution"):
            os.mkdir("solution")

        if self.equation.ic_type == ICType.EQUILIBRIUM:
            self._write_equilibrium_initial_condition()

        elif self.equation.ic_type == ICType.DIRECT:
            self._write_direct_initial_condition()

    def write(self):
        self._write_json_config()
        self._write_initial_condition()

class RungeKutta(_JsonData):
    _json_values = [
        "type",
        "number_stages",
        "max_step",
        "max_time",
        "print_frequency",
        "plot_every_n_steps",
        "plot_frequency",
        "cfl",
    ]
    __slots__ = _json_values
    _default = "runge_kutta.json"
    _type = "runge_kutta"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type = self._type

def make_solver(solver_config):
    """
    Make a solver from a dictionary representation of the json config

    Parameters
    ----------
    solver_config: dict
        The dictionary representation of the json config

    Returns
    -------
    Solver: A solver object
    """
    type = solver_config["type"]
    if type == "runge_kutta":
        return RungeKutta(**solver_config)
    else:
        raise Exception("Unknown solver type")

class ICType(Enum):
    DIRECT = 1
    EQUILIBRIUM = 2


class Domain(_JsonData):
    _json_values = [
        "number_cells",
        "length",
    ]
    __slots__ = _json_values
    _default = "domain.json"


class Advection(_JsonData):
    _json_values = [
        "type",
        "velocity"
    ]
    __slots__ = _json_values
    _default = "advection_eq.json"
    _type = "advection"
    ic_type = ICType.DIRECT

    def __init__(self, **kwargs):
        self.type = self._type
        super().__init__(**kwargs)
    
class Burgers(_JsonData):
    _json_values = [
        "type",
    ]

    __slots__ = _json_values
    _default = "burgers_eq.json"
    _type = "burgers"
    ic_type = ICType.DIRECT

    def __init__(self, **kwargs):
        self.type = self._type
        super().__init__(**kwargs)

class Boltzmann(_JsonData):
    _json_values = [
        "type",
        "dv",
        "n_vel_increments",
    ]
    __slots__ = _json_values
    _default = "boltzmann_eq.json"
    _type = "boltzmann"
    ic_type = ICType.EQUILIBRIUM

    def __init__(self, **kwargs):
        self.type = self._type
        super().__init__(**kwargs)

def make_equation(config):
    """
    Make an equation from the dictionary representation of the json config

    Parameters
    ----------
    config: dict
        The dictionary representation of the json config

    Returns
    -------
    equation:
        The equation
    """
    type = config.get("type", None)
    if type == "advection":
        return Advection(**config)
    elif type == "burgers":
        return Burgers(**config)
    elif type == "boltzmann":
        return Boltzmann(**config)
    else:
        raise Exception("Unknown equation")

def prep():
    prep_script = sys.argv[1]
    config = Config()
    namespace = {
        "config": config,
        "RungeKutta": RungeKutta,
        "Domain": Domain,
        "Advection": Advection,
        "Burgers": Burgers,
        "Boltzmann": Boltzmann,
        "FlowState": FlowState,
        "GasModel": GasModel,
    }
    with open(prep_script, "r") as f:
        exec(f.read(), namespace)
    config.write()


if __name__ == "__main__":
    prep()
