#! /usr/bin/python3

import json
import os
import sys

from enum import Enum
import numpy as np

kB = 1.380649e-23 # J / K

def maxwellian_distribution(flow_state, vel, gas, volume):
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

    Returns
    -------
    list<float>
        The distribution function
    """
    mass = gas.mass
    T = flow_state.T
    norm = (mass/(2*np.pi*kB*T))**(1./2.)
    exponent = -mass/(2*kB*T)*(vel - flow_state.v)**2
    number_particles = flow_state.rho / mass * volume
    return number_particles * norm * np.exp(exponent)

class FlowState:
    _values = ["rho", "T", "v"]
    __slots__ = _values 

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            self.__setattr__(key, value)

    def to_dict(self):
        flow_state = {}
        for key in self._values:
            if hasattr(self, key):
                flow_state[key] = self.__getattribute__(key)
        return flow_state


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



class BoundaryType(Enum):
    Dirichlet = "dirichlet"
    Neumann = "neumann"
    Periodic = "periodic"

def string_to_boundary_type(string):
    if string == BoundaryType.Dirichlet.value:
        return BoundaryType.Dirichlet
    if string == BoundaryType.Neumann.value:
        return BoundaryType.Neumann
    if string == BoundaryType.Periodic.value:
        return BoundaryType.Periodic
    raise Exception(f"Unknown boundary type {string}")

class BoundaryCondition(_JsonData):
    _json_values = ["type", "value"]
    __slots__ = _json_values
    _default = "boundary_condition.json"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if type(self.type) != BoundaryType:
            self.type = string_to_boundary_type(self.type)

    def to_dict(self):
        return {
            "type": self.type.value,
            "value": self.value
        }


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
        self._transform_equilibrium_dirichlet_boundaries()

        json_values = self.to_dict()
        json_values["solver"] = self.solver.to_dict()
        json_values["domain"] = self.domain.to_dict()
        json_values["equation"] = self.equation.to_dict()
        if hasattr(self, "gas_model"):
            json_values["gas_model"] = self.gas_model.to_dict()

        if not os.path.exists("config"):
            os.mkdir("config")
        with open("config/config.json", "w") as file:
            json.dump(json_values, file, indent=4)

    def _transform_equilibrium_dirichlet_boundaries(self):
        if self.equation.field_type != FieldType.EQUILIBRIUM:
            return

        left = self.domain.left_boundary
        right = self.domain.right_boundary
        if left.type == BoundaryType.Dirichlet or right.type == BoundaryType.Dirichlet:
            min_v = self.equation.min_v
            max_v = self.equation.max_v
            nv = self.equation.n_vel_increments
            dv = (max_v - min_v) / nv
            vels = np.linspace(min_v + dv/2, max_v-dv/2, nv)
            dx = self.domain.length / self.domain.number_cells

        if left.type == BoundaryType.Dirichlet:
            dist = []
            for vel in vels:
                dist.append(maxwellian_distribution(left.value, vel, self.gas_model, dx))
            self.domain.left_boundary.value = dist

        if right.type == BoundaryType.Dirichlet:
            dist = []
            for vel in vels:
                dist.append(maxwellian_distribution(right.value, vel, self.gas_model, dx))
            self.domain.right_boundary.value = dist

    def _write_direct_initial_condition(self):
        dx = self.domain.length / self.domain.number_cells
        with open("solution/phi_0.beq", "w") as ic:
            for i in range(self.domain.number_cells):
                x = (i + 0.5) * dx
                value = self.initial_condition(x)
                ic.write(f"{value}\n")

    def _write_equilibrium_initial_condition(self):
        dx = self.domain.length / self.domain.number_cells
        min_v = self.equation.min_v
        max_v = self.equation.max_v
        nv = self.equation.n_vel_increments
        dv = (max_v - min_v) / nv
        vels = np.linspace(min_v + dv/2, max_v-dv/2, nv)
        with open("solution/phi_0.beq", "w") as ic:
            for vel in vels:
                for i in range(self.domain.number_cells):
                    x = (i + 0.5) * dx
                    flow_state = self.initial_condition(x)
                    dist = maxwellian_distribution(flow_state, vel, self.gas_model, dx)
                    ic.write(f"{dist}\n")

    def _write_initial_condition(self):
        if not os.path.exists("solution"):
            os.mkdir("solution")

        if self.equation.field_type == FieldType.EQUILIBRIUM:
            self._write_equilibrium_initial_condition()

        elif self.equation.field_type == FieldType.DIRECT:
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

class FieldType(Enum):
    DIRECT = 1
    EQUILIBRIUM = 2



class Domain(_JsonData):
    _json_values = [
        "number_cells",
        "length",
    ]

    _python_values = [
        "left_boundary",
        "right_boundary"
    ]

    __slots__ = _json_values + _python_values
    _default = "domain.json"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.left_boundary = BoundaryCondition()
        self.right_boundary = BoundaryCondition()

    def to_dict(self):
        domain = super().to_dict()
        domain["left_boundary"] = self.left_boundary.to_dict()
        domain["right_boundary"] = self.right_boundary.to_dict()
        return domain


class Advection(_JsonData):
    _json_values = [
        "type",
        "velocity"
    ]
    __slots__ = _json_values
    _default = "advection_eq.json"
    _type = "advection"
    field_type = FieldType.DIRECT

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
    field_type = FieldType.DIRECT

    def __init__(self, **kwargs):
        self.type = self._type
        super().__init__(**kwargs)

class Boltzmann(_JsonData):
    _json_values = [
        "type",
        "min_v",
        "max_v",
        "n_vel_increments",
    ]
    __slots__ = _json_values

    _default = "boltzmann_eq.json"
    _type = "boltzmann"
    field_type = FieldType.EQUILIBRIUM

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
        "BoundaryCondition": BoundaryCondition,
        "BoundaryType": BoundaryType,
    }
    with open(prep_script, "r") as f:
        exec(f.read(), namespace)
    config.write()


if __name__ == "__main__":
    prep()
