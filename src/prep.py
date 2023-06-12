#! /usr/bin/python3

import json
import os
import sys

class _JsonData:
    def __init__(self):
        beq = os.environ["BEQ"]
        with open(f"{beq}/resources/defaults/{self._default}", "r") as default_json:
            defaults = json.load(default_json)

        for key in self._json_values:
            setattr(self, key, defaults[key])

    def to_dict(self):
        json_values = {}
        for value in self._json_values:
            json_values[value] = getattr(self, value)
        return json_values


class _Config(_JsonData):
    _json_values = [
        "title"
    ]

    _python_values = [
        "initial_condition",
        "solver",
        "domain",
    ]

    __slots__ = _json_values + _python_values

    _default = "config.json"

    def _write_json_config(self):
        if not os.path.exists("config"):
            os.mkdir("config")
        json_values = self.to_dict()
        json_values["solver"] = self.solver.to_dict()
        json_values["domain"] = self.domain.to_dict()
        json_values["equation"] = self.equation.to_dict()
        with open("config/config.json", "w") as file:
            json.dump(json_values, file, indent=4)

    def _write_initial_condition(self):
        if not os.path.exists("solution"):
            os.mkdir("solution")
        dx = self.domain.length / self.domain.number_cells
        with open("solution/phi_0.beq", "w") as ic:
            for i in range(self.domain.number_cells):
                x = (i + 0.5) * dx
                value = self.initial_condition(x)
                ic.write(f"{value}\n")

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
        "plot_frequency",
        "cfl",
    ]
    __slots__ = _json_values
    _default = "runge_kutta.json"
    _type = "runge_kutta"

    def __init__(self):
        super().__init__()
        self.type = self._type

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

    def __init__(self):
        self.type = self._type
    

def prep():
    prep_script = sys.argv[1]
    config = _Config()
    namespace = {
        "config": config,
        "RungeKutta": RungeKutta,
        "Domain": Domain,
        "Advection": Advection,
    }
    with open(prep_script, "r") as f:
        exec(f.read(), namespace)
    config.write()


if __name__ == "__main__":
    prep()
