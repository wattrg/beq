#! /usr/bin/python3

import json
import os
import sys

class _Config:
    _json_values = [
        "number_cells",
        "length",
        "velocity",
        "max_step",
        "max_time",
        "print_frequency",
        "plot_frequency",
        "cfl",
    ]

    _python_values = [
        "initial_condition"
    ]

    __slots__ = _json_values + _python_values

    def __init__(self):
        beq = os.environ["BEQ"]
        with open(f"{beq}/resources/defaults.json", "r") as default_json:
            defaults = json.load(default_json)

        for key in defaults:
            setattr(self, key, defaults[key])

    def _write_json_config(self):
        if not os.path.exists("config"):
            os.mkdir("config")
        json_values = {}
        for value in self._json_values:
            json_values[value] = getattr(self, value)
        with open("config/config.json", "w") as file:
            json.dump(json_values, file, indent=4)

    def _write_initial_condition(self):
        if not os.path.exists("solution"):
            os.mkdir("solution")
        dx = self.length / self.number_cells
        with open("solution/phi_0.beq", "w") as ic:
            for i in range(self.number_cells):
                x = (i + 0.5) * dx
                value = self.initial_condition(x)
                ic.write(f"{value}\n")

    def write(self):
        self._write_json_config()
        self._write_initial_condition()


def prep():
    prep_script = sys.argv[1]
    config = _Config()
    namespace = {"config": config}
    with open(prep_script, "r") as f:
        exec(f.read(), namespace)
    config.write()


if __name__ == "__main__":
    prep()
