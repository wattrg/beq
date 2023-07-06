#! /usr/bin/python3

from enum import Enum
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import glob
import sys
import json

from prep import Config, kB

NA = 6.022e23 # particles / mol
R_UNIV = 8.314463 # J / K / mol

HELP = """BEQ post processing help

beq post action [options]

Available actions are:
    - animate: animate the data directly
    - phi_to_macroscopic: Convert data to macroscopic properties
    - animate_macroscopic: Animate the macroscopic properties 
"""

def main():
    if len(sys.argv) < 2:
        print(HELP)
        return

    action = get_action(sys.argv[1])
    if action == Action.ANIMATE_DIRECT:
        repeat = "--repeat" in sys.argv
        direct_animation(repeat)
    elif action == Action.PHI_TO_MACROSCOPIC:
        config = read_config()
        phi_to_macroscopic(config)
    elif action == Action.ANIMATE_MACROSCOPIC:
        variable = sys.argv[2]
        # if the macroscopic variables haven't been computed yet
        # we need to do that first
        if not os.path.exists("plot"):
            config = read_config()
            phi_to_macroscopic(config)

        plot_macroscopic(variable)
    elif action == Action.ANIMATE_DISTRIBUTION:
        config = read_config()
        cell_i = int(sys.argv[2])
        animate_distribution(config, cell_i)

def read_config():
    with open("config/config.json", "r") as config_file:
        config = json.load(config_file)
    return Config(config)

class Action(Enum):
    ANIMATE_DIRECT = 1
    PHI_TO_MACROSCOPIC = 2
    ANIMATE_MACROSCOPIC = 3
    ANIMATE_DISTRIBUTION = 4

def get_action(action_string):
    if action_string == "animate_direct":
        return Action.ANIMATE_DIRECT
    elif action_string == "phi_to_macroscopic":
        return Action.PHI_TO_MACROSCOPIC
    elif action_string == "animate_macroscopic":
        return Action.ANIMATE_MACROSCOPIC
    elif action_string == "animate_distribution":
        return Action.ANIMATE_DISTRIBUTION
    else:
        return Action.ANIMATE_DIRECT

def animate_distribution(config, cell_i, interval=20, repeat=False):
    data = read_phi()
    nc = config.domain.number_cells
    min_v = config.equation.min_v
    max_v = config.equation.max_v
    nv = config.equation.n_vel_increments
    dv = (max_v - min_v) / nv
    vels = np.linspace(min_v+dv/2, max_v-dv/2, nv)

    distributions = []
    for time_i in range(len(data)):
        phis = data[time_i].reshape((nv, nc))
        distributions.append(phis[:, cell_i])
    
    animate(distributions, repeat, interval, vels)

def phi_to_macroscopic(config):
    """
    Convert the distribution function to macroscopic properties for a simulation
    """
    data = read_phi()
    length = config.domain.length
    nc = config.domain.number_cells
    dx = length / nc
    volume = dx
    nv = config.equation.n_vel_increments
    mass = config.gas_model.mass

    # compute gas properties
    R = kB / mass # J / (kg * K)

    density = np.zeros(nc)
    momentum = np.zeros(nc)
    energy = np.zeros(nc)
    temperature = np.zeros(nc)
    pressure = np.zeros(nc)
    velocity = np.zeros(nc)

    variables = {
        "density": density,
        "momentum": momentum,
        "energy": energy,
        "temperature": temperature,
        "pressure": pressure,
        "velocity": velocity
    }
    if not os.path.exists("plot"):
        os.mkdir("plot")


    nv = config.equation.n_vel_increments
    max_v = config.equation.max_v
    min_v = config.equation.min_v
    dv = (max_v - min_v) / nv
    vels = np.linspace(min_v+dv/2, max_v-dv/2, nv)

    for time_i in range(len(data)):
        phis = data[time_i].reshape((nv, nc))
        for cell_i in range(nc):
            phi = phis[:, cell_i]
            density[cell_i] = mass * moment_of_distribution(vels, phi, 0) / volume
            momentum[cell_i] = mass * moment_of_distribution(vels, phi, 1) / volume
            velocity[cell_i] = momentum[cell_i] / density[cell_i]
            thermal_energy_x = mass * moment_of_distribution(vels, phi, 2, velocity[cell_i]) / (2 * volume)
            temperature[cell_i] = thermal_energy_x/(density[cell_i]*(1./2.)*R)
            energy[cell_i] = (3./2.) * R * temperature[cell_i]
            pressure[cell_i] = density[cell_i] * R * temperature[cell_i]

        if not os.path.exists(f"plot/{time_i}"):
            os.mkdir(f"plot/{time_i}")

        for var in variables:
            with open(f"plot/{time_i}/{var}.beq", "w") as f:
                np.savetxt(f, variables[var])

def moment_of_distribution(vels, phi, order, avg=0):
    """ 
    Compute a moment of the distribution function 

    Parameters
    ----------
    config: Config
        The configuration of the simulation
    phi: list<float>
        The distribution function
    order: int
        The order of the moment to calculate

    Returns
    -------
    float: The moment of the distribution function
    """
    return np.trapz((vels-avg)**order * phi, vels-avg)

def plot_macroscopic(var, repeat=False):
    folders = next(os.walk("plot"))[1]
    number_solutions = len(folders)

    data = []
    for i in range(number_solutions):
        data.append(np.loadtxt(f"plot/{i}/{var}.beq"))

    animate(data, repeat, 20)
    
def read_phi():
    """
    Read all the values of phi from a simultion.
    """
    files = glob.glob("solution/phi_*.beq")
    number_solutions = len(files)

    data = []
    for i in range(number_solutions):
        data.append(np.loadtxt(f"solution/phi_{i}.beq"))
    return data

def direct_animation(repeat=False):
    data = read_phi()
    animate(data, repeat, 20)
    
def animate(data, repeat, interval, x_data=None):
    fig, ax = plt.subplots()
    initial_data = data[0]

    if x_data is None:
        x_data = np.arange(0, len(initial_data), 1)

    line, = ax.plot(x_data, initial_data)
    title = ax.set_title(f"plot index 0")

    def init():
        return line,
    
    def animate(i):
        line.set_data(x_data, data[i])
        title.set_text(f"plot index {i}")
        return line, title

    _ = FuncAnimation(fig, animate, repeat=repeat, init_func=init,
                      frames=len(data), interval=interval)
    plt.show()

if __name__ == "__main__":
    main()    
