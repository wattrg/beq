#! /usr/bin/python3

from enum import Enum
import os
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import glob
import sys
import json

from prep import Config

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
        with open("config/config.json", "r") as config_file:
            config = json.load(config_file)
        config = Config(config)
        phi_to_macroscopic(config)
    elif action == Action.ANIMATE_MACROSCOPIC:
        variables = sys.argv[2:]
        plot_macroscopic(variables)
    elif action == Action.ANIMATE_DISTRIBUTION:
        with open("config/config.json", "r") as config_file:
            config = json.load(config_file)
        config = Config(config)
        cell_i = int(sys.argv[2])
        animate_distribution(config, cell_i)

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
    nv = config.equation.n_vel_increments * 2
    dv = config.equation.dv
    n_vel_inc = config.equation.n_vel_increments
    max_v = (n_vel_inc + 0.5) * dv
    vels = np.arange(0.5*dv, max_v, dv)
    vels = np.append(vels, -vels)

    distributions = []
    for time_i in range(len(data)):
        phis = data[time_i].reshape((nc, nv))
        distributions.append(phis[cell_i, :])
    
    animate(distributions, repeat, interval, vels)

def phi_to_macroscopic(config):
    """
    Convert the distribution function to macroscopic properties for a simulation
    """
    data = read_phi()
    length = config.domain.length
    nc = config.domain.number_cells
    dx = length / nc
    nv = config.equation.n_vel_increments * 2
    mass = config.gas_model.mass
    molar_mass = mass * NA

    # compute gas properties
    R = R_UNIV / molar_mass
    Cv = (3 / 2) * R

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


    for time_i in range(len(data)):
        phis = data[time_i].reshape((nc, nv))
        for cell_i in range(nc):
            phi = phis[cell_i]
            # moments of the distribution function
            density[cell_i] = mass * moment_of_distibution(config, phi, 0) / dx
            momentum[cell_i] = mass * moment_of_distibution(config, phi, 1) / dx
            energy[cell_i] = mass * moment_of_distibution(config, phi, 2) / 2 / dx

            # derived quantities
            temperature[cell_i] = energy[cell_i] / Cv            
            velocity[cell_i] = momentum[cell_i] / mass
            pressure = density[cell_i] * R * temperature[cell_i]


        if not os.path.exists(f"plot/{time_i}"):
            os.mkdir(f"plot/{time_i}")

        for var in variables:
            with open(f"plot/{time_i}/{var}.beq", "w") as f:
                np.savetxt(f, variables[var])


def moment_of_distibution(config, phi, order):
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
    dv = config.equation.dv
    nv = config.equation.n_vel_increments
    max_v = (nv + 0.5) * dv
    vels = np.arange(0.5*dv, max_v, dv)
    vels = np.append(vels, -vels)
    nv = len(vels)

    total = 0
    for i in range(nv):
        total += vels[i]**order * phi[i] * dv
    return total

def plot_macroscopic(var, repeat=False):
    folders = next(os.walk("plot"))[1]
    number_solutions = len(folders)

    data = []
    for i in range(number_solutions):
        data.append(np.loadtxt(f"plot/{i}/{var[0]}.beq"))

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

    def init():
        return line,
    
    def animate(i):
        line.set_data(x_data, data[i])
        return line,

    _ = FuncAnimation(fig, animate, repeat=repeat, init_func=init,
                      frames=len(data), interval=interval, blit=True)
    plt.show()

if __name__ == "__main__":
    main()    
