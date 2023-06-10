#! /usr/bin/python3

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import glob
import sys

files = glob.glob("solution/phi_*.beq")
number_solutions = len(files)

fig, ax = plt.subplots()
initial_data = np.loadtxt("solution/phi_0.beq")
line, = ax.plot(initial_data)

def init():
    return line,

def animate(i):
    data = np.loadtxt(f"solution/phi_{i}.beq")
    x_data = np.arange(0, len(data), 1)
    line.set_data(x_data, data)
    return line,

if "--repeat" in sys.argv:
    repeat = True
else:
    repeat = False
anim = FuncAnimation(fig, animate, repeat=repeat, init_func=init, frames=number_solutions, interval = 20, blit=True)
plt.show()
