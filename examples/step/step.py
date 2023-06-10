import math

n = 100
L = 1
u = 1.0
flow_time = L / u
dx = L / n
sigma = 0.4
x0 = 2.0

def initial_condition(x):
    if x < L / 3 or x > 2 * L / 3:
        return 1.0
    else:
        return 0.0

config.number_cells = n
config.length = L
config.velocity = u
config.max_step = 10000
config.max_time = 10 * flow_time
config.print_frequency = 100
config.plot_frequency = 1
config.cfl = 1.0
config.initial_condition = initial_condition
