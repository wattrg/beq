import math

n = 1000
L = 10
dx = L / n
sigma = 0.4
x0 = 2.0

def initial_condition(x):
    exponent = 0.5 * (x - x0) / sigma;
    return 1.0 / (sigma * math.sqrt(2*math.pi)) * math.exp(-exponent*exponent);

config.number_cells = n
config.length = L
config.velocity = 100
config.max_step = 1000
config.max_time = 0.1
config.print_frequency = 100
config.plot_frequency = 10
config.cfl = 1.0
config.initial_condition = initial_condition
