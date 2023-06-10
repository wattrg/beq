import math

n = 1000
L = 10
dx = L / n
u = 100
flow_time = L / u
sigma = 0.4
x0 = 2.0

def initial_condition(x):
    exponent = 0.5 * (x - x0) / sigma;
    return 1.0 / (sigma * math.sqrt(2*math.pi)) * math.exp(-exponent*exponent);

config.number_cells = n
config.length = L
config.velocity = u
config.max_step = 1000
config.max_time = flow_time
config.print_frequency = 100
config.plot_frequency = 10
config.cfl = 1.0
config.initial_condition = initial_condition
