import math

n = 1000
L = 10
dx = L / n
flow_time = 1.5 * L / 1
sigma = 0.2
x0 = L / 3

def initial_condition(x):
    exponent = 0.5 * (x - x0) / sigma;
    return 1.0 / (sigma * math.sqrt(2*math.pi)) * math.exp(-exponent*exponent);


config.domain = Domain()
config.domain.number_cells = n
config.domain.length = L

config.equation = Burgers()

config.solver = RungeKutta()
config.solver.max_step = 10000
config.solver.max_time = flow_time
config.solver.print_frequency = 100
config.solver.plot_frequency = flow_time / 50
config.solver.cfl = 1.0
config.initial_condition = initial_condition
