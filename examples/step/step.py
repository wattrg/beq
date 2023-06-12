import math

config.title = "Step"

n = 100
L = 1
u = 1.0
flow_time = L / u
dx = L / n

def initial_condition(x):
    if x < L / 3 or x > 2 * L / 3:
        return 1.0
    else:
        return 0.0

config.solver = RungeKutta()
config.solver.max_step = 100000
config.solver.max_time = 1 * flow_time
config.solver.print_frequency = 10
config.solver.plot_frequency = 10
config.solver.cfl = 1.0

config.domain = Domain()
config.domain.number_cells = n
config.domain.length = L

config.equation = Advection()
config.equation.velocity = u

config.initial_condition = initial_condition
