import math

n = 100
L = 10
dx = L / n
u = 100
flow_time = L / u
sigma = 0.4
x0 = L / 2

def initial_condition(x):
    exponent = 0.5 * (x - x0) / sigma;
    return 1.0 / (sigma * math.sqrt(2*math.pi)) * math.exp(-exponent*exponent);


config.domain = Domain()
config.domain.number_cells = n
config.domain.length = L
config.domain.left_boundary = BoundaryCondition(type=BoundaryType.Periodic)
config.domain.right_boundary = BoundaryCondition(type=BoundaryType.Periodic)

config.equation = Advection()
config.equation.velocity = u

config.solver = RungeKutta()
config.solver.max_step = 10000
config.solver.max_time = flow_time
config.solver.print_frequency = 100
config.solver.plot_frequency = flow_time / 20
config.solver.cfl = 1.0
config.initial_condition = initial_condition
