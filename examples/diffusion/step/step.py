import math

config.title = "Step"

n = 100
L = 1
u = 1.0
flow_time = L / u
dx = L / n

def initial_condition(x):
    if x < L /2:
        return 0.0
    else:
        return 1.0

config.solver = RungeKutta()
config.solver.max_step = 100000
config.solver.max_time = 1 * flow_time
config.solver.print_frequency = 10
config.solver.plot_frequency = flow_time / 20
config.solver.cfl = 1.0

config.domain = Domain()
config.domain.number_cells = n
config.domain.length = L
config.domain.left_boundary = BoundaryCondition(type=BoundaryType.Dirichlet, value=[0])
config.domain.right_boundary = BoundaryCondition(type=BoundaryType.Dirichlet, value=[0.5])

config.equation = Diffusion()
config.equation.velocity = u

config.initial_condition = initial_condition
config.gas_model = GasModel("air")
