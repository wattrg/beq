import numpy as np
config.title = "boltzmann quiescent flow"

v = 200
L = 10
flow_time = 1.0 * L / v
n = 100
n_vel_increments = 200

def initial_condition(x):
    vi = v + 0.0*x
    T = 273 + 100.0*np.exp(-(x-L/2)**2/(2*0.1**2))    
    return FlowState(rho = 1.3, T = T, v = vi)

config.domain = Domain()
config.domain.number_cells = n
config.domain.length = L
config.domain.left_boundary = BoundaryCondition(
    # type=BoundaryType.Dirichlet, value=initial_condition(0)
    type = BoundaryType.Neumann
    # type = BoundaryType.Periodic
)
config.domain.right_boundary = BoundaryCondition(
    # type=BoundaryType.Dirichlet, value=initial_condition(L)
    type = BoundaryType.Neumann
    # type = BoundaryType.Periodic
)

config.gas_model = GasModel("air")

config.equation = Boltzmann()
config.equation.min_v = -3000
config.equation.max_v = 3000
config.equation.n_vel_increments = n_vel_increments

config.solver = RungeKutta()
config.solver.max_step = 30000
config.solver.max_time = flow_time
config.solver.print_frequency = 1000
config.solver.plot_frequency = flow_time / 100
# config.solver.plot_every_n_steps = 1
config.solver.cfl = 0.5

config.initial_condition = initial_condition
