config.title = "boltzmann quiescent flow"

v = 10
L = 10
flow_time = 1 # L / v
n = 12

n_vel_increments = 300

def initial_condition(x):
    T = 273 + 1.0 * x    
    return FlowState(rho = 1.3, T = T, v = v)

config.domain = Domain()
config.domain.number_cells = n
config.domain.length = L
config.domain.left_boundary = BoundaryCondition(
    type=BoundaryType.Dirichlet, value=initial_condition(0)
)
config.domain.right_boundary = BoundaryCondition(
    type=BoundaryType.Dirichlet, value=initial_condition(L)
)

config.gas_model = GasModel("air")

config.equation = Boltzmann()
config.equation.min_v = -4000
config.equation.max_v = 4000
config.equation.n_vel_increments = n_vel_increments

config.solver = RungeKutta()
config.solver.max_step = 1000000
config.solver.max_time = flow_time
config.solver.print_frequency = 100
config.solver.plot_frequency = flow_time / 20
config.solver.cfl = 0.5

config.initial_condition = initial_condition
