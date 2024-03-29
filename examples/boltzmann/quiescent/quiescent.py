config.title = "boltzmann quiescent flow"

v = 10
L = 1
flow_time = 2 * L / v
n = 10

n_vel_increments = 20

def initial_condition(x):
    return FlowState(rho = 1.3, T = 273, v = v)


config.domain = Domain()
config.domain.number_cells = n
config.domain.length = L
config.domain.left_boundary = BoundaryCondition(
    # type = BoundaryType.Periodic,
    type = BoundaryType.Neumann,
    # type=BoundaryType.Dirichlet, value = initial_condition(0)
)
config.domain.right_boundary = BoundaryCondition(
    # type = BoundaryType.Periodic,
    type = BoundaryType.Neumann,
    # type=BoundaryType.Dirichlet, value = initial_condition(L)
)

config.gas_model = GasModel("air")

config.equation = Boltzmann()
config.equation.min_v = -2000
config.equation.max_v = 2000
config.equation.n_vel_increments = n_vel_increments
config.equation.collision_operator = BGK()

config.solver = RungeKutta()
config.solver.max_step = 10000
config.solver.max_time = flow_time
config.solver.print_frequency = 1000
config.solver.plot_frequency = flow_time / 20
config.solver.cfl = 0.1

config.initial_condition = initial_condition
