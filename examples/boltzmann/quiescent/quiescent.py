config.title = "boltzmann quiescent flow"

v = 1
L = 10
flow_time = L / v
n = 100

def initial_condition(x):
    return FlowState(T = 300, v = 1)


config.domain = Domain()
config.domain.number_cells = n
config.domain.length = L

config.gas_model = GasModel("air")

config.equation = Boltzmann()

config.solver = RungeKutta()
config.solver.max_step = 10000
config.solver.max_time = flow_time
config.solver.print_frequency = 100
config.solver.plot_frequency = flow_time / 20
config.solver.cfl = 1.0

config.initial_condition = initial_condition
