import numpy as np
config.title = "shock wave"

kB = 1.380649e-23

L = 100
n = 200
n_vel_increments = 500

gm = GasModel("argon")
R = kB / gm.mass
M = 5
T = 300
rho = 1.0
v = M*np.sqrt(gm.gamma*R*T)

flow_time = 10.0 * L / v
print(f"flow time = {flow_time}")

def temp_ratio(mach, gamma):
    a = 2*gamma*mach**2 - (gamma-1)
    b = (gamma-1)*mach**2 + 2
    c = (gamma+1)**2 * mach**2
    return a*b/c

def density_ratio(mach, gamma):
    a = (gamma+1) * mach**2
    b = (gamma-1) * mach**2 + 2
    return a/b

def mach_ratio(mach, gamma):
    a = (gamma-1) * mach**2 + 2
    b = 2 * gamma * mach**2 - (gamma-1)
    return np.sqrt(a/b) / mach

def initial_condition(x):
    if x < L/2:
        return FlowState(rho=rho, T=T, v=v)
    else:
        M_2 = mach_ratio(M, gm.gamma)*M
        T_2 = T * temp_ratio(M, gm.gamma)
        v_2 = M_2 * np.sqrt(gm.gamma*R*T_2)
        rho_2 = rho * density_ratio(M, gm.gamma)
        return FlowState(rho=rho_2, T=T_2, v = v_2)

config.domain = Domain()
config.domain.number_cells = n
config.domain.length = L
config.domain.left_boundary = BoundaryCondition(
    type=BoundaryType.Dirichlet, value=initial_condition(0)
)
config.domain.right_boundary = BoundaryCondition(
    type=BoundaryType.Dirichlet, value=initial_condition(L)
)

config.gas_model = gm

config.equation = Boltzmann()
config.equation.min_v = -5000
config.equation.max_v = 5000
config.equation.n_vel_increments = n_vel_increments

config.solver = RungeKutta()
config.solver.max_step = 1000000
config.solver.max_time = flow_time
config.solver.print_frequency = 1000
config.solver.plot_frequency = flow_time / 100
# config.solver.plot_every_n_steps = 1
config.solver.cfl = 0.5

config.initial_condition = initial_condition
