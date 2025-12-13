import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.utils as utils
from   ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as const
from   ast2000tools.relativity import RelativityExperiments

LorentzFactor   = lambda v : 1/(1-v**2)

KgToMeters      = lambda kg : kg * (const.G / (const.c**2))

SecondsToMeters = lambda s : s*const.c

"""
PART 1)
"""

print("4)\n")

E_over_m = lambda M,r,v_shell : ((1-(2*M)/r)**0.5)*LorentzFactor(v_shell)

M   = KgToMeters(2.8113e7 * const.m_sun)
r_0 = 1 * const.AU
v_0 = 0.214

print(f"M = {M:e} meters")
print(f"r = {r_0:e} meters\n")

E_m_0 = E_over_m(M,r_0,v_0)
print(f"E/m natural units: {E_m_0}")
print(f"E/m SI units     : {E_m_0 * const.c**2}\n")

print("7)\n")

t_1 = 20.0138
t_2 = 43.3633

t_n_1 = 1020.71
t_n   = 1444.33

d_tau = SecondsToMeters(20.3634)
d_t_1 = SecondsToMeters(t_2-t_1)
d_t_n = SecondsToMeters(t_n - t_n_1)

print(f"d_tau = {d_tau:e} meters")
print(f"d_t_1 = {d_t_1:e} meters")
print(f"d_t_n = {d_t_n:e} meters\n")

r = lambda M,d_tau,d_t_shell,E_m,r_0 : 2*M * (1 - (d_tau/d_t_shell) * E_m * ((1 - (2*M)/r_0)**0.5))**(-1)

dr_1 = r(M,d_tau,d_t_1,E_m_0,r_0)
dr_n = r(M,d_tau,d_t_n,E_m_0,r_0)

print(f"dr_1 = {dr_1 / const.AU:e} AU")
print(f"dr_n = {dr_n / const.AU:e} AU\n")

print(f"dr_1 = {dr_1 / (2*M):e} r_s")
print(f"dr_n = {dr_n / (2*M):e} r_s\n")

"""
PART 2)
"""

print("3)\n")

x ,y = np.loadtxt("black_hole_descent_frame_1.txt")
x_l ,y_l = np.loadtxt("black_hole_descent_frame_1_with_light_travel.txt")

diffs = lambda x : x[2:] - x[1:-1]

dy = diffs(y)
dy_l = diffs(y_l)

ax = plt.axes()


ax.plot(x[2:],dy,color = "blue",label = "no light travel")
ax.plot(x_l[2:],dy_l,color = "green",label = "light travel")
ax.legend()

plt.title("Time differences for frame 1")
plt.xlabel("Ray index")
plt.ylabel("Time difference [s]")
print("plotting...\n")
plt.show()

print("6) \n")

x ,y = np.loadtxt("black_hole_descent_frame_2.txt")
x_l ,y_l = np.loadtxt("black_hole_descent_frame_2_with_light_travel.txt")

diffs = lambda x : x[2:] - x[1:-1]

dy = diffs(y)
dy_l = diffs(y_l)

ax = plt.axes()


ax.plot(x[2:],dy,color = "red",label = "no light travel")
ax.plot(x_l[2:],dy_l,color = "orange",label = "light travel")
ax.legend()

plt.title("Time differences for frame 2  ")
plt.xlabel("Ray index")
plt.ylabel("Time difference [s]")
print("plotting...\n")
plt.show()