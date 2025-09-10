
import numpy as np
import matplotlib.pyplot as plt
# AST imports
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission

seed = utils.get_seed('bmthune')
Mission = SpaceMission(seed)
system = Mission.system

#Henter data fra system
a = system.semi_major_axes
e = system.eccentricities

angles = np.linspace(0,2*np.pi, 1000)


r = np.zeros((len(a),2, len(angles)))
for i in range(len(a)):        
        r[i][0] = ((a[i]*(1-e[i]**2))/(1+e[i]*np.cos(angles))) * np.cos(angles)
        r[i][1] = ((a[i]*(1-e[i]**2))/(1+e[i]*np.cos(angles))) * np.sin(angles)

fig, ax = plt.subplots()
# 1 = 1, 2= 2, 3= 5 , 4 = 4, 5= ?, 6= ? 7=3, 
colours = ["red", 'darkorange', 'mediumblue', 'limegreen', 'purple', 'darkviolet', 'gold']
#plot
for i in range(len(r)):
    p = r[i]
    print(p[0],p[1])
    ax.plot(p[0],p[1], color = colours[i])
sol = plt.Circle((0, 0), 1, color = 'gold')
plt.axis('equal')
ax.add_patch(sol)
plt.show()