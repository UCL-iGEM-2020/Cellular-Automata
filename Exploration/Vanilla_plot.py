import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#from figure_gif2 import figure_gif
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits import mplot3d
import math

plotting_data = np.loadtxt('test1.txt')

print(plotting_data)

t_obs =  plotting_data[0,:]
cells_per_update = plotting_data[1,:]
total_lactate = plotting_data[2,:]
I = plotting_data[3,:]
j_max_list = plotting_data[4,:]
MT = plotting_data[5,:]
MaxI = plotting_data[6,:]

fig1 = plt.figure(figsize=(12,8))
ax1=fig1.add_subplot(2,1,1)
ax3=fig1.add_subplot(2,1,2)

fig2 = plt.figure(figsize=(12,8))
ax4=fig2.add_subplot(2,1,1)
ax5=fig2.add_subplot(2,1,2)

fig3 = plt.figure(figsize=(12,8))
ax6=fig3.add_subplot(2,1,1)
ax7=fig3.add_subplot(2,1,2)

ax1.set_xlabel('Time in hours')
ax1.set_ylabel('cells')
ax1.plot(t_obs,cells_per_update,'b-',label='Active cells')

ax3.set_xlabel('Time in hours')
ax3.set_ylabel('Total Lactate in the grid - mmol')
ax3.plot(t_obs,total_lactate,'b-')

ax4.plot(t_obs, I, 'b-', label='Total Current output')
ax4.set_xlabel("Hours")
ax4.set_ylabel("Total current density (A/m2)")
ax4.set_title("Total Current output")

ax5.plot(t_obs, j_max_list, 'b-', label='Maximum direct current density')
ax5.set_xlabel("Hours")
ax5.set_ylabel("Maximum direct transfer current density (A/m2)")
ax5.set_title("Maximum current density")

ax6.plot(t_obs, MT, 'b-', label='Maximum current density')
ax6.set_xlabel("Hours")
ax6.set_ylabel("Mediator Transfer current density (A/m2)")
ax6.set_title("Mediator Transfer current density")

ax7.plot(t_obs, MaxI, 'b-')
ax7.set_xlabel("Hours")
ax7.set_ylabel("Maximum total current density (A/m2)")
ax7.set_title("Maximum total current density")
plt.show()
