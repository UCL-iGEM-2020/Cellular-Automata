import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math

plotting_data = np.loadtxt('/Users/Oliver_Leblanc/Desktop/RUNS/RUN2/1week4.txt')

t_obs =  plotting_data[0,:]
cells_per_update = plotting_data[1,:]
total_lactate = plotting_data[2,:]
I = plotting_data[3,:]
j_max_list = plotting_data[4,:]
MT = plotting_data[5,:]
MaxI = plotting_data[6,:]
thickness = plotting_data[7,:]
total_ox_conc = plotting_data[8,:]
total_red_conc = plotting_data[9,:]
total_med_conc = plotting_data[10,:]
active_cells = plotting_data[11,:]
quiescent_cells = plotting_data[12,:]
dead_cells = plotting_data[13,:]

t_obs = np.delete(t_obs,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
cells_per_update = np.delete(cells_per_update,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
total_lactate = np.delete(total_lactate,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
I = np.delete(I,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
j_max_list = np.delete(j_max_list,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
MT = np.delete(MT,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
MaxI = np.delete(MaxI,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
thickness = np.delete(thickness,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
total_ox_conc = np.delete(total_ox_conc,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
total_red_conc = np.delete(total_red_conc,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
total_med_conc = np.delete(total_med_conc,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
active_cells = np.delete(active_cells,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
quiescent_cells = np.delete(quiescent_cells,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
dead_cells = np.delete(dead_cells,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])

fig1 = plt.figure(figsize=(12,8))
ax1=fig1.add_subplot(1,1,1)

fig2 = plt.figure(figsize=(12,8))
ax2=fig2.add_subplot(1,1,1)

fig3 = plt.figure(figsize=(12,8))
ax3=fig3.add_subplot(1,1,1)

fig4 = plt.figure(figsize=(12,8))
ax4=fig4.add_subplot(1,1,1)

fig5 = plt.figure(figsize=(12,8))
ax5=fig5.add_subplot(1,1,1)

fig6 = plt.figure(figsize=(12,8))
ax6=fig6.add_subplot(1,1,1)

fig7 = plt.figure(figsize=(12,8))
ax7 = fig7.add_subplot(1,1,1)

fig8 = plt.figure(figsize=(12,8))
ax8 = fig8.add_subplot(1,1,1)

fig9 = plt.figure(figsize=(12, 8))
ax9 = fig9.add_subplot(1,1,1)

fig10 = plt.figure(figsize=(12, 8))
ax10 = fig10.add_subplot(1,1,1)

fig11 = plt.figure(figsize=(12, 8))
ax11 = fig11.add_subplot(1,1,1)

ax1.set_xlabel('Hours')
ax1.set_ylabel('cells')
ax1.plot(t_obs,cells_per_update, "#29cce6", label='Active cells')

ax2.set_xlabel('Hours')
ax2.set_ylabel('Total Lactate in the grid - mmol')
ax2.plot(t_obs,total_lactate, "#29cce6")
ax2.set_title("Lactate concentration")


ax3.plot(t_obs, I,"#29cce6", label='Total Current output')
ax3.set_xlabel("Hours")
ax3.set_ylabel("Total current density (A/m2)")
ax3.set_title("Total Current output")

ax4.plot(t_obs, j_max_list, "#29cce6", label='Maximum direct current density')
ax4.set_xlabel("Hours")
ax4.set_ylabel("Maximum current density (A/m2)")
ax4.set_title("Maximum Direct Transfer current density")

ax5.plot(t_obs, MT, "#29cce6", label='Maximum current density')
ax5.set_xlabel("Hours")
ax5.set_ylabel("Mediator Transfer current density (A/m2)")
ax5.set_title("Mediator Transfer current density")

ax6.plot(t_obs, MaxI, "#29cce6")
ax6.set_xlabel("Hours")
ax6.set_ylabel("Maximum total current density (A/m2)")
ax6.set_title("Maximum total current density")


ax7.plot(t_obs, thickness, "#29cce6")
ax7.set_xlabel("Hours")
ax7.set_ylabel("Average biofilm thickness (m)")
ax7.set_title("Biofilm Thickness")


ax8.plot(t_obs, total_ox_conc, "#29cce6", label='Oxidized mediator')
ax8.plot(t_obs, total_red_conc, "#043854", label='Reduced mediator')
ax8.legend()
ax8.set_xlabel("Hours")
ax8.set_ylabel("Mediator concentration (mmol)")
ax8.set_title("Mediator Concentration")

ax9.plot(t_obs, total_med_conc, "#29cce6")
ax9.set_xlabel("Hours")
ax9.set_ylabel("Mediator concentration (mmol)")
ax9.set_title("Total mediator Concentration")

ax10.plot(t_obs, active_cells, label='Active cells')
ax10.plot(t_obs, quiescent_cells, label='Quiescent cells')
ax10.plot(t_obs, dead_cells, label='Dead cells')
ax10.legend()
ax10.set_xlabel("Hours")
ax10.set_ylabel("Agents")
ax10.set_title("Agent type over time")

#ax12.plot(t_obs, I)
#ax12.plot(t_lim, I_t(t_lim))
#ax12.set_xlabel("Hours")
#ax12.set_ylabel("rate of desalination (litres/hour)")
#ax12.set_title("rate of desalination over time")
 
plt.show()





