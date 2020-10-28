import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from figure_gif2 import figure_gif
from mpl_toolkits.axes_grid1 import make_axes_locatable

def diffusion_update(diff_grid,diff_calc,grid,dt,h):
    D = 3.51*10**(-1)
    for i in range(grid):
        for j in range(grid):
            if i==0 and j==0:
                diff_calc[i,j] = ((D*dt*(diff_grid[i+1,j] + diff_grid[i,grid-1] + diff_grid[i,j+1] + diff_grid[grid-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
            elif i==grid-1 and j==grid-1:
                diff_calc[i,j] = ((D*dt*(diff_grid[0,j] + diff_grid[i,grid-1] + diff_grid[i,0] + diff_grid[i-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
            elif i==0 and j==grid-1:
                diff_calc[i,j] = ((D*dt*(diff_grid[i+1,j] + diff_grid[i,j-1] + diff_grid[i,0] + diff_grid[grid-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
            elif i==grid-1 and j==0:
                diff_calc[i,j] = ((D*dt*(diff_grid[0,j] + diff_grid[i,grid-1] + diff_grid[i,j+1] + diff_grid[i-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
            elif i==0: 
                diff_calc[i,j] = ((D*dt*(diff_grid[i+1,j] + diff_grid[i,j-1] + diff_grid[i,j+1] + diff_grid[grid-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
            elif j==0:
                diff_calc[i,j] = ((D*dt*(diff_grid[i+1,j] + diff_grid[i,grid-1] + diff_grid[i,j+1] + diff_grid[i-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
            elif i==grid-1:
                diff_calc[i,j] = ((D*dt*(diff_grid[0,j] + diff_grid[i,j-1] + diff_grid[i,j+1] + diff_grid[i-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
            elif j==grid-1:
                diff_calc[i,j] = ((D*dt*(diff_grid[i+1,j] + diff_grid[i,j-1] + diff_grid[i,0] + diff_grid[i-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
            else:
                diff_calc[i,j] = ((D*dt*(diff_grid[i+1,j] + diff_grid[i,j-1] + diff_grid[i,j+1] + diff_grid[i-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method

    return diff_calc

def diffusion_plot_relative(grid,t):
    
    initial_conc=10
    dt=0.5
    h=1
    
    diff_grid = np.zeros((grid,grid)) #grid with lactate concentrations
    diff_calc = np.zeros((grid,grid)) #grid that stores changes
    
    diff_grid[round(grid/2),round(grid/2)]=initial_conc
    diff_calc[round(grid/2),round(grid/2)]=initial_conc

    
    anim=figure_gif(figsize=(7.5,7.5))
    ax=anim.add_subplot(1,1,1)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
        
    for i in range(0,round((t)*1/dt)):  
        ax.clear() 
        diff_plot=ax.imshow(diff_grid)
        cb=plt.colorbar(diff_plot, cax=cax, orientation='vertical')
        cb.outline.set_visible(False)
        ax.set_title(f'generation number = {i:.2f}')
        anim.add_frame()
        diff_calc = np.zeros((grid,grid))
        diff_calc = diffusion_update(diff_grid,diff_calc,grid,dt,h)
        diff_grid = diff_calc
        
    anim.show()
    
diffusion_plot_relative(50,25)
