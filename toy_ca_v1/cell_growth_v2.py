import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from figure_gif2 import figure_gif
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

def diffusion_update(diff_grid,diff_calc,grid,dt,h):
    D = 3.51*10**(-10)
    for i in range(grid):
        for j in range(grid):
            if i==0 and j==0:
                diff_calc[i,j] = ((D*dt*(diff_grid[i+1,j] + diff_grid[i,grid-1] + diff_grid[i,j+1] + diff_grid[grid-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
            elif i==grid-1 and j==grid-1:
                diff_calc[i,j] = ((D*dt*(diff_grid[0,j] + diff_grid[i,j-1] + diff_grid[i,0] + diff_grid[i-1,j] - 4*diff_grid[i,j]))/h**2) + diff_grid[i,j] #finite difference method
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

def calc_poisson(average_doubling_time,number_of_updates)
    time = number_of_updates*0.01
    poisson = (((math.e)**(-average_doubling_time))*((average_doubling_time)**(time))/math.factorial(time))
    

    return poisson

def UpdateGrid_v2(M, grid, diff_grid, initial_conc, division_threshold):
    can_divide=[]
    # Generate list of cells that can divide
    for i in range(grid):
        for j in range(grid):
            #leftover_choices=[0,1,2,3]
            if M[i,j] == 1:
                if diff_grid[i,j] > (division_threshold):
                    coordinates=[]
                    if i==0 and j==0:
                        if M[i+1,j] == 0 or M[grid-1,j] == 0 or M[i,j+1] == 0 or M[i,grid-1] == 0:
                            coordinates.append(i)
                            coordinates.append(j)
                            can_divide.append(coordinates)
                    elif i==grid-1 and j==grid-1:
                        if M[0,j] == 0 or M[i-1,j] == 0 or M[i,0] == 0 or M[i,j-1] == 0:
                            coordinates.append(i)
                            coordinates.append(j)
                            can_divide.append(coordinates)
                    elif i==0 and j==grid-1:
                        if M[i+1,j] == 0 or M[grid-1,j] == 0 or M[i,0] == 0 or M[i,j-1] == 0:
                            coordinates.append(i)
                            coordinates.append(j)
                            can_divide.append(coordinates)
                    elif i==grid-1 and j==0:
                        if M[0,j] == 0 or M[i-1,j] == 0 or M[i,j+1] == 0 or M[i,j-1] == 0:
                            coordinates.append(i)
                            coordinates.append(j)
                            can_divide.append(coordinates)
                    elif i==(grid-1):
                        if M[0,j] == 0 or M[i-1,j] == 0 or M[i,j+1] == 0 or M[i,j-1] == 0:
                            coordinates.append(i)
                            coordinates.append(j)
                            can_divide.append(coordinates)
                    elif i == 0:
                        if M[i+1,j] == 0 or M[grid-1,j] == 0 or M[i,j+1] == 0 or M[i,j-1] == 0:
                            coordinates.append(i)
                            coordinates.append(j)
                            can_divide.append(coordinates)
                    elif j == (grid-1):
                        if M[i+1,j] == 0 or M[i-1,j] == 0 or M[i,0] == 0 or M[i,j-1] == 0:
                            coordinates.append(i)
                            coordinates.append(j)
                            can_divide.append(coordinates)
                    elif j == 0:
                        if M[i+1,j] == 0 or M[i-1,j] == 0 or M[i,j+1] == 0 or M[i,grid-1] == 0:
                            coordinates.append(i)
                            coordinates.append(j)
                            can_divide.append(coordinates)
                    else:
                        if M[i+1,j] == 0 or M[i-1,j] == 0 or M[i,j+1] == 0 or M[i,j-1] == 0:
                            coordinates.append(i)
                            coordinates.append(j)
                            can_divide.append(coordinates)

    #choose a random subset of cells to divide in that timestep          
    for cell in round(poisson*(len(can_divide)):
        divider = random.choice(can_divide)
        i=divider[0]
        j=divider[1]
        if len(can_divide) == 0:
            break
        can_divide.remove(divider)
        leftover_choices=[0,1,2,3]
        choice=random.choice(leftover_choices)
            if choice == 0:
                if i==(grid-1):
                    if M[0,j] == 0:
                        M[0,j] = 1
                        diff_grid[i,j]=diff_grid[i,j]-(division_threshold)
                    else:
                        leftover_choices.remove(0)
                        if len(leftover_choices) == 0:
                            continue
                        choice=random.choice(leftover_choices)
                elif M[i+1,j] == 0:
                    M[i+1,j] = 1
                    diff_grid[i,j]=diff_grid[i,j]-(division_threshold)
                else:
                    leftover_choices.remove(0)
                    if len(leftover_choices) == 0:
                        continue
                    choice=random.choice(leftover_choices)
            if choice == 1:
                if i == 0:
                    if M[grid-1,j] == 0:
                        M[grid-1,j] = 1
                        diff_grid[i,j]=diff_grid[i,j]-(division_threshold)
                    else:
                        leftover_choices.remove(1)
                        if len(leftover_choices) == 0:
                            continue
                        choice=random.choice(leftover_choices)
                elif M[i-1,j] == 0:
                    M[i-1,j] = 1
                    diff_grid[i,j]=diff_grid[i,j]-(division_threshold)
                else:
                    leftover_choices.remove(1)
                    if len(leftover_choices) == 0:
                        continue
                    choice=random.choice(leftover_choices)
            if choice == 2:
                if j == (grid-1):
                    if M[i,0] == 0:
                        M[i,0] = 1
                        diff_grid[i,j]=diff_grid[i,j]-(division_threshold)
                    else:
                        leftover_choices.remove(2)
                        if len(leftover_choices) == 0:
                            continue
                        choice=random.choice(leftover_choices)
                elif M[i,j+1] == 0:
                    M[i,j+1] = 1
                    diff_grid[i,j]=diff_grid[i,j]-(division_threshold)
                else:
                    leftover_choices.remove(2)
                    if len(leftover_choices) == 0:
                        continue
                    choice=random.choice(leftover_choices)
            if choice == 3:
                if j == 0:
                    if M[i,grid-1] == 0:
                        M[i,grid-1] = 1
                        diff_grid[i,j]=diff_grid[i,j]-(division_threshold)
                    else:
                        leftover_choices.remove(3)
                        if len(leftover_choices) == 0:
                            continue
                        choice=random.choice(leftover_choices)
                elif M[i,j-1] == 0:
                    M[i,j-1] = 1
                    diff_grid[i,j]=diff_grid[i,j]-(division_threshold)
                else:
                    leftover_choices.remove(3)
                    if len(leftover_choices) == 0:
                        continue
                    choice=random.choice(leftover_choices)
    return M

def life_v2(grid, ngen,t):
    
    #Setting up initial lactate diffusion conditions and plotting
    
    initial_conc=100
    survival_threshold=0.005
    division_threshold=0.01

    dt=0.01 # can't be bigger than 0.01780626781 secs for h = 5*(10**(-6))
    h=5*(10**(-6))
    
    diff_grid = np.zeros((grid,grid)) #grid with lactate concentrations
    diff_calc = np.zeros((grid,grid)) #grid that stores changes
    
    diff_grid[round(grid/2),round(grid/2)]=initial_conc
    diff_calc[round(grid/2),round(grid/2)]=initial_conc

    
    anim1=figure_gif(figsize=(7.5,7.5))
    ax1=anim1.add_subplot(1,1,1)
    
    divisor = make_axes_locatable(ax1)
    cax = divisor.append_axes('right', size='5%', pad=0.05)
        
    for i in range(0,round((t)*1/dt)):  
        ax1.clear() 
        diff_plot=ax1.imshow(diff_grid)
        cb=plt.colorbar(diff_plot, cax=cax, orientation='vertical')
        cb.outline.set_visible(False)
        ax1.set_title(f'Number of timesteps = {i:.2f}')
        anim1.add_frame()
        diff_calc = np.zeros((grid,grid))
        diff_calc = diffusion_update(diff_grid,diff_calc,grid,dt,h)
        diff_grid = diff_calc

    
    #Initialize agent and number of neighbours
    M = np.zeros((grid,grid))
    M[round(grid/2),round(grid/2)]=1
    nneigh = np.zeros((grid,grid))
    
    anim2=figure_gif(figsize=(15,7.5))
    ax2=anim2.add_subplot(1,1,1)
    
    fig1 = plt.figure(figsize=(12,8))
    ax3=fig1.add_subplot(1,1,1)
    
    cells_per_gen=[]
    t_obs = np.linspace(0,ngen,ngen)

    # loop over the results, drawing results
    # into the figure and saving it as an
    # animation frame
    for i in range(ngen): 
        ax2.clear() 
        ax2.imshow(M, cmap='Greys')
        ax2.set_title(f'generation number = {i:.2f}')
        cells_per_gen.append(np.sum(M))
        ax1.clear() 
        diff_plot=ax1.imshow(diff_grid)
        cb=plt.colorbar(diff_plot, cax=cax, orientation='vertical')
        cb.outline.set_visible(False)
        ax1.set_title(f'Number of timesteps = {i:.2f}')
        anim1.add_frame()
        M = UpdateGrid_v2(M, grid, diff_grid, initial_conc)
        diff_calc = np.zeros((grid,grid))
        diff_calc = diffusion_update(diff_grid,diff_calc,grid,dt,h)
        diff_grid = diff_calc
        anim2.add_frame()

    anim1.show()
    anim2.show()  
    ax3.plot(t_obs,cells_per_gen,'b-',label='Living cells')
    plt.show()
    ax2.set_label("Living cells per generation")

life_v2(50,100,10)
