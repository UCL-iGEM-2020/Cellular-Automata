import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from figure_gif2 import figure_gif
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

def diffusion_update(grid,M,diff_grid,diff_calc,dt,h):
    D = 3.51*10**(-10)
    
    diff_biofilm = np.zeros((grid,grid))
    D_biofilm = D/2
    
    top = diff_grid[0:-2,1:-1]
    left = diff_grid[1:-1,0:-2]
    bottom = diff_grid[2:,1:-1]
    right = diff_grid[1:-1,2:]
    center = diff_grid[1:-1,1:-1]    
    
    delta_calc = (top + left + bottom + right - 4 * center) / h**2
    calc = diff_grid[1:-1,1:-1]
    diff_calc[1:-1,1:-1] = calc + dt * (D * delta_calc)
    diff_biofilm[1:-1,1:-1] = calc + dt * (D_biofilm * delta_calc)
    
    bool_biofilm=M>0
    diff_calc[bool_biofilm]=diff_biofilm[bool_biofilm]
    
    diff_calc[0,:] = 1*diff_grid[1,:]
    diff_calc[-1,:] = 1*diff_grid[-2,:]
    diff_calc[:,0] = 1*diff_grid[:,1]
    diff_calc[:,-1] = 1*diff_grid[:,-2]
    
    return diff_calc    

def compute_ratio(M,grid):
 # computing the ratio of cells that are active to cells that can divide

    can=0

    for i in range(grid):
            for j in range(grid):
                if M[i,j] == 1:
                    if i==0 and j==0:
                        if M[i+1,j] == 0 or M[grid-1,j] == 0 or M[i,j+1] == 0 or M[i,grid-1] == 0:
                            can+=1
                            continue
                    elif i==grid-1 and j==grid-1:
                        if M[0,j] == 0 or M[i-1,j] == 0 or M[i,0] == 0 or M[i,j-1] == 0:
                            can+=1
                            continue
                    elif i==0 and j==grid-1:
                        if M[i+1,j] == 0 or M[grid-1,j] == 0 or M[i,0] == 0 or M[i,j-1] == 0:
                            can+=1
                            continue
                    elif i==grid-1 and j==0:
                        if M[0,j] == 0 or M[i-1,j] == 0 or M[i,j+1] == 0 or M[i,j-1] == 0:
                            can+=1
                            continue
                    elif i==(grid-1):
                        if M[0,j] == 0 or M[i-1,j] == 0 or M[i,j+1] == 0 or M[i,j-1] == 0:
                            can+=1
                            continue
                    elif i == 0:
                        if M[i+1,j] == 0 or M[grid-1,j] == 0 or M[i,j+1] == 0 or M[i,j-1] == 0:
                            can+=1
                            continue
                    elif j == (grid-1):
                        if M[i+1,j] == 0 or M[i-1,j] == 0 or M[i,0] == 0 or M[i,j-1] == 0:
                            can+=1
                            continue
                    elif j == 0:
                        if M[i+1,j] == 0 or M[i-1,j] == 0 or M[i,j+1] == 0 or M[i,grid-1] == 0:
                            can+=1
                            continue
                    else:
                        if M[i+1,j] == 0 or M[i-1,j] == 0 or M[i,j+1] == 0 or M[i,j-1] == 0:
                            can+=1
                            continue
                    
    ratio_active_can=(np.count_nonzero(M == 1))/can
    
    # Correct multiple singular cell counting

    return ratio_active_can
                    
def cell_behaviour(M,M2,diff_grid,time_growth_update,grid,h,ratio_active_can):

    maximum_growth_rate=0.827
    lactate_half_saturation=13.2 #mM
    mmol_per_cell=(0.00000000000049/(19.1*10**-3)) #https://onlinelibrary.wiley.com/doi/epdf/10.1002/bit.21101
    maximum_mmol=mmol_per_cell*((time_growth_update/3600)*maximum_growth_rate)
    n=(3600/time_growth_update)
    


    for i in range(grid):
            for j in range(grid):
                if M[i, j] == 1 or M[i,j] == 2:
                    
                    # Computing growth rate
                    lactate_concentration=(diff_grid[i,j])/h
                    growth_rate=maximum_growth_rate*(lactate_concentration/(lactate_concentration+lactate_half_saturation))
                                        
                    # Preventing compound growth
                    prevent_compound=((1+growth_rate/n)**n)/(1+growth_rate)

                    # Computing probability of division
                    deterministic_divisions=(growth_rate*(time_growth_update/3600))
                    divisions=((growth_rate*(time_growth_update/3600))*(ratio_active_can/prevent_compound))

                    # Computing lactate consumption
                    specific_mmol=mmol_per_cell*deterministic_divisions
                    
                    # Cells compute their chance of becoming quiescent
                    prob_quiescence=math.e**(-(diff_grid[i,j]*(1/maximum_mmol)))
                    if random.random() <= prob_quiescence:
                        M2[i,j] == 2
                        print("Quiescence!")
                        continue

                    # If cells are quiescent, check if there's enough lactate for them to go back to normal
                    if M[i,j] == 2:
                        if random.random() >= prob_quiescence:
                            print("Recovery!")
                            M2[i,j] == 1

                    # If cell stays quiescent, there is a probability of it dying this timestep
                        else:
                            if random.random() <= (time_growth_update/3600*24):  # Arbitrary, Pedro has to figure out
                                print("Death!")
                                M2[i,j] == 3
                                continue

                    # Notice "continue" statements above as cell will not consume lactate or divide if quiescent
                    diff_grid[i,j] -= np.random.normal(specific_mmol,specific_mmol/4,1)

                    # Making a list of empty squares around the divider
                    coordinates=[]

                    if i==(grid-1):
                        if M[0,j] == 0:
                            coordinates.append(0)
                    else:
                        if M[i+1,j] == 0:
                            coordinates.append(0)
                    if i == 0:
                        if M[grid-1,j] == 0:
                            coordinates.append(1)
                    else:
                        if M[i-1,j] == 0:
                            coordinates.append(1)
                    if j == (grid-1):
                         if M[i,0] == 0:
                            coordinates.append(2)
                    else:
                        if M[i,j+1] == 0:
                            coordinates.append(2)
                    if j == 0:
                        if M[i,grid-1] == 0:
                            coordinates.append(3)
                    else:
                        if M[i,j-1] == 0:
                            coordinates.append(3)


                    # Cell has a chance of dividing if it has empty squares around it
                    # If cell is meant to divide more than once, we compute its chances that many times
                    for chances in range(math.ceil(divisions)):
                        if len(coordinates) > 0:
                            singular_division=divisions/math.ceil(divisions)
                            if random.random() <= singular_division:
                                choice=random.choice(coordinates)
                                # Choosing one square at random and dividing
                                if choice == 0:
                                    coordinates.remove(0)
                                    if i==(grid-1):
                                        M2[0,j] = 1
                                    else: 
                                        M2[i+1,j] = 1

                                if choice == 1:
                                    coordinates.remove(1)
                                    if i == 0:
                                        M2[grid-1,j] = 1
                                    else:
                                        M2[i-1,j] = 1

                                if choice == 2:
                                    coordinates.remove(2)
                                    if j == (grid-1):
                                        M2[i,0] = 1
                                    else: 
                                        M2[i,j+1] = 1
                                if choice == 3:
                                    coordinates.remove(3)
                                    if j == 0:
                                        M2[i,grid-1] = 1
                                    else:
                                        M2[i,j-1] = 1


    # Now copying all of the living cells from M to M2
    bool_M=M>0
    M2[bool_M]=M[bool_M]
                        
    return M2, diff_grid




def life_v6(grid, n_updates_growth, time_growth_update):
    
    # Current intial concentration is arbitrary, we should probably ensure we maintain critical flow rate / concentration. Where to find that?
    # https://www.sciencedirect.com/science/article/pii/S0960852412007985?via%3Dihub
    
    initial_conc=1/(grid*grid)
    dt=0.01 # can't be bigger than 0.01780626781 secs for h = 5*(10**(-6))
    h=5*(10**(-6))

    
    diff_grid = np.zeros((grid,grid)) #grid with lactate concentrations
    diff_calc = np.zeros((grid,grid)) #grid that stores changes
    
    diff_grid[0:grid-1,0:grid-1]=initial_conc
    diff_calc[0:grid-1,0:grid-1]=initial_conc
    
    #Initialize agent
    M = np.zeros((grid,grid))
    M[(round(grid/2)-2):(round(grid/2)+2),(round(grid/2)-2):(round(grid/2)+2)]=1
    
    
    #anim2=figure_gif(figsize=(15,7.5))
    #ax2=anim2.add_subplot(1,1,1)
    #ax2.imshow(M, cmap='Greys')
    #ax2.set_title('generation number = 0')
    
    fig1 = plt.figure(figsize=(12,8))
    ax1=fig1.add_subplot(1,1,1)
    cells_per_update=[0]
    t_obs=[0]
    t=0
    
    for timesteps in range(n_updates_growth):
        #ax2.clear() 
        #ax2.imshow(M, cmap='Greys')
        #ax2.set_title(f'timestep = {(t):.2f}')
        #anim2.add_frame()
        t += 1
        print(np.sum(diff_grid))
        cells_per_update.append(np.count_nonzero(M == 1))
        t_obs.append(t)
        ratio_active_can = compute_ratio(M,grid)
        M2 = np.zeros((grid,grid))
        M2, diff_grid = cell_behaviour(M,M2,diff_grid,time_growth_update,grid,h,ratio_active_can)
        M = M2
        for i in range(round(time_growth_update/dt)):
            diff_calc = np.zeros((grid,grid))
            diff_calc = diffusion_update(grid,M,diff_grid,diff_calc,dt,h)
            diff_grid = diff_calc
    
    print("Total time in hours =",(n_updates_growth*time_growth_update)/3600)
    print("Total cells, final =", np.count_nonzero(M == 1))

    cells_per_update.pop(0)
    t_obs.pop(0)
    
    #anim2.show()  
    ax1.plot(t_obs,cells_per_update,'b-',label='Living cells')
    plt.show()

life_v6(50,250,100)

def expected_growth(initial_cells,n_updates_growth,time_growth_update):
    max_growth_rate=0.827
    hours=(n_updates_growth*time_growth_update)/3600
    expected_cells=initial_cells*(1+max_growth_rate)**hours
    
    return expected_cells    
