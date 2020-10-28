import random
import numpy as np
import math
import matplotlib.pyplot as plt
import sys 
import os

def cells_per_agent(layers):
    D = 1.855 * 10 ** (-9) # m^2 / sec https://link.springer.com/article/10.1007/s10953-005-6987-3
    h = layers*5*(10**(-6))
    n_cells = layers**3
    dt = (h**2)/(6*D) - 0.1 * ((h**2)/(6*D))

    print("Cells per agent: ", n_cells)
    print("timestep: ", dt)

    return h, n_cells, dt


def lactate_input(bulk_mmol,boundary_cubes,bulk_cubes,production_per_timestep,stirred_cubes,diff_grid):
    
    # Check lactate at the boundary
    boundary_lactate=np.average(diff_grid[0, :, :])
    
    #Compute total lactate + added lactate and distribute throughout grid
    bulk_mmol=(((boundary_lactate*boundary_cubes)+(bulk_mmol*bulk_cubes)+production_per_timestep)/stirred_cubes)
    
    
    return bulk_mmol

def diffusion_update(diff_grid, dt, h, bulk_mmol, grid, M):
    
    D = 1.855 * 10 ** (-9) # m^2 / sec https://link.springer.com/article/10.1007/s10953-005-6987-3

    D_biofilm = D/2 # m^2 / sec
    
    left = np.zeros((grid, grid, grid))
    right = np.zeros((grid, grid, grid))
    top = np.zeros((grid, grid, grid))
    bottom = np.zeros((grid, grid, grid))
    front = np.zeros((grid, grid, grid))
    back = np.zeros((grid, grid, grid))

    #Toroidal boundary conditiions
    top[:, 0, :] = diff_grid[:, -1, :]
    top[:, 1:, :] = diff_grid[:, :-1, :]

    left[:, :, 0] = diff_grid[:, :, -1]
    left[:, :, 1:] = diff_grid[:, :, :-1]

    bottom[:, :-1, :] = diff_grid[:, 1:, :]
    bottom[:, -1, :] = diff_grid[:, 0, :]

    right[:, :, :-1] = diff_grid[:, :, 1:]
    right[:, :, -1] = diff_grid[:, :, 0]

    front[:-1, :, :] = diff_grid[1:, :, :]
    front[-1, :, :] = diff_grid[-2, :, :] #Neumann boundary condition

    back[1:, :, :] = diff_grid[:-1, :, :]
    back[0, :, :] = bulk_mmol # Set up concentration gradient with bulk liquid

    diff_grid = ((D * dt * (top + left + bottom + right + front + back - 6 * diff_grid)) / h ** 2) + diff_grid
    diff_biofilm = ((D_biofilm * dt * (top + left + bottom + right + front + back - 6 * diff_grid)) / h ** 2) + diff_grid

    bool_biofilm = M > 0
    diff_grid[bool_biofilm] = diff_biofilm[bool_biofilm]

    return diff_grid

def mediator_diffusion(grid, med_grid, dt, h, M, bc):
    '''
    Function similar to the lactate diffusion function with a different boundary condition.
    To be used to update the oxidized and reduced mediator concentrations.
    The mediator is flavin mononucleotide (FMN)

    Concentrations are in moles
    '''
    D_med = 0.43 * 10 ** (-9)  # m^2 per second
    D_med_biofilm = D_med/2

    left = np.zeros((grid, grid, grid))
    right = np.zeros((grid, grid, grid))
    top = np.zeros((grid, grid, grid))
    bottom = np.zeros((grid, grid, grid))
    back = np.zeros((grid, grid, grid))
    front = np.zeros((grid, grid, grid))

    top[:, 0, :] = med_grid[:, -1, :]
    top[:, 1:, :] = med_grid[:, :-1, :]

    left[:, :, 0] = med_grid[:, :, -1]
    left[:, :, 1:] = med_grid[:, :, :-1]

    bottom[:, :-1, :] = med_grid[:, 1:, :]
    bottom[:, -1, :] = med_grid[:, 0, :]

    right[:, :, :-1] = med_grid[:, :, 1:]
    right[:, :, -1] = med_grid[:, :, 0]


    front[:-1, :, :] = med_grid[1:, :, :]
    front[-1, :, :] = med_grid[-2, :, :] #+ 2*h*(r*h**2)/D_med_biofilm #2*(h*r)/D_med_biofilm #BV flux


    back[1:, :, :] = med_grid[:-1, :, :]
    back[0, :, :] = bc

    med_biofilm = ((D_med_biofilm * dt * (top + left + bottom + right + front + back - 6 * med_grid)) / h ** 2) + med_grid

    med_grid = ((D_med * dt * (top + left + bottom + right + front + back - 6 * med_grid)) / h ** 2) + med_grid

    # Update the grid to include concentrations for biofilm diffusion if there are cells present
    bool_biofilm = M > 0
    med_grid[bool_biofilm] = med_biofilm[bool_biofilm]

    return med_grid


def compute_ratio(M,grid):
 # computing the ratio of cells that are active to cells that can divide
    
    left = np.zeros((grid, grid, grid))
    right = np.zeros((grid, grid, grid))
    top = np.zeros((grid, grid, grid))
    bottom = np.zeros((grid, grid, grid))
    front = np.zeros((grid, grid, grid))
    back = np.zeros((grid, grid, grid))

    #Toroidal boundary conditiions
    top[:, 0, :] = M[:, -1, :]
    top[:, 1:, :] = M[:, :-1, :]

    left[:, :, 0] = M[:, :, -1]
    left[:, :, 1:] = M[:, :, :-1]

    bottom[:, :-1, :] = M[:, 1:, :]
    bottom[:, -1, :] = M[:, 0, :]

    right[:, :, :-1] = M[:, :, 1:]
    right[:, :, -1] = M[:, :, 0]

    front[:-1, :, :] = M[1:, :, :]
    front[-1, :, :] = M[-2, :, :] #Neumann boundary condition

    back[1:, :, :] = M[:-1, :, :]
    back[0, :, :] = M[1, :, :] #Neumann boundary condition

    bool_top= top > 0
    bool_bottom= bottom > 0
    bool_right= right > 0
    bool_left= left > 0
    bool_front= front > 0
    bool_back= back > 0

    top[bool_top] = 1
    bottom[bool_bottom] = 1
    right[bool_right] = 1
    left[bool_left] = 1
    front[bool_front] = 1 
    back[bool_back] = 1 

    total = top + bottom + right + left + front + back
    bool_active = M != 1
    total[bool_active] = 10
    bool_total = total < 6
    can = np.count_nonzero(bool_total == True)
    
    if can == 0:
        ratio_active_can=0
    else:
        ratio_active_can=(np.count_nonzero(M == 1))/can

    return ratio_active_can

def cell_behaviour(M,M2,diff_grid,time_growth_update,grid,h,ratio_active_can, ox_grid,n_cells):

    maximum_growth_rate=0.827
    lactate_half_saturation = 13.2 #mM
    mediator_half_saturation = 0.0001 #mM
    mmol_per_agent=(0.00000000000049/(19.1*10**-3))*n_cells #https://onlinelibrary.wiley.com/doi/epdf/10.1002/bit.21101
    maximum_mmol=mmol_per_agent*((time_growth_update/3600)*maximum_growth_rate)
    n=(3600/time_growth_update)
    uptake_grid = np.zeros((grid,grid,grid))
    nernst = CalculateNernst()

    for k in range(grid):
            for i in range(grid):
                for j in range(grid):
                    if M[k,i,j] == 1 or M[k,i,j] == 2:

                        # Computing growth rate
                        lactate_concentration=(diff_grid[k,i,j])/((h**3)*1000) #mM
                        mediator_conc =(ox_grid[k,i,j])/((h**3)) #mM

                        if k == grid-1:
                            growth_rate=maximum_growth_rate*(lactate_concentration/(lactate_concentration+lactate_half_saturation))*nernst

                        else:
                            growth_rate=maximum_growth_rate*(lactate_concentration/(lactate_concentration+lactate_half_saturation))*(mediator_conc/(mediator_conc+mediator_half_saturation))


                        # Preventing compound growth
                        prevent_compound=((1+growth_rate/n)**n)/(1+growth_rate)

                        # Computing probability of division
                        deterministic_divisions=(growth_rate*(time_growth_update/3600))
                        divisions=((growth_rate*(time_growth_update/3600))*(ratio_active_can/prevent_compound))
                        
                        # Computing lactate consumption
                        specific_mmol=mmol_per_agent*deterministic_divisions*n_cells

                        #Store specific lactate consumption in a grid so it can be used later
                        uptake_grid[k, i, j] = specific_mmol

                        # Cells compute their chance of becoming quiescent
                        prob_quiescence=math.e**(-(diff_grid[k,i,j]*(1/maximum_mmol)))
                        
                        if M[k,i,j] == 1:
                            if random.random() <= prob_quiescence:
                                M2[k,i,j] = 2
                                #print("Quiescence!")
                                continue

                        # If cells are quiescent, check if there's enough lactate for them to go back to normal
                        if M[k,i,j] == 2:
                            if growth_rate >= (0.827*0.1):
                                #print("Recovery!")
                                M2[k,i,j] = 1

                        # If cell stays quiescent, there is a probability of it dying this timestep
                            else:
                                if random.random() <= (time_growth_update/(3600*24)):  # Arbitrary, Pedro has to figure out
                                    #print("Death!")
                                    M2[k,i,j] = 3
                                    continue
                                else:
                                    continue

                        # Notice "continue" statements above as cell will not consume lactate or divide if quiescent
                        diff_grid[k,i,j] -= np.random.normal(specific_mmol,specific_mmol/4,1)

                        # Making a list of empty squares around the divider
                        coordinates=[]
                        
                        if k != (grid-1):
                            if M[k+1,i,j] == 0:
                                coordinates.append(0)
                        
                        if i==(grid-1):
                            if M[k,0,j] == 0:
                                coordinates.append(1)
                        else:
                            if M[k,i+1,j] == 0:
                                coordinates.append(2)
                                
                        
                        if i == 0:
                            if M[k,grid-1,j] == 0:
                                coordinates.append(3)
                        else:
                            if M[k,i-1,j] == 0:
                                coordinates.append(4)
                                
                        
                        if j == (grid-1):
                             if M[k,i,0] == 0:
                                coordinates.append(5)
                        else:
                            if M[k,i,j+1] == 0:
                                coordinates.append(6)
                                
                                
                        if j == 0:
                            if M[k,i,grid-1] == 0:
                                coordinates.append(7)
                        else:
                            if M[k,i,j-1] == 0:
                                coordinates.append(8)                              
                                
                        if k != 0:
                            if M[k-1,i,j] == 0:
                                coordinates.append(9)
                            
                        # Cell has a chance of dividing if it has empty squares around it
                        # If cell is meant to divide more than once, we compute its chances that many times
                        for chances in range(math.ceil(divisions)):
                            if len(coordinates) > 0:
                                singular_division=divisions/math.ceil(divisions)
                                if random.random() <= singular_division:
                                    if coordinates[0] == 0:
                                        coordinates.remove(0)
                                        M2[k+1,i,j] = 1
                                        continue
                                        
                                    if coordinates[0] == 9:
                                        coordinates.remove(9)
                                        M2[k-1,i,j] = 1
                                        continue
                                        
                                    choice=random.choice(coordinates)
                                    # Choosing one square at random and dividing
                                    if choice == 1:
                                        coordinates.remove(1)
                                        M2[k,0,j] = 1
                                    if choice == 2:
                                        coordinates.remove(2)
                                        M2[k,i+1,j] = 1
                                    if choice == 3:
                                        coordinates.remove(3)
                                        M2[k,grid-1,j] = 1                                    
                                    if choice == 4:
                                        coordinates.remove(4)
                                        M2[k,i-1,j] = 1
                                    if choice == 5:
                                        coordinates.remove(5)
                                        M2[k,i,0] = 1
                                    if choice == 6:
                                        coordinates.remove(6)
                                        M2[k,i,j+1] = 1
                                    if choice == 7:
                                        coordinates.remove(7)
                                        M2[k,i,grid-1] = 1
                                    if choice == 8:
                                        coordinates.remove(8)
                                        M2[k,i,j-1] = 1
    # Now copying all of the living cells from M to M2    
    
    bool_M2=M2 == 0
    M2[bool_M2]=M[bool_M2]
                            
    return M2, diff_grid, uptake_grid

def MaxDirectCurrent(grid, M, h, layers):
    '''
    Calculates current output after each cell growth update for direct transfer
    '''
    n = 4  # Number of electrons produced per lactate
    F = 96485  # Faraday's constant C
    L = h  # Biofilm thickness
    u_max = 0.827 / 3600  # units are s^-1
    Y = 19.1 #grams of dry cell per mol of lactate

    y_s = n * F  # unit = C/mol of lactate
    m = 0.49 * 10 ** (-12)  # mass of a cell in dry cell weight according to bionumbers = 0.28 pg

    # Calculate area using total anode surface being modelled
    A = (grid*h)**2
    V = A * L # Volume of monolayer on anode

    #Calculate concentration of active cells
    bool = M[grid-1] == 1
    X = (m * (layers**2) * np.sum(M[grid-1][bool]))/V
    q_max = u_max/Y  # mol lactate per g of biomass

    j_max = y_s * q_max * X * L

    return j_max

def TotalMaxCurrent(grid, M, h, n_cells):
    '''
    Calculates maximum current density in biofilm
    '''
    n = 4  # Number of electrons produced per lactate
    F = 96485  # Faraday's constant
    u_max = 0.827 / 3600  # units are s^-1
    Y = 19.1  # grams of dry cell per mol of lactate

    y_s = n * F  # unit = C/mol of lactate
    m = 0.49 * 10 ** (-12)  # mass of a cell in dry cell weight according to bionumbers = 0.28 pg

    # Volume of each column
    A = h ** 2
    V = A * (grid * h)

    active_cells = np.count_nonzero(M == 1, axis=0)
    total, n, L = BiofilmThickness(M, h)

    X = (m * active_cells * n_cells) / V  # Concentration of cells in each column
    q_max = u_max / Y  # mol lactate per g of biomass

    j_max = y_s * q_max * X * L  #Calculate j_max for each column
    total_j_max = np.sum(j_max)
    #print("Maximum current density", total_j_max)

    return total_j_max


def CalculateNernst():
    F = 96485  # Faraday's constant
    R = 8.31
    T = 273
    E = -0.136  # V #E_OM
    E_ka = -0.155  # V

    nernst = 1/(1 + np.exp((-F/(R*T))*(E-E_ka)))

    return nernst


def MediatorConcentration(M, ox_grid, red_grid, uptake_grid, time_growth_update, diff_grid, h, n_cells):

    # Calculate mediator secretion
    rate = 1.441 #mmol/gDW/h
    m = 0.00000000000049 #g of each cell
    mmol_secreted = (rate / 3600) * m * n_cells # mmol/s
    mediator_secretion = (mmol_secreted/h**3) * 1000 #mM
    lactate_half_saturation = 13.2 #mM
    lactate_concentration = (diff_grid/h**3)*1000 #mM
    conc_secreted = mediator_secretion * (lactate_concentration / (lactate_concentration + lactate_half_saturation))

    #Calculate amount of moles secreted
    mol_secreted = conc_secreted*(h**3)*(10**-6) #mol/s

    bool_mediator = M[:-1] == 1  # only alive cells will reduce the mediators

    ox_grid[:-1][bool_mediator] = ox_grid[:-1][bool_mediator] + (mol_secreted[:-1][bool_mediator] * time_growth_update)

    # Calculate rate of reduction by the biofilm
    f = 1  #fraction of electrons used for mediator transfer
    # uptake grid is in mmol
    R_M = (4 / 2) * f * (uptake_grid * (10 ** -3)) * time_growth_update  # in mol

    #Check if rate is more than concentration available
    check = ox_grid < R_M
    R_M[check] = ox_grid[check]

    negative_check = R_M < 0
    R_M[negative_check] = 0

    ox_grid[:-1][bool_mediator] -= R_M[:-1][bool_mediator]  # mol
    red_grid[:-1][bool_mediator] += R_M[:-1][bool_mediator]


    return ox_grid, red_grid


def BiofilmThickness(M, h):
    '''
    This function calculates the thickness of each stack in the grid to implement an agent-based behaviour
    '''
    # Count number of cells per column
    total = np.count_nonzero(M, axis = 0) # Total number of cells per column

    n = total > 0 #gives cells with thickness greater than zero
    thickness = total * h # individual thickness of grid

    return total, n, thickness

def ChargeTransfer(red_grid, ox_grid, h):
    k = 1.6*10**(-5) #m/s #heterogeneous reaction rate
    a = 0.5

    n_e = 2 #Number of electrons transferred by mediator
    R = 8.31
    T = 303.15 #Kelvin = 30 degrees
    F = 96485

    charge = n_e*F/(R*T)

    #Potentials
    E = 0.3 #0.300 #V polarized electrode potential
    E_0 = -0.219 #V standard redox potential for mediator redox reaction

    #Assume oxidation and reduction occur reversibly at the anode surface
    r_ox = (red_grid[-1]/h**3)*(k)*np.exp((1-a)*charge*(E-E_0)) #rate of oxidation at the electrode surface
    r_red = (ox_grid[-1]/h**3)*(k)*np.exp(-a*charge*(E-E_0))


    r_ox_grid = r_ox - r_red
    r_red_grid = r_red - r_ox

    return r_red_grid, r_ox_grid

def MediatorTransfer(red_grid, ox_grid, grid, D_med, M, h, r_red_grid, r_ox_grid, n_cells):
    n_e = 2 #Number of electrons transferred by mediator
    F = 96485

    total, n, thickness = BiofilmThickness(M, h)

    # For obtaining the index of the top of each stack
    z = np.reshape(total, (1, grid ** 2))  # index of z axis
    z = z.astype(int)  # change data type to integer
    col = [np.arange(grid)] * grid  # create grid with numbers for column indexing
    row = np.transpose(col)  # create grid with numbers for row indexing
    row = np.reshape(row, (1, grid ** 2))
    col = np.reshape(col, (1, grid ** 2))

    # Calculate concentration gradient of reduced form
    bottom = red_grid[-1]/(h**3)  # anode surface
    top = red_grid[-z, row, col]/(h**3)
    top = np.reshape(top, np.shape(bottom))

    gradient = (bottom - top)

    flux = -(D_med * gradient[n]) / thickness[n]
    #print(thickness[n])

    j_M = F*n_e*flux
    j_M_r = F*n_e*(r_ox_grid[n])

    boolean = j_M_r < j_M

    j_M[boolean] = j_M_r[boolean]

    j = np.sum(j_M)

    # print("Flux:", np.average(flux))
    # print("Anode reaction rate:", np.average(r_ox_grid))
    # print("Total mediator current:", j, "\n")
    return j

def ConductiveTransfer(M, h):
    #Constants
    k_bio = 0.5 #*(10**-3)*100 #mS/cm

    #Initialize potentials
    E_OM = -0.136 #V
    E_anode = 0.300 #V

    total, n, L = BiofilmThickness(M,h)
    #No current production if there's an inactive cell
    j = -k_bio*(E_OM - E_anode)/L
    j_total = np.sum(j[0][0])

    return j_total

def plot_3D_CA(M,grid):
    xdata=[]
    ydata=[]
    zdata=[]
    color=[]
        
    for k in range(grid):
        for i in range(grid):
            for j in range(grid):
                if M[k,i,j] == 1:
                    xdata.append(k)
                    ydata.append(i)
                    zdata.append(j)
                    color.append("black")
                if M[k,i,j] == 2:
                    xdata.append(k)
                    ydata.append(i)
                    zdata.append(j)
                    color.append("blue")
                if M[k,i,j] == 3:
                    xdata.append(k)
                    ydata.append(i)
                    zdata.append(j)
                    color.append("red")
    
    
    return xdata, ydata, zdata, color

def life_v7(grid, hours_of_simulation, time_growth_update, layers):
    # PARAMETRES

    # lactate diffusion, retrieved from https://www.sciencedirect.com/science/article/pii/S0960852412007985?via%3Dihub

    h, n_cells, dt = cells_per_agent(layers)    
    
    # lactate input and chamber configuration
    
    lactate_production_coculture = 8.2132 # mmol of lactate per grams of AFDCW per hour
    coculture_steady_state = 0.4 # grams of AFDCW
    anode_surface = 0.01178 # m^2
    anode_volume = 0.09817 # litres
    volume_of_chamber = 250 / 1000 # litres
    volume_of_grid = ((grid*h)**3) * 1000 # litres
    
    n_grids= anode_surface / ((grid*h)**2)
    total_volume = volume_of_chamber - anode_volume # litres
    n_cubes= total_volume / (h**3)
    volume_bulk_liquid = volume_of_chamber - anode_volume - (n_grids*volume_of_grid) # litres
    
    boundary_cubes = grid*grid*n_grids
    bulk_cubes = n_cubes - grid*grid*grid*n_grids
    stirred_cubes = boundary_cubes + bulk_cubes
    
    accumulation_time = 168 * 1 #hours
    accumulated_mmol = lactate_production_coculture*coculture_steady_state*accumulation_time #mmol
    initial_mmol = accumulated_mmol / n_cubes #mmol
    
    total_production = lactate_production_coculture * coculture_steady_state #mmol per hour
    production_per_timestep = total_production / (3600/dt) #mmol per timestep

    # Set up lactate diffusion grid
    
    diff_grid = np.zeros((grid,grid,grid)) #grid with lactate concentrations
    diff_calc = np.zeros((grid,grid,grid)) #grid that stores changes
    
    diff_grid[0:grid,0:grid,0:grid] = initial_mmol
    diff_calc[0:grid,0:grid,0:grid] = initial_mmol

    # Set up mediator diffusion grid

    conc = 1.25*10**-19 #1 mM = 1*10^-3 mol/m^3
    ox_grid = np.zeros((grid, grid, grid))  # mol
    ox_grid[:] = 0
    red_grid = np.zeros((grid, grid, grid))  # mol
    red_grid[:] = conc #From Picoroneai 13 paper in nernst monod fitting paper

    D_med = 0.43 * 10 ** (-9)  # m^2 per second
    D_med_biofilm = D_med / 2

    # Set up agent grid
    M = np.zeros((grid,grid,grid))
    #M[-1, :] = 1
    M[(grid-1),(round(grid/2)),(round(grid/2))] = 1
    
    # Setting up plots and animation
    
    #anim=figure_gif(figsize=(15,7.5))
    #ax2=anim.add_subplot(1,1,1,projection='3d')
    
    #fig1 = plt.figure(figsize=(12,8))
    #ax1=fig1.add_subplot(2,1,1)
    #ax3=fig1.add_subplot(2,1,2)

    #fig2 = plt.figure(figsize=(12,8))
    #ax4=fig2.add_subplot(2,1,1)
    #ax5=fig2.add_subplot(2,1,2)

    #fig3 = plt.figure(figsize=(12,8))
    #ax6=fig3.add_subplot(2,1,1)
    #ax7=fig3.add_subplot(2,1,2)

    #fig4 = plt.figure(figsize=(12,8))
    #ax8 = fig4.add_subplot(2,1,1)
    #ax9 = fig4.add_subplot(2,1,2)

    #fig5 = plt.figure(figsize=(12, 8))
    #ax10 = fig5.add_subplot(1,1,1)

    #fig6 = plt.figure(figsize=(20, 20))
    #ax11 = fig6.add_subplot(3,1,1)


    cells_per_update=[]
    total_lactate=[]
    total_red_conc = []
    total_ox_conc = []
    total_med_conc = []
    t_obs=[]
    thickness = []
    j_max_list = []
    MT = []
    I = []
    MaxI = []
    active_cells = []
    quiescent_cells = []
    dead_cells = []

    t=0

    bulk_mmol = initial_mmol

    counter=1

    r_red_grid, r_ox_grid = ChargeTransfer(red_grid, ox_grid, h)

    
    # Starting simulation for loop
    
    n_updates=(hours_of_simulation*3600)/time_growth_update
    hour_update=time_growth_update/3600

    print("Concentration of reduced form:", np.sum(red_grid))
    print("Total mediator concentration", np.sum(ox_grid) + np.sum(red_grid))
        
    for timesteps in range(math.ceil(n_updates)):
        t += hour_update
        cells_per_update.append(np.count_nonzero(M == 1))
        t_obs.append(t)

        ratio_active_can = compute_ratio(M,grid)
        M2 = np.zeros((grid,grid,grid))
        M2, diff_grid, uptake_grid = cell_behaviour(M,M2,diff_grid,time_growth_update,grid,h,ratio_active_can, ox_grid, n_cells)
        total_lactate.append(np.sum(diff_grid))
        M = M2
        #Current output update

        C_red = red_grid*1000 #converts amount to mmol
        C_ox = ox_grid*1000

        total_red_conc.append(np.sum(C_red))
        total_ox_conc.append(np.sum(C_ox))
        total_med_conc.append(np.sum(C_red)+np.sum(C_ox))
        
        total, n, L = BiofilmThickness(M, h)
        thickness.append(np.average(L))

        j = MaxDirectCurrent(grid, M, h, layers)
        j_max_list.append(j)

        # print("Anode reaction rate:", np.average(r_ox_grid))
        m_current = MediatorTransfer(red_grid, ox_grid, grid, D_med_biofilm, M, h, r_red_grid, r_ox_grid, n_cells)
        MT.append(m_current)

        total = j + m_current
        I.append(total)

        max_current = TotalMaxCurrent(grid, M, h, n_cells)
        MaxI.append(max_current)

        ox_grid, red_grid = MediatorConcentration(M, ox_grid, red_grid, uptake_grid, time_growth_update, diff_grid, h, n_cells)  # moles
        active_cells.append(np.count_nonzero(M == 1))
        quiescent_cells.append(np.count_nonzero(M == 2))
        dead_cells.append(np.count_nonzero(M == 3))


        if t >= counter:
            print("Number of hours computed: ", counter, "/", hours_of_simulation)
            print(diff_grid[1,1,1])
            counter += 1
        for i in range(round(time_growth_update/dt)):
            bulk_mmol = lactate_input(bulk_mmol,boundary_cubes,bulk_cubes,production_per_timestep,stirred_cubes,diff_grid)
            diff_calc = np.zeros((grid,grid,grid))
            diff_calc = diffusion_update(diff_grid, dt, h, bulk_mmol, grid, M)
            diff_grid = diff_calc

            ox_grid = mediator_diffusion(grid, ox_grid, dt, h, M, conc) #mol
            red_grid = mediator_diffusion(grid, red_grid, dt, h, M, 0) #mol

            r_red_grid, r_ox_grid = ChargeTransfer(red_grid, ox_grid, h)

            #Check if rate is more than concentration available
            check = red_grid[-1] < r_ox_grid * dt * h ** 2

            r_ox_grid[check] = red_grid[-1][check] / (dt * h ** 2)
            r_red_grid[check] = -red_grid[-1][check] / (dt * h ** 2)

            negative_check = r_ox_grid < 0
            r_ox_grid[negative_check] = 0
            r_red_grid[negative_check] = 0

            ox_grid[-1] += r_ox_grid * dt * (h ** 2)
            red_grid[-1] += r_red_grid * dt * (h ** 2)


        # print("Concentration of reduced form:", np.sum(red_grid))
        # print("Total mediator concentration", np.sum(ox_grid) + np.sum(red_grid))


    # Plotting results
    
    print("Total time in hours =", hours_of_simulation)
    print("Total agents, final =", np.count_nonzero(M == 1))
    
    t_obs_array = np.array(t_obs)
    cells_array = np.array(cells_per_update)
    lactate_array = np.array(total_lactate)
    total_current_array = np.array(I)
    direct_current_array = np.array(j_max_list)
    mediator_current_array = np.array(MT)
    maximum_current_array = np.array(MaxI)
    biofilm_thickness = np.array(thickness)
    total_ox = np.array(total_ox_conc)
    total_red = np.array(total_red_conc)
    total_med = np.array(total_med_conc)
    a_cells = np.array(active_cells)
    q_cells = np.array(quiescent_cells)
    d_cells = np.array(dead_cells)


    length = len(t_obs)

    plotting_data = np.concatenate((t_obs_array,cells_array,lactate_array,total_current_array,direct_current_array,mediator_current_array,maximum_current_array,biofilm_thickness,total_ox,total_red,total_med,a_cells,q_cells,d_cells)).reshape(14,length)


    np.savetxt(os.path.join("/Users/Oliver_Leblanc/Desktop/RUNS/RUN2", '1week2.txt'), plotting_data)
    np.save(os.path.join("/Users/Oliver_Leblanc/Desktop/RUNS/RUN2", '1week2M'), M)
    np.save(os.path.join("/Users/Oliver_Leblanc/Desktop/RUNS/RUN2", '1week2diff'), diff_grid)
    np.save(os.path.join("/Users/Oliver_Leblanc/Desktop/RUNS/RUN2", '1week2ox'), ox_grid)
    np.save(os.path.join("/Users/Oliver_Leblanc/Desktop/RUNS/RUN2", '1week2red'), red_grid)




    #anim.show()  

    #ax1.set_xlabel('Time in hours')
    #ax1.set_ylabel('cells')
    #ax1.plot(t_obs,cells_per_update,'b-',label='Active cells')

    #ax3.set_xlabel('Time in hours')
    #ax3.set_ylabel('Total Lactate in the grid - mmol')
    #ax3.plot(t_obs,total_lactate,'b-')

    #ax4.plot(t_obs, I, 'b-', label='Total Current output')
    #ax4.set_xlabel("Hours")
    #ax4.set_ylabel("Total current density (A/m2)")
    #ax4.set_title("Total Current output")

    #ax5.plot(t_obs, j_max_list, 'b-', label='Maximum direct current density')
    #ax5.set_xlabel("Hours")
    #ax5.set_ylabel("Maximum direct transfer current density (A/m2)")
    #ax5.set_title("Maximum current density")

    #ax6.plot(t_obs, MT, 'b-', label='Maximum current density')
    #ax6.set_xlabel("Hours")
    #ax6.set_ylabel("Mediator Transfer current density (A/m2)")
    #ax6.set_title("Mediator Transfer current density")

    #ax7.plot(t_obs, MaxI, 'b-')
    #ax7.set_xlabel("Hours")
    #ax7.set_ylabel("Maximum total current density (A/m2)")
    #ax7.set_title("Maximum total current density")
    #plt.show()



life_v7(10,168,10,5)
