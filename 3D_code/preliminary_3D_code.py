import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from figure_gif2 import figure_gif
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math


def diffusion_update(diff_grid, dt, h, initial_conc, grid, M):
    D = 3.51 * 10 ** (-10)

    D_biofilm = D/2

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
    back[0, :, :] = initial_conc #Dirichlet boundary condition

    diff_grid = ((D * dt * (top + left + bottom + right + front + back - 6 * diff_grid)) / h ** 2) + diff_grid
    diff_biofilm = ((D_biofilm * dt * (top + left + bottom + right + front + back - 6 * diff_grid)) / h ** 2) + diff_grid

    bool_biofilm = M > 0
    diff_grid[bool_biofilm] = diff_biofilm[bool_biofilm]

    return diff_grid


def lactate_consumption(M, diff_grid, net_cells):
    y = 19.1*10**-3  # grams of dry cell per mmol of lactate, from https://onlinelibrary.wiley.com/doi/epdf/10.1002/bit.21101
    mass = 0.466 * 10 ** (-12)  # Derived from weight of an E.coli cell
    mmol_total = (net_cells*mass)/y # mmoles consumed in total by biomass
    mmol_consumption = (mmol_total/np.sum(M)) #mmoles of lactate consumed per cell

    k, i, j = np.where(M == 1) #Returns coordinates of cells that are alive
    diff_grid[k, i, j] = ((diff_grid[k, i, j]*5*(10**(-6))) - mmol_consumption)/5*(10**(-6))

    return diff_grid


def monod(diff_grid,time_growth_update, M):

    #Calculate growth rate
    u_max = 0.827 #maximum growth rate
    Ks = 13.2 #lactate_half_saturation constant (mM)
    S = (np.average(diff_grid)) #average lactate_concentration (mM)?
    u = u_max*(S/(S + Ks)) #specific growth rate

    initial_cells = np.sum(M) #Total number of cells at the start
    final_cells = (1 + u)*np.sum(M) #Total number of cells after growth update
    initial_t = math.log10(initial_cells)/math.log10(1+u)
    final_t = math.log10(final_cells)/math.log10(1+u)

    t = np.linspace(int(initial_t), int(final_t), int(((3600/time_growth_update)+1)))

    obs = []
    total_before = initial_cells
    for i in t:
        total_after = (1+u)**i
        net = total_after-total_before
        obs.append(net)
        total_before = total_after

    return obs


def UpdateGrid_v5_3d(M, grid, time_growth_update, b):
    ###This function will still have the problem with the torus conditions at the borders

    can_divide=[] #coordinates of cells that can divide
    empty_cells = []

    i_start = 0
    i_end = 0
    j_start = 0
    j_end = 0
    k_start = 0
    k_end = 0

    # Generate list of cells that can divide
    for k in range(grid):
        for i in range(grid):
            for j in range(grid):
                coordinates = []
                if M[k, i, j] == 1:
                    if (i > 0) and (i < grid - 1) and (j > 0) and (j < grid - 1) and (k > 0) and (k < grid -1):
                        i_start = i - 1
                        i_end = i + 2
                        j_start = j - 1
                        j_end = j + 2
                        k_start = k - 1
                        k_end = k + 2
                        step = 2
                    else:
                        step = 1
                        if i == 0:
                            i_start = i
                        if i == grid-1:
                            i_end = i + 1
                        if j == 0:
                            j_start = j
                        if j == grid - 1:
                            j_end = j + 1
                        if k == 0:
                            k_start = k
                        if k == grid - 1:
                            k_end = k + 1

                    for x in range(i_start, i_end, step):
                        if x != i and M[k,x,j,] == 0:
                            coordinates.append([k,x,j])

                    for y in range(j_start, j_end, step):
                        if y != i and M[k,i,y] == 0:
                            coordinates.append([k,i,y])

                    for z in range(k_start, k_end, step):
                        if z != k and M[z,i,j] == 0:
                            coordinates.append([z,i,j])
                    #print("The values of x and y are:",x,y)

                    if len(coordinates) >= 1:
                        empty_cells.append(coordinates)
                        can_divide.append([k,i,j])


    divisions = np.random.normal(b,0.1,1)
    divisions = int(round(divisions[0]))

    if divisions == -1:
        divisions = 0

    if divisions > len(can_divide):
        divisions = len(can_divide)

    for value in range(divisions):
        if len(can_divide) == 0:
            break
        divider = random.choice(can_divide)
        index = int(can_divide.index(divider))
        can_divide.remove(divider)

        # Choosing one square at random and dividing
        choice = random.choice(empty_cells[index])
        i = choice[1]
        j = choice[2]
        k = choice[0]
        M[k, i, j] = 1

    #To be added: cells growing towards electrode surface
    return M


def MaxCurrent(grid, M, L):
    '''
    Calculates current output after each cell growth update
    '''
    n = 12  # Number of electrons produced per lactate
    F = 96485  # Faraday's constant
    #L = 5 * 10 ** (-6)  # Biofilm thickness
    u_max = 0.827 / 3600  # units are s^-1
    Y = 19.1 #grams of dry cell per mol of lactate

    y_s = n * F  # unit = C/mol of lactate
    m = 0.466 * 10 ** (-12)  # mass of a cell in dry cell weight according to bionumbers = 0.28

    # Calculate area using total number of live cells
    # Assumes cells are all clustered together
    A = (grid * 5 * (10 ** (-6))) ** 2 #Area isn't calculated correctly
    V = A * grid* 5 * (10 ** (-6)) #The whole volume being modelled
    X = (m * np.sum(M[0])) / V
    q_max = u_max/Y  # mol lactate per g of biomass

    j_max = y_s * q_max * X * L

    return j_max

def DirectTransfer(j_max):
    F = 96485  # Faraday's constant
    R = 8.31
    T = 273 #Kelvin
    E = 0.600  #mV #From conduction biofilm paper (Marcus et al.)
    E_ka = -0.155 # mV

    nernst = 1/(1/(1 + np.exp((-F/(R*T))*(E-E_ka))))
    j = j_max*nernst

    return j

def mediator_diffusion(D, diff_grid, dt, h, grid):

    left = np.zeros((grid, grid, grid))
    right = np.zeros((grid, grid, grid))
    top = np.zeros((grid, grid, grid))
    bottom = np.zeros((grid, grid, grid))
    front = np.zeros((grid, grid, grid))
    back = np.zeros((grid, grid, grid))

    top[:, 0, :] = diff_grid[:, -1, :]
    top[:, 1:, :] = diff_grid[:, :-1, :]

    left[:, :, 0] = diff_grid[:, :, -1]
    left[:, :, 1:] = diff_grid[:, :, :-1]

    bottom[:, :-1, :] = diff_grid[:, 1:, :]
    bottom[:, -1, :] = diff_grid[:, 0, :]

    right[:, :, :-1] = diff_grid[:, :, 1:]
    right[:, :, -1] = diff_grid[:, :, 0]

    front[:-1, :, :] = diff_grid[1:, :, :]
    front[-1, :, :] = diff_grid[-2, :, :] #Neumann boundary conditions

    back[1:, :, :] = diff_grid[:-1, :, :]
    back[0, :, :] = 0 #diff_grid[-1, :, :] #Dirichlet boundary condition

    #print("Top:\n{0} \n Left: \n{1} \n Bottom: \n{2} \n Right: \n{3} \n Center:\n {4}".format(top, left, bottom, right, diff_grid))

    diff_grid = ((D * dt * (top + left + bottom + right + front + back - 6 * diff_grid)) / h ** 2) + diff_grid

    return diff_grid


def MediatorTransfer(red_grid, ox_grid, grid, D_med, M, h):
    k = 1.6*10**(-5)
    a = 0.5

    n_e = 2
    R = 8.31
    T = 273
    F = 96485

    total, n, thickness = BiofilmThickness(M, h)

    z = np.reshape(total, (1, grid ** 2)) #index of z axis
    z = z.astype(int) #changes type to integer
    col = [np.arange(grid)] * grid #create grid with numbers for colum indexing
    row = np.transpose(col) #create grid with numbers for colum indexing
    row = np.reshape(row, (1, grid ** 2))
    col = np.reshape(col, (1, grid ** 2))

    bottom = ox_grid[0]
    top = ox_grid[z, row, col]
    top = np.reshape(top, np.shape(bottom))

    gradient = top - bottom

    flux = (D_med*gradient)/thickness

    j_M = F*n_e*flux
    j = np.sum(j_M)
    #print(flux)
    print("Total current:", j)
    #print("Flux:", flux)
    #print("Biofilm thickness: {0}\n\n".format(thickness))

    return j

def BiofilmThickness(M, h):
    # Count number of cells per column
    total = M.sum(axis = 0)
    n = total > 0 #gives cells with thickness greater than zero
    total[n] -= 1 #gives index in z-axis for each stack of cells

    thickness = total * h #individual thickness of grid
    return total, n, thickness

def MediatorSecretion(ox_grid, M):
    Y = 0.0029 #mM of flavin per cell (arbritrary)

    #Set a value for the rate of mediator production
    secretion = Y*M #This may need to be agent-based

    #Increase mediator concentration in each cell by that amount (oxidized concentration)
    ox_grid += secretion

    return ox_grid


def life_v5(grid, n_updates_growth, time_growth_update):
    # Setting up initial lactate diffusion conditions and plotting
    # Current intial concentration is arbitrary, we should probably ensure we maintain critical flow rate / concentration. Where to find that?
    # https://www.sciencedirect.com/science/article/pii/S0960852412007985?via%3Dihub

    initial_conc = 10 #mM
    maximum_growth_rate = 0.827
    dt = 0.01  # can't be bigger than 0.01780626781 secs for h = 5*(10**(-6))
    h = 5.0 * (10 ** (-6))

    diff_grid = np.zeros((grid, grid, grid))  # grid with lactate concentrations
    diff_calc = np.zeros((grid, grid, grid))  # grid that stores changes

    #Initial lactate conditions
    diff_grid[:] = initial_conc
    diff_calc[:] = initial_conc

    #Initialize grids for mediator diffusion
    med_conc = 100  # mM
    ox_grid = np.zeros((grid, grid, grid))
    ox_grid[-10:] = med_conc
    red_grid = np.zeros((grid, grid, grid))

    dt_med = 0.009
    D_med = 0.43 * 10**(-9) #m^2 per second
    D_med_biofilm = D_med/2


    # Initialize agent and number of neighbours
    M = np.zeros((grid, grid, grid))
    midpoint = round(grid/2)
    #M[(midpoint - 5):(midpoint + 5), (midpoint - 5):(midpoint + 5), (midpoint - 5):(midpoint + 5)] = 1

    M[0:2, :] = 1

    print("initial cells", np.sum(M))

    # anim2=figure_gif(figsize=(15,7.5))
    # ax2=anim2.add_subplot(1,1,1)
    # ax2.imshow(M, cmap='Greys')
    # ax2.set_title('generation number = 0')

    fig1 = plt.figure(figsize=(12, 8))
    ax1 = fig1.add_subplot(1, 1, 1)

    minute = 60

    D = 3.51 * 10 ** (-10)
    D_biofilm = D / 2

    for i in range(round(minute / dt)):
        diff_grid = diffusion_update(diff_grid, dt, h, initial_conc, grid, M)
        #ox_grid = MediatorSecretion(ox_grid, M)
        #ox_grid = mediator_diffusion(D_med_biofilm, ox_grid, dt, h, grid)
        # red_grid = diffusion_update(D_med, red_grid, dt, h, med_conc)

    cells_per_update = [0]
    j_max = [0]
    DT = [0]
    MT = [0]
    thickness = [0]
    I = [0]

    lactate_conc = [np.sum(diff_grid)]
    t_obs = [0]
    t = 0

    hours = round((n_updates_growth * time_growth_update) / 3600)
    if hours == 0:
        hours = 1
    timestep = 0
    for a in range(hours):
        obs = monod(diff_grid, time_growth_update, M)
        for b in obs:
            # ax2.clear()
            # ax2.imshow(M, cmap='Greys')
            # ax2.set_title(f'generation number = {(t):.2f}')
            # anim2.add_frame()
            cells_per_update.append(np.sum(M))
            t_obs.append(t)
            net_cells = cells_per_update[-1] - cells_per_update[-2]
            diff_grid = lactate_consumption(M, diff_grid, net_cells)
            M = UpdateGrid_v5_3d(M, grid, time_growth_update, b)

            # total, n, L = BiofilmThickness(M, h)
            # thickness.append(L)

            # j = MaxCurrent(grid, M, L)
            # j_max.append(j)
            # I_0 = DirectTransfer(j)
            # DT.append(I_0)
            #
            # m_current = MediatorTransfer(red_grid, ox_grid, grid, D_med_biofilm, M, h)
            #
            # MT.append(m_current)
            # total = I_0 + m_current
            # I.append(total)

            total_lactate = np.sum(diff_grid)
            lactate_conc.append(total_lactate)

            t += 1
            for c in range(round(time_growth_update / dt)):
                diff_grid = diffusion_update(diff_grid, dt, h, initial_conc, grid, M)
                # ox_grid = mediator_diffusion(D_med_biofilm, ox_grid, dt, h, grid)
                # red_grid = mediator_diffusion(D_med_biofilm, red_grid, dt, h, grid)

    print("Total time in hours =", (n_updates_growth * time_growth_update) / 3600)
    print("Total cells, final =", np.sum(M))

    cells_per_update.pop(0)
    t_obs.pop(0)
    # I.pop(0)
    lactate_conc.pop(0)
    # j_max.pop(0)
    # DT.pop(0)
    #MT.pop(0)
    #thickness.pop(0)


    ax1.plot(t_obs, cells_per_update, 'b-', label='Living cells')
    plt.xlabel("Number of generations")
    plt.ylabel("Number of living cells")
    plt.title("Living cells")
    plt.show()

    plt.figure()
    plt.plot(t_obs, lactate_conc, 'b-')
    plt.xlabel("Number of generations")
    plt.ylabel("Lactate concentration")
    plt.title("Substrate concentration")
    plt.show()

    # plt.figure()
    # plt.plot(t_obs, I, 'b-', label='Total Current output')
    # plt.xlabel("Number of generations")
    # plt.ylabel("Total current density (A/m2)")
    # plt.title("Total Current output")
    # plt.show()

    # plt.figure()
    # plt.plot(t_obs, j_max, 'b-', label='Maximum current density')
    # plt.xlabel("Number of generations")
    # plt.ylabel("Maximum current density (A/m2)")
    # plt.title("Maximum current density")
    # plt.show()
    #
    # plt.figure()
    # plt.plot(t_obs, DT, 'b-', label='Maximum current density')
    # plt.xlabel("Number of generations")
    # plt.ylabel("Direct Transfer current density (A/m2)")
    # plt.title("Direct Transfer current density")
    # plt.show()
    #
    # plt.figure()
    # plt.plot(t_obs, thickness, 'b-')
    # plt.xlabel("Number of generations")
    # plt.ylabel("Average biofilm thickness (m)")
    # plt.title("Biofilm Thickness")
    # plt.show()

    # plt.figure()
    # plt.plot(t_obs, MT, 'b-')
    # plt.xlabel("Number of generations")
    # plt.ylabel("Mediator Transfer current density (A/m2)")
    # plt.title("Mediator Transfer current density")
    # plt.show()


def main():
    life_v5(20, 144, 200)


if __name__ == "__main__":
    main()
