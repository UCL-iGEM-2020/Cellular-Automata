# Desalination rate
import scipy.integrate
import numpy as np  # not needed because will be in main code
import matplotlib.pyplot as plt  # also not needed because will be in main code

E_0 = 0.219
sep = 0.014  # specific energy production kwh/m3
t = 168  # typical time needed to desalinate in our conditions (h)


plotting_data = np.loadtxt('/Users/Oliver_Leblanc/Desktop/RUNS/RUN2/1week4.txt')

t_obs =  plotting_data[0,:]
I = plotting_data[3,:]


t_obs = np.delete(t_obs,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
I = np.delete(I,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])

I_fit = np.polyfit(t_obs, I, 2)  # fitting the points to a polynomial function
I_t = np.poly1d(I_fit)  # getting the object representing the polynomial
I_t = lambda x:-1.935*(10**(-8))*(x**2) - 3.349*(10**(-6))*x + 0.04029
t_lim = np.linspace(0, 168, 169)
I_int = scipy.integrate.quad(I_t,0,t)
print(I_int)

# final put together of the equation:
sep = 0.014  # specific energy production kwh/m3
t = 168  # typical time needed to desalinate in our conditions (h)
I_int = 6.6908753856
ndr = E_0/sep*I_int/t

print(ndr)

# lets see what it looks like, might need to adjust the polynomial degree
fig1 = plt.figure(figsize=(12,8))
ax1=fig1.add_subplot(2,1,1)
ax2=fig1.add_subplot(2,1,2)
ax1.set_xlabel('Time in hours')
ax1.set_ylabel('Current density')
ax2.set_xlabel('Time in hours')
ax2.set_ylabel('Current density')
ax1.plot(t_obs,I, "#29cce6", label='Actual')
ax2.plot(t_lim,I_t(t_lim), "#29cce6", label='fit')
plt.show()
