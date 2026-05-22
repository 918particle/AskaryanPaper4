import numpy as np
import matplotlib.pyplot as plt

#constants
f_0 = 0.19 #GHz
gamma = 0.02 #GHz
E_0 = 1 #No units
R_0 = 1 #No units
omega_C = 2*np.pi*1.0#GHz
omega_0 = 2*np.pi*0.6 #GHz
epsilon = omega_0/omega_C #No units

#Derived constants
z_1 = complex(2*omega_C-2*np.pi*gamma,2*np.pi*f_0)
z_2 = complex(omega_0-2*np.pi*gamma,2*np.pi*f_0)
z_3 = complex(-omega_0-2*np.pi*gamma,2*np.pi*f_0)

#time
t_min = -10 #Nanoseconds
t_max = 100 #Nanoseconds
dt = 0.2 #Nanoseconds
times = np.arange(t_min,t_max,dt)

#signal
signal = np.zeros(times.shape)
i=0
for t in times:
	if(t>=0):
		c1 = np.exp(z_1*t)/z_1 - 1/z_1
		c2 = np.exp(z_2*t)/z_2 - 1/z_2
		c3 = np.exp(z_3*t)/z_3
		signal[i]=E_0*R_0*(2*np.exp(-2*omega_C*t)*c1.real-
			(1+0.5*epsilon)*np.exp(-omega_0*t)*c2.real-
			(1-0.5*epsilon)*np.exp(omega_0*t)*c3.real)
	else:
		signal[i]=E_0*R_0*(1-0.5*epsilon)*(-1/z_3).real*np.exp(omega_0*t)
	i+=1

#graph
plt.plot(times,signal)
plt.xlim(t_min,t_max)
plt.ylim(-0.25,0.25)
plt.xlabel("Time (ns)")
plt.ylabel("Signal (V)")
plt.savefig('test.png')