# muplots
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

x = np.loadtxt("M3sig0.txt", delimiter=",")
x = np.transpose(x)

muval = np.linspace(-5, 5, 41)

meas = np.loadtxt("allfive3_sig0.txt", delimiter=",")
meas = meas.reshape(41, 3)
meas = np.transpose(meas)

colour = np.loadtxt("allcolour3_sig0.txt", delimiter = ",")
colour = colour.reshape(41, 3)

plt.figure(1)
plt.plot(x[0], x[1], 'r', label = "M1 real")
plt.plot(x[0], x[2], 'r--', label = "M1 imag")
plt.plot(x[0], x[3], 'b', label = "M2 real")
plt.plot(x[0], x[4], 'b--', label = "M2 imag")
for mi in range(41):
	plt.plot(muval, meas[1], 'ko', markersize = 2)#linewidth = 5)#, s = 10)#, edgecolor = 'k') #color = colour[mi])
plt.ylim([-5, 5])
plt.xlabel("$\mu_3$")
plt.ylabel("$M^*$")

plt.legend(bbox_to_anchor=(0.2, -0.3, 0.6, 0.2), loc='lower left', ncol=2, mode="expand", borderaxespad=0)
plt.subplots_adjust(bottom = 0.25, top = 0.96)

plt.savefig("M3sig0.pdf")

plt.show() 

