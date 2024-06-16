# this used to plot functions from numerical solution

import numpy as np
import matplotlib.pyplot as plt
N = 200000
x = np.loadtxt('abss.txt', delimiter=",")
x = x.reshape(len(x)//2, 2)
x = np.transpose(x)

#cr = np.loadtxt('cr.txt')

plt.figure(1)
plt.plot(x[0], x[1], label = "sig")
#plt.plot(x[0], x[2], label = "grad")
#plt.plot(x[0], x[3], label = "")
#plt.plot(x[0], x[4], label = "gamatot")
plt.ylim([-10, 10])
plt.legend()
plt.axhline(y = 0, linestyle = '--', color = 'r')
plt.axvline(x = 0.0, linestyle = '--', color = 'k')
plt.show()
