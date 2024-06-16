# for plotting trajectories from simulations
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


Nplot = 200
p = 3
v1 = int(sys.argv[1])
v2 = int(sys.argv[2])
v3 = int(sys.argv[3])
r = int(sys.argv[4])

	
x = np.loadtxt("trajxt{0}_{1}_{2}_{3}_{4}.txt".format(p, v1, v2, v3, r), delimiter=",")
x = np.transpose(x)

yfile = './trajyt{0}_{1}_{2}_{3}_{4}.txt'.format(p, v1, v2, v3, r)


if (os.stat(yfile).st_size != 0):
	
	y = np.loadtxt("trajyt{0}_{1}_{2}_{3}_{4}.txt".format(p, v1, v2, v3, r), delimiter=",")
	y = np.transpose(y)

	plt.figure(1)
	plt.subplot(1, 2, 1)
	for i in range(Nplot):
		plt.plot(x[0], x[i+1])
	#plt.tick_params(axis = x, pad = 2000000)
	#plt.ylim([-1,10])
	plt.subplot(1, 2, 2)
	for i in range(Nplot):
		plt.plot(y[0], y[i+1])
	
else:
	plt.figure(1)
	for i in range(Nplot):
		plt.plot(x[0], x[i+1])
	#plt.ylim([-1, 1])
	#plt.xlim([0, 1.2])
		
plt.savefig("trajplot{0}_{1}_{2}_{3}_{4}.pdf".format(p, v1, v2, v3, r))
	
plt.show()
