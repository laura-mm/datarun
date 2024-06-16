#phase23

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math
lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']
code = ['A', 'B', 'C', 'D', 'E', 'G', 'N', 'O']
xlims = [[0, 5], [0, 5], [0, 4], [0, 8], [0, 3], [0, 6], [0, 6], [0, 3]]
ylims = [[0, 5], [0, 10], [0, 5], [0, 10], [0, 4], [0, 3], [0, 4], [0, 8]]
g2 = [-1, -1, -0.8, -0.8, -0.8, 0, -0.88, -0.79]
g3 = [-0.5, -0.5, -0.4, -0.4, -0.4, 0, -0.39, -0.44]
m2 = [-2, -2, -3.5, -4, -3.37, -2, -3.38, -3.38]
m3 = [4, -4, 6, 5, 6, -4, 6, 6]
s2maxval = [5, 5, 4, 8, 3, 6, 6, 3]
s3maxval = [5, 10, 5, 10, 4, 3, 4, 8]


colour = np.loadtxt("allcolour23.txt", delimiter = ",")
colour = colour.reshape(8, 21, 21, 3) #ci s2i s3i col 


plt.subplots(2, 4, figsize = (12, 7.5))

for ci in range(8):
		
	plt.subplot(2, 4, ci + 1)
	
	s2max = s2maxval[ci]
	s2 = np.linspace(0, s2max, 21)
	s3max = s3maxval[ci]
	s3 = np.linspace(0, s3max, 21)

	
	plt.title('$\gamma_2 = {0}$, $\gamma_3 = {1}$ \n $\mu_2 = {2}$, $\mu_3 = {3}$'.format(g2[ci], g3[ci], m2[ci], m3[ci]), fontsize = 15)
	
	for s2i in range(21):
		for s3i in range(21):
			plt.scatter(s2[s2i], s3[s3i], c = colour[ci][s2i][s3i].reshape(1,-1), marker = 'o', edgecolor = 'k', linewidth = 0.5, s = 40)
		
	filec = 'planec{0}.txt'.format(code[ci])
	if (os.path.exists(filec)):
		c = np.loadtxt(filec, delimiter=",")
		c = c.reshape(len(c)//3, 3)
		c = np.transpose(c)
		plt.plot(c[0], c[1], 'm', linewidth = 3)
		
	filem = 'planem{0}.txt'.format(code[ci])
	if (os.path.exists(filem)):
		m = np.loadtxt(filem, delimiter=",")
		m = m.reshape(len(m)//3, 3)
		m = np.transpose(m)
		plt.plot(m[0], m[1], 'y', linewidth = 3)
		
	filea = 'planeda{0}.txt'.format(code[ci])
	if (os.path.exists(filea)):
		a = np.loadtxt(filea, delimiter=",")
		a = a.reshape(len(a)//3, 3)
		a = np.transpose(a)
		plt.plot(a[0], a[1], 'c', linewidth = 3)
		
	fileb = 'planedb{0}.txt'.format(code[ci])
	if (os.path.exists(fileb)):
		b = np.loadtxt(fileb, delimiter=",")
		b = b.reshape(len(b)//3, 3)
		b = np.transpose(b)
		plt.plot(b[0], b[1], 'c', linewidth = 3)
	
	if (ci > 3):	
		plt.xlabel("$\sigma_2$", fontsize = 15)
	if (ci % 4 == 0):
		plt.ylabel("$\sigma_3$", fontsize = 15)
		
	plt.xlim([-s2max/20, s2max*21/20])
	plt.ylim([-s3max/20, s3max*21/20])
	#plt.xlim([0, 10])
	#plt.ylim([0, 10])
	#plt.xticks(fontsize=20)
	#plt.yticks(fontsize=20)
	plt.subplots_adjust(left = 0.05, right = 0.98, top = 0.92, bottom = 0.15, hspace = 0.4)
	
plt.plot(-1, -1, 'm', label = 'linear instability point')
plt.plot(-1, -1, 'c', label = 'divergance point')
plt.plot(-1, -1, 'y', label = 'maximum $\sigma_3$ point')	




plt.legend(bbox_to_anchor=(-3, -0.4, 3.5, 0.2), loc='lower left', ncol=6, mode="expand", borderaxespad=0, fontsize = 15)
			
plt.savefig("phasesim23.pdf")		
#plt.show()

