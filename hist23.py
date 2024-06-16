# histograms

import numpy as np
import matplotlib.pyplot as plt
import math
#plt.rc('text', usetex=True)
#plt.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')
#plt.rc('xtick',labelsize=15)
#plt.rc('ytick',labelsize=15)
#plt.rc('axes', labelsize=20)
#plt.rc('axes', titlesize=20)





col = ['m', 'g', 'b']
lab = ['$\gamma_2 = -1, \gamma_3 = -1/2$', '$\gamma_2 = \gamma_3 = 0$', '$\gamma_2 = \gamma_3 = 1$']
tit = ['$\sigma_2 = \sigma_3 = 10^{-0.9}$', '$\sigma_2 = \sigma_3 = 10^{-0.6}$', '$\sigma_2 = \sigma_3 = 10^{-.03}$']
lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$', '$(g)$', '$(h)$', '$(i)$', '$(j)$', '$(k)$', '$(l)$',]

# sigma, z1, z2, top, middle, bottom, M, q, help for g
# sigma, z1, phi, M, q, help for measfp

#gival = [0, 3, 5]
#sigmaval = {0.1, 1, 10}


fp = np.loadtxt("hist23fp.txt", delimiter = ",")
fp = fp.reshape(3, 3, len(fp)//9)

binn = 80


plt.figure(figsize=(12,8))

for gi in range (3):

	for si in range (3):

		plt.subplot(3, 3, (3*gi) + 1 + si)
		gaus_m = fp[gi][si][1]
		gaus_sd = fp[gi][si][0]
		#x1 = (g[1][i[si]]*gaus_sd) + gaus_m
		x = np.linspace(0, 10, 10001)
		
		data = np.loadtxt("hist23_1_{0}_{1}.txt".format(gi, si), delimiter = ",")
		print(data.shape)
		
		
		plt.hist(data, binn, density = 1, edgecolor = col[gi], color = col[gi])
		
		#plt.hist(data[gi][si], binn, density = 1, edgecolor = col[gi], color = col[gi])
		plt.plot(x, np.exp((-(x-gaus_m)**2)/(2*gaus_sd**2))/(gaus_sd*np.sqrt(2*math.pi)), 'k')
		plt.xlim(0, 1)
		plt.ylim(0, 10)
		if gi == 2:
			plt.xlabel(r'$x^*$', fontsize = 15)
		if si == 0:
			plt.ylabel(r'$p(x^*)$', fontsize = 15)
		#plt.xticks(np.arange(0, 4, 1))
		#plt.yticks(np.arange(0, 5, 1))
		plt.text(0.8, 8, lett[3*gi + si], fontsize = 15)
		if gi == 0:
			plt.title(tit[si], fontsize = 15)
		
		
	#if gi == 2:
plt.plot(0, 0, col[0], label = lab[0])
plt.plot(0, 0, col[1], label = lab[1])
plt.plot(0, 0, col[2], label = lab[2])
	
	#if (si == 0 & gi == 2):
plt.legend(bbox_to_anchor=(-2.4, -0.6, 3.4, 0.2), loc='lower left', ncol=3, mode="expand", borderaxespad=0, fontsize = 15)
plt.subplots_adjust(left = 0.07, right = 0.98, bottom = 0.15, top = 0.94, wspace = 0.2, hspace = 0.3)
plt.savefig('hist23.pdf')


plt.show()

