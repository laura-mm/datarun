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





col = ['m', 'r', 'y', 'g', 'c', 'b']
tit = ['$\gamma_2 = -1$', '$\gamma_2 = -0.8$', '$\gamma_2 = -0.5$', '$\gamma_2 = 0$', '$\gamma_2 = 0.5$', '$\gamma_2 = 1$']
#tit = ['$\sigma_2 = \sigma_3 = 10^{-0.9}$', '$\sigma_2 = \sigma_3 = 10^{-0.6}$', '$\sigma_2 = \sigma_3 = 10^{-.03}$']
lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$', '$(g)$', '$(h)$', '$(i)$', '$(j)$', '$(k)$', '$(l)$']
muval = np.linspace(-4, 2, 25)


# sigma, z1, z2, top, middle, bottom, M, q, help for g
# sigma, z1, phi, M, q, help for measfp

#binn = 80



mi = 24

plt.figure(figsize=(12,8))
#plt.title("$\mu_2 = {0}$".format(muval[mi]))

for gi in range (6):
	plt.subplot(2, 3, gi + 1)

	fp = np.loadtxt("hist2fp_{0}_{1}.txt".format(mi, gi), delimiter = ",")
	fp = fp.reshape(21, 3)

	for si in range (21):
		gaus_m = fp[si][1]
		gaus_sd = fp[si][0]
		#x1 = (g[1][i[si]]*gaus_sd) + gaus_m
		x = np.linspace(0, 10, 10001)
		
	#data = np.loadtxt("hist23_1_{0}_{1}.txt".format(gi, si), delimiter = ",")
	#print(data.shape)
		
		
	#plt.hist(data, binn, density = 1, edgecolor = col[gi], color = col[gi])
		
	#plt.hist(data[gi][si], binn, density = 1, edgecolor = col[gi], color = col[gi])
		if (fp[si][2] < 50):
			plt.plot(x, np.exp((-(x-gaus_m)**2)/(2*gaus_sd**2))/(gaus_sd*np.sqrt(2*math.pi)), label = "$z_1={0}, m={1}".format(fp[si][2], fp[si][1]), color = col[gi]) #"$z_1 = {0}$".format(fp[si][2]), color = col[gi])
		else:
			plt.plot(x, np.exp((-(x-gaus_m)**2)/(2*gaus_sd**2))/(gaus_sd*np.sqrt(2*math.pi)), color = col[gi])
		plt.xlim(0, 10)
		plt.ylim(0, 1)
		plt.legend(loc='upper right', ncol=1, borderaxespad=0, fontsize = 5)

#		if gi == 2:
#			plt.xlabel(r'$x^*$', fontsize = 15)
#		if si == 0:
#			plt.ylabel(r'$p(x^*)$', fontsize = 15)
		#plt.xticks(np.arange(0, 4, 1))
		#plt.yticks(np.arange(0, 5, 1))
#		plt.text(0.8, 8, lett[3*gi + si], fontsize = 15)
#		if gi == 0:
	plt.title(tit[gi], fontsize = 15)
		
		
	#if gi == 2:
	#plt.plot(0, 0, col[0], label = lab[0])
	#plt.plot(0, 0, col[1], label = lab[1])
	#plt.plot(0, 0, col[2], label = lab[2])
	
	#if (si == 0 & gi == 2):
#plt.subplots_adjust(left = 0.07, right = 0.98, bottom = 0.15, top = 0.94, wspace = 0.2, hspace = 0.3)
plt.savefig('hist2_{0}.pdf'.format(mi))


plt.show()

