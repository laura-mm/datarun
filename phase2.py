import numpy as np
import matplotlib.pyplot as plt
import math

col = ['m', 'r', 'y', 'g', 'c', 'b']
col_ = ['m--', 'r--', 'y--', 'g--', 'c--', 'b--']
colo = ['m.', 'r.', 'y.', 'g.', 'c.', 'b.']
lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']

gammalab2 = ['$\gamma_2 = -1$', '$\gamma_2 = -0.8$', '$\gamma_2 = -0.5$', '$\gamma_2 = 0$', '$\gamma_2 = 0.5$', '$\gamma_2 = 1$']

muval = np.linspace(-4, 2, 25)
x = np.linspace(-1.0, 1.0, 21)
sigmaval = np.power(10, x)

# sigma, z1, phi, M, q, help for measfp
# phi, M, q, diversity, dsq, h, for meas



# phase diagram for all gamma without simulation results
######################################################


p4 = [0, 1, 3, 5] # for the paper only these 4

plt.figure(1, figsize = (6, 6))
for gii in range(4):
	
	gi = p4[gii]
	
	z = np.loadtxt("bunin2_{0}.txt".format(gi), delimiter = ",")
	z = z.reshape(len(z)//2, 2)
	z = np.transpose(z)
	plt.semilogy(z[0], z[1], col[gi], linestyle = ':')
	
	m = np.loadtxt("mu2_{0}.txt".format(gi), delimiter = ",")
	m = m.reshape(len(m)//2, 2)
	m = np.transpose(m)
	plt.semilogy(m[0], m[1], col_[gi])
	
	plt.plot(0, 0, col[gi], label = gammalab2[gi])
	
	plt.ylim([10**(-2), 10**2])
	plt.xlim([-4, 3.5])
	plt.xlabel("$\mu_2$", fontsize = 15)
	plt.ylabel("$\sigma_2$", fontsize = 15)
	plt.legend(bbox_to_anchor=(0, -0.3, 1, 0.2), loc='lower left', ncol=2, mode="expand", borderaxespad=0, fontsize = 15)
	plt.subplots_adjust(bottom = 0.25, top = 0.96)
	
plt.savefig("phase2.pdf")

# phase diarams with simulations
###############################################################

colour = np.loadtxt("allcolour2.txt", delimiter = ",")
colour = colour.reshape(25, 6, 21, 3)
colour = np.transpose(colour, (1, 0, 2, 3))

plt.subplots(2, 3, figsize = (12, 8))
for gi in range(6):
	
	plt.subplot(2, 3, gi + 1)
	
	for mi in range(25):
		for si in range(21):
			plt.scatter(muval[mi], sigmaval[si], c=colour[gi][mi][si].reshape(1,-1), marker = 'o', edgecolor = 'k', linewidth = 0.5, s = 40)
	
	#plt.imshow(colour[gi], extent = [-4, 1.5, 0.1, 10], origin = 'lower')
	
	z = np.loadtxt("bunin2_{0}.txt".format(gi), delimiter = ",")
	z = z.reshape(len(z)//2, 2)
	z = np.transpose(z)
	plt.semilogy(z[0], z[1], col[gi], linestyle = ':', linewidth = 3)
	
	m = np.loadtxt("mu2_{0}.txt".format(gi), delimiter = ",")
	m = m.reshape(len(m)//2, 2)
	m = np.transpose(m)
	plt.semilogy(m[0], m[1], col_[gi], linewidth = 3)
	
	plt.ylim([10**(-1.2), 10**1.2])
	plt.xlim([-4.5, 2.5])
	if (gi == 3 or gi == 4 or gi == 5):
		plt.xlabel("$\mu_2$", fontsize = 15)
	if (gi == 0 or gi == 3):
		plt.ylabel("$\sigma_2$", fontsize = 15)
		
for gi in range(6):
	plt.plot(0, 0, col[gi], label = gammalab2[gi])
	



plt.legend(bbox_to_anchor=(-2.4, -0.4, 3.4, 0.2), loc='lower left', ncol=6, mode="expand", borderaxespad=0, fontsize = 15)



plt.subplots_adjust(bottom = 0.2, top = 0.95, right = 0.95, left = 0.08, wspace = 0.2, hspace = 0.14)
plt.savefig("phasesim2.pdf")



#plt.show()
