import numpy as np
import matplotlib.pyplot as plt
import math

col = ['m', 'r', 'y', 'g', 'c', 'b']
col_ = ['m--', 'r--', 'y--', 'g--', 'c--', 'b--']
colo = ['m.', 'r.', 'y.', 'g.', 'c.', 'b.']
lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']

gammalab3 = ['$\gamma_3 = -0.5$', '$\gamma_3 = -0.4$', '$\gamma_3 = -0.3$', '$\gamma_3 = 0$', '$\gamma_3 = 0.5$', '$\gamma_3 = 1$']

muval = np.linspace(-4, 2, 25)
x = np.linspace(-1.0, 1.0, 21)
sigmaval = np.power(10, x)

points = [] # gi, 012, mi
for gi in range(6):
	p = np.loadtxt("points3_{0}.txt".format(gi), delimiter = ",")
	p = p.reshape(25, 3)
	p = np.transpose(p)
	points.append(p)
	

	

#meas = np.loadtxt("allfive3.txt", delimiter=",")
#meas = meas.reshape(25, 6, 21, 3)
#meas = np.transpose(meas, (0, 1, 3, 2)) # mu, gam, meas, sig


allmeasfp = []
for mi in range(25):
	measfp = []
	for gi in range(6):
		g = np.loadtxt("measfp3_{0}_{1}.txt".format(mi, gi), delimiter = ",")
		g = g.reshape(len(g)//6, 6)
		g = np.transpose(g)
		measfp.append(g)
	allmeasfp.append(measfp)
	
	
# all mu on same plot
##################################################
	
plt.subplots(2, 3, figsize = (12, 8))
for gi in range(6):
	
	plt.subplot(2, 3, gi + 1)
	
	for mi in range(25):
		plt.semilogx(allmeasfp[mi][gi][0], allmeasfp[mi][gi][3], col[gi])
		#plt.semilogx(sigmaval, meas[mi][gi][1], 'k.')#, label = gammalab[gi])
	plt.xlim([sigmaval[0], sigmaval[-1]])
	plt.ylim([0, 1])
	if (gi == 3 or gi == 4 or gi == 5):
		plt.xlabel("$\sigma_3$", fontsize = 15)
	if (gi == 0 or gi == 3):
		plt.ylabel("$M^*$", fontsize = 15)
		
for gi in range(6):
	plt.plot(0, 0, col[gi], label = gammalab3[gi])
	
plt.legend(bbox_to_anchor=(-2.4, -0.4, 3.4, 0.2), loc='lower left', ncol=6, mode="expand", borderaxespad=0, fontsize = 15)

plt.subplots_adjust(bottom = 0.2, top = 0.95, right = 0.95, left = 0.08, wspace = 0.2, hspace = 0.14)
plt.savefig("allMplots3.pdf")


# examples of M for the paper
###############################################################
mival = [8, 14, 17, 18]
Mlim = [1, 1, 2, 2]

plt.subplots(1, 4, figsize=(12,4))

for mi in range (4):

	plt.subplot(1, 4, 1 + mi) # M
	for gi in range(6):
		plt.semilogx(allmeasfp[mival[mi]][gi][0], allmeasfp[mival[mi]][gi][3], col[gi])
		#plt.semilogx(sigmaval, meas[mival[mi]][gi][1], colo[gi])
		plt.axvline(x = points[gi][0][mival[mi]], linestyle = '--', color = col[gi])
		plt.axvline(x = points[gi][1][mival[mi]], linestyle = ':', color = col[gi])
		plt.axvline(x = points[gi][2][mival[mi]], linestyle = ':', color = col[gi])
	plt.xlim([sigmaval[0], sigmaval[-1]])
	plt.title("$\mu_3 = {0}$".format(muval[mival[mi]]), fontsize = 15)
	plt.ylim([0, Mlim[mi]])
	plt.xlabel("$\sigma_3$", fontsize = 15)
	if (mi == 0):
		plt.ylabel("$M^*$", fontsize = 15)
		
	for gi in range(6):
		plt.plot(0, 0, col[gi], label = gammalab3[gi])
	
	plt.subplots_adjust(bottom = 0.25, top = 0.92, hspace = 0.36, left = 0.05, right = 0.97)
	
plt.legend(bbox_to_anchor=(-3.6, -0.35, 4.6, 0.2), loc='lower left', ncol=6, mode="expand", borderaxespad=0, fontsize = 15)
	
plt.savefig("Mplots3.pdf")

