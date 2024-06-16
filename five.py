# fiveplots for higher order
import numpy as np
import matplotlib.pyplot as plt
import math

col = ['m', 'r', 'y', 'g', 'c', 'b']
col_ = ['m--', 'r--', 'y--', 'g--', 'c--', 'b--']
colo = ['m.', 'r.', 'y.', 'g.', 'c.', 'b.']
lett = ['$(a)$', '$(b)$', '$(c)$', '$(d)$', '$(e)$', '$(f)$']

gammalab2 = ['$\gamma_2 = -1$', '$\gamma_2 = -0.8$', '$\gamma_2 = -0.5$', '$\gamma_2 = 0$', '$\gamma_2 = 0.5$', '$\gamma_2 = 1$']
gammalab3 = ['$\gamma_3 = -0.5$', '$\gamma_3 = -0.4$', '$\gamma_3 = -0.3$', '$\gamma_3 = 0$', '$\gamma_3 = 0.5$', '$\gamma_3 = 1$']
gammalab4 = ['$\gamma_4 = -1/3$', '$\gamma_4 = -0.3$', '$\gamma_4 = -0.2$', '$\gamma_4 = 0$', '$\gamma_4 = 0.5$', '$\gamma_4 = 1$']
gammalab = [gammalab2, gammalab3, gammalab4]

muval = np.linspace(-4, 2, 25)

x = np.linspace(-1.0, 1.0, 21)
sigmaval = np.power(10, x)

plt.subplots(3, 3, figsize = (12,10))

for pi in range(3):
	
	pi += 2

	points = [] # gi, 012, mi
	for gi in range(6):
		p = np.loadtxt("points{0}_{1}.txt".format(pi, gi), delimiter = ",")
		p = p.reshape(25, 3)
		p = np.transpose(p)
		points.append(p)
	
	if (pi != 3):	
		meas = np.loadtxt("allfive{0}.txt".format(pi), delimiter=",")
		meas = meas.reshape(25, 6, 21, 3)
		meas = np.transpose(meas, (0, 1, 3, 2)) # mu, gam, meas, sig
	
	mi = 20

	measfp = []
	for gi in range(6):
		g = np.loadtxt("measfp{0}_{1}_{2}.txt".format(pi, mi, gi), delimiter = ",")
		g = g.reshape(len(g)//6, 6)
		g = np.transpose(g)
		measfp.append(g)

	plt.subplot(3, 3, pi - 1) # phi
	for gi in range(6):
		plt.semilogx(measfp[gi][0], measfp[gi][2], col[gi])
		if (pi != 3):
			plt.semilogx(sigmaval, meas[mi][gi][0], colo[gi], label = gammalab[pi-2][gi])           
		plt.axvline(x = points[gi][0][mi], linestyle = '--', color = col[gi])
		plt.axvline(x = points[gi][1][mi], linestyle = ':', color = col[gi])
		plt.axvline(x = points[gi][2][mi], linestyle = ':', color = col[gi])
	plt.xlim([sigmaval[0], sigmaval[-1]])
	if (pi == 2):
		plt.ylabel("$\phi$", fontsize = 15)
	if (pi == 3):
		plt.ylim([0.477, 1.025])
	plt.title("$\mu_{0} = {1}$".format(pi, muval[mi]), fontsize = 15)

	plt.subplot(3, 3, pi + 2) # M
	for gi in range(6):
		plt.semilogx(measfp[gi][0], measfp[gi][3], col[gi])
		if (pi != 3):
			plt.semilogx(sigmaval, meas[mi][gi][1], colo[gi], label = gammalab[pi-2][gi])
		plt.axvline(x = points[gi][0][mi], linestyle = '--', color = col[gi])
		plt.axvline(x = points[gi][1][mi], linestyle = ':', color = col[gi])
		plt.axvline(x = points[gi][2][mi], linestyle = ':', color = col[gi])
	plt.xlim([sigmaval[0], sigmaval[-1]])
	if (pi == 2):
		plt.ylim([0, 1])
		plt.ylabel("$M^*$", fontsize = 15)
	if (pi == 3):
		plt.ylim([0.25, 0.5])
		plt.yticks(np.linspace(0.25, 0.5, 6))
	if (pi == 4):
		plt.ylim([0.35, 0.55])
		plt.yticks(np.linspace(0.35, 0.55, 5))
	
	plt.subplot(3, 3, pi + 5) # q
	for gi in range(6):
		plt.semilogx(measfp[gi][0], np.divide(measfp[gi][4]*measfp[gi][3], abs(measfp[gi][3])), col[gi], label = gammalab[pi-2][gi])
		#plt.semilogx(measfp[gi][0], measfp[gi][4], col[gi], label = gammalab[pi-2][gi])
		if (pi != 3):
			plt.semilogx(sigmaval, meas[mi][gi][2], colo[gi])
		plt.axvline(x = points[gi][0][mi], linestyle = '--', color = col[gi])
		plt.axvline(x = points[gi][1][mi], linestyle = ':', color = col[gi])
		plt.axvline(x = points[gi][2][mi], linestyle = ':', color = col[gi])
	plt.xlim([sigmaval[0], sigmaval[-1]])
	plt.xlabel("$\sigma_{0}$".format(pi), fontsize = 15)
	if (pi == 2):
		plt.ylabel("$q$", fontsize = 15)
		plt.ylim([0, 1])
	if (pi == 3):
		plt.ylim([0.1, 0.6])
		plt.yticks(np.linspace(0.1, 0.6, 6))
	if (pi == 4):
		plt.ylim([0.2, 0.6])
		plt.yticks(np.linspace(0.2, 0.6, 5))
	plt.legend(bbox_to_anchor=(0, -0.7, 1, 0.2), loc='lower left', ncol=2, mode="expand", borderaxespad=0, fontsize = 15)
	
plt.subplots_adjust(left = 0.06, right = 0.98, top = 0.96, wspace = 0.15, hspace = 0.2, bottom = 0.17)


plt.savefig("fivep_{0}.pdf".format(mi))
plt.show()	
	
