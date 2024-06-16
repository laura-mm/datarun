from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math
import mpld3
import plotly
import plotly.tools as tls
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import numpy
numpy.set_printoptions(threshold=sys.maxsize)


# codes for plots plotly surfaces
# F gam2 = 0, gam3 = 0
# G gam2 = -1, gam3 = -0.5
# H gam2 = 1, gam3 = 1
# I gam2 = 0.5, gam3 = 0.5
# J gam2 = 0.4, gam3 = 0.4
# K gam2 = -1, gam3 = 1
# L gam2 = 1, gam3 = -0.5
# M gam2 = -0.8, gam3 = -0.4

crit3 = np.loadtxt("Pcrit233G.txt", delimiter = ",")
crit3 = crit3.reshape(101, 1501, 5) 
#crit3 = np.transpose(crit3, (2, 0, 1)) # (mu2, mu3, sig2, sig3, maxsig), mu2, mu3
#crit3m = np.transpose(crit3, (2, 0, 1)) # (mu2, mu3, sig2, sig3, maxsig), mu2, mu3

div2 = np.loadtxt("Pdiv232G.txt", delimiter = ",")
div2 = div2.reshape(101, 25000, 3) # mu2, z1, (mu2, mu3, sigd)
div2e = np.transpose(div2, (2, 0, 1)) # (mu2, mu3, sigd), mu2, z1
#print(div2[1])
#divline2 = np.transpose(div2, (1, 0, 2))

div3 = np.loadtxt("Pdiv233G.txt", delimiter = ",")
div3 = div3.reshape(101, 25001, 3) # mu2, z1, (mu2, mu3, sig3)
#div3 = np.transpose(div3, (2, 0, 1)) # (mu2, mu3, sig3), mu2, z1

crit2b = np.loadtxt("Pcrit232G.txt", delimiter = ",")
crit2b = crit2b.reshape(101, 3) # mu2, (mu2, mu3, sigc)
crit2b = np.transpose(crit2b) # (mu2, mu3, sigc), mu2

crit2ai = [-4 for i in range(101)]
crit2a = [crit2b[0], crit2ai, crit2b[2]]

crit2 = [crit2a, crit2b]
crit2 = np.transpose(crit2, (1, 2, 0))

extra3 = np.loadtxt("Pextra3G.txt", delimiter = ",")
extra3 = extra3.reshape(101, 5) # mu2, (mu2, lcstop, lmstop, jdstart, jdstop)
extra3 = np.transpose(extra3)

#print(div2e[2][0])
div2e = np.where(np.isinf(div2e), -np.Inf, div2e)
extra2 = [np.nanargmax(div2e[2][i]) for i in range(101)]
#print(extra2[0])
#print(div2e[2][0][extra2[0]])

cmin = min(extra3[1]) 
print(cmin)
mmin = min(extra3[2])
print(mmin)
jdmin = min(filter(lambda i: i > 0, (extra3[4] - extra3[3])))
print(jdmin)
d2min = min(extra2)
print(d2min)
d2min = 1000
#jdmin = 1000

div2d = [[div2[i][int(j*extra2[i]/d2min)] for j in range(int(d2min) + 1)] for i in range(101)]
div2d = np.transpose(div2d, (2, 0, 1)) # (mu2, mu3, sigd), mu2, z1

div3d = []
for i in range(101):
	if extra3[3][i] == -1:
		temp = [[div3[i][0][0], 0.0, (j*10.0/(int(jdmin)))] for j in range(int(jdmin) + 1)] # for antisymmetric
		#temp = [[div3[i][0][0], 0.0, math.inf] for j in range(int(jdmin) + 1)] # for gam3 = -0.4
		
	else:
		temp = [div3[i][int(extra3[3][i] + round(j*(extra3[4][i] - extra3[3][i])/jdmin))] for j in range(int(jdmin) + 1)]
	div3d.append(temp)
div3d = np.transpose(div3d, (2, 0, 1))

print(div3d)

if (cmin != -1 and cmin != 0):
	crit3c = [[[crit3[i][int(round(l*extra3[1][i]/cmin))][0], crit3[i][int(round(l*extra3[1][i]/cmin))][1], crit3[i][int(round(l*extra3[1][i]/cmin))][3]] for l in range(int(cmin) + 1)] for i in range(101)] # mu2, mu3, (mu2, mu3, sig3)
	crit3c = np.transpose(crit3c, (2, 0, 1))
else:
	crit3c = np.transpose(crit3, (2, 0, 1))
	

if (mmin != -1 and mmin != 0):
	crit3m = [[[crit3[i][int(round(l*extra3[2][i]/mmin))][0], crit3[i][int(round(l*extra3[2][i]/mmin))][1], crit3[i][int(round(l*extra3[2][i]/mmin))][4]] for l in range(int(mmin) + 1)] for i in range(101)] # mu2, mu3, (mu2, mu3, sig3)
	crit3m = np.transpose(crit3m, (2, 0, 1))
else:
	crit3m = np.transpose(crit3, (2, 0, 1))

surfacecolorc = np.ones(crit3c[0].shape)
surfacecolorm = np.ones(crit3m[0].shape)
surfacecolord = np.ones(div2d[0].shape)
surfacecolorc2 = np.ones(crit2[0].shape)
surfacecolord3 = np.ones(div3d[0].shape)


fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'surface'}, {'type': 'surface'}]], subplot_titles=("σ<sub>3</sub> = 0", "σ<sub>2</sub> = 0"))


fig.add_trace(go.Surface(z=div2d[2], x=div2d[0], y=div2d[1], opacity=0.5, surfacecolor = surfacecolord, colorscale = 'ice', showscale=False), row=1, col=1)
fig.add_trace(go.Surface(z=crit2[2], x=crit2[0], y=crit2[1], opacity=0.5, surfacecolor = surfacecolorc2, colorscale = 'reds', showscale=False), row=1, col=1)
fig.add_trace(go.Scatter3d(x = [-2, -2], y = [4, 4], z = [0, 10], mode = 'lines', line = dict(width = 5, color = 'black', showscale=False)), row=1, col=1)
fig.add_trace(go.Scatter3d(x = [-2, -2], y = [-4, -4], z = [0, 10], mode = 'lines', line = dict(width = 5, color = 'black', showscale=False)), row=1, col=1)
#fig = go.Figure(data=[go.Scatter3d(x=[0,0], y=[0,0], z=[0,10], mode='markers')], row=1, col=1)
#fig.update_scenes(zaxis_type="log")
#fig.update_layout(scene = dict(zaxis = dict(range=[-2,2])))
fig.update_layout(scene = dict(xaxis = dict(range=[-4,6])))
fig.update_layout(scene = dict(yaxis = dict(range=[-4,6])))
fig.update_layout(scene = dict(zaxis = dict(range=[0,10])))
fig.update_layout(scene = dict(xaxis_title= "μ<sub>2</sub>", xaxis_title_font_size=20, yaxis_title="μ<sub>3</sub>", yaxis_title_font_size=20, zaxis_title="σ<sub>2</sub>", zaxis_title_font_size=20))
fig.update_layout(scene = dict(aspectmode='cube'))





if (cmin != -1 and cmin != 0):
	fig.add_trace(go.Surface(z=crit3c[2], x=crit3c[0], y=crit3c[1], opacity=0.5, surfacecolor = surfacecolorc, colorscale = 'reds', showscale=False), row=1, col=2)
else:
	fig.add_trace(go.Surface(z=crit3c[3], x=crit3c[0], y=crit3c[1], opacity=0.5, surfacecolor = surfacecolorc, colorscale = 'reds', showscale=False), row=1, col=2)
#if (mmin != -1 and mmin != 0):
#	fig.add_trace(go.Surface(z=crit3m[2], x=crit3m[0], y=crit3m[1], opacity=0.5, surfacecolor = surfacecolorm, colorscale = 'greens', showscale=False), row=1, col=2)
#else:
#	fig.add_trace(go.Surface(z=crit3m[4], x=crit3m[0], y=crit3m[1], opacity=0.5, surfacecolor = surfacecolorm, colorscale = 'greens', showscale=False), row=1, col=2)
fig.add_trace(go.Surface(z=div3d[2], x=div3d[0], y=div3d[1], opacity=0.5, surfacecolor = surfacecolord3, colorscale = 'ice', showscale=False), row=1, col=2)
fig.add_trace(go.Scatter3d(x = [-2, -2], y = [4, 4], z = [0, 10], mode = 'lines', line = dict(width = 5, color = 'black', showscale=False)), row=1, col=2)
fig.add_trace(go.Scatter3d(x = [-2, -2], y = [-4, -4], z = [0, 10], mode = 'lines', line = dict(width = 5, color = 'black', showscale=False)), row=1, col=2)
#fig.update_scenes(zaxis_type="log")
#fig.update_layout(scene = dict(zaxis = dict(range=[-2,2])))
fig.update_layout(scene2 = dict(xaxis = dict(range=[-4,6])))
fig.update_layout(scene2 = dict(yaxis = dict(range=[-4,6])))
fig.update_layout(scene2 = dict(zaxis = dict(range=[0,10])))
fig.update_layout(scene2 = dict(xaxis_title="μ<sub>2</sub>", xaxis_title_font_size=20, yaxis_title="μ<sub>3</sub>", yaxis_title_font_size=20, zaxis_title="σ<sub>3</sub>", zaxis_title_font_size=20))
fig.update_layout(scene2 = dict(aspectmode='cube'))


fig.update_layout(title = "γ<sub>2</sub> = -1, γ<sub>3</sub> = -0.5", title_font_size=20)
fig.update_layout(showlegend=False)
fig.update_annotations(font_size=20)
fig.update_scenes(xaxis_tickfont_size=15)
fig.update_scenes(yaxis_tickfont_size=15)
fig.update_scenes(zaxis_tickfont_size=15)

fig.show()
fig.write_html("./plotlylG.html")
fig.write_image("./plotlylG.png")



