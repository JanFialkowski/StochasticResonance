import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sp
eps = 0.05
a = 1.001
D = 0.02
sigma = 0.4
R = 100
N = 500
Phi = np.pi/2-0.1
Tmax = 100
dt = 0.01
def fhnderiv(u):
	us = u[:,0]
	vs = u[:,1]
	udifs = us-us[:,np.newaxis]
	vdifs = vs - vs[:,np.newaxis]
	udifs *= sigma/(2*R)*Rangetrix
	vdifs *= sigma/(2*R)*Rangetrix
	du = np.zeros_like(u)
	du[:,0] = us-us**3/3.-vs + buu*udifs.sum(axis=1)+buv*vdifs.sum(axis=1)
	du[:,1] = us + a + bvu*udifs.sum(axis=1)+buu*vdifs.sum(axis=1)
	return du

u = np.array([[2*np.cos(T),2*np.sin(T)] for T in np.random.rand(N)*2*np.pi])
T = np.arange(0,Tmax,dt)
print T
uv = np.zeros((len(T),u.shape[0],u.shape[1]))
uv[0] = u
buu = np.cos(Phi)
buv = np.sin(Phi)
bvu = np.sin(-Phi)
dummy = np.zeros(len(u[:]))
index  = np.arange(len(u[:]))
dummy[index<=R] = 1
dummy[len(index)-index<=R] = 1
Rangetrix = sp.circulant(dummy).T
del dummy, index
for i in xrange(len(T)-1):
	deriv = fhnderiv(uv[i])
	rand = np.random.standard_normal(N)
	newu = np.zeros_like(uv[i])
	newu[:,0] = uv[i,:,0]+dt*deriv[:,0]
	newu[:,1] = uv[i,:,1] +dt*deriv[:,1] + np.sqrt(dt*2*D)*rand
	uv[i+1]=newu
	print T[i]
showing = uv[-500:,:,0]
plt.imshow(showing, origin = "lower")
plt.show()
