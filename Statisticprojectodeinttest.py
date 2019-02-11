import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sp
from scipy.integrate import odeint
import pickle

eps = 0.05
a = 1.001
D = 0.0001
sigma = 0.4
sigma2 = 0.01
N = 500
R = int(N*0.12)
Phi = np.pi/2-0.1
Tmax = 200
dt = 0.001
dz = 25
Multiplex = True

def fhnderiv(u):
	us = u[:,0]
	vs = u[:,1]
	udifs = (us-us[:,np.newaxis])
	vdifs = (vs - vs[:,np.newaxis])
	udifs[:N,:N] *= Rangetrix[:N,:N]
	vdifs[:N,:N] *= Rangetrix[:N,:N]
	udifs[N:,N:] *= Rangetrix[N:,N:]
	vdifs[N:,N:] *= Rangetrix[N:,N:]
	du = np.zeros_like(u)
	du[:N,0] = (us[:N]-us[:N]**3/3.-vs[:N] + buu*udifs[:N,:N].sum(axis=1)+buv*vdifs[:N,:N].sum(axis=1)+np.diagonal(udifs[N:,:N]*sigma2))/eps
	du[N:,0] = (us[N:]-us[N:]**3/3.-vs[N:] + buu*udifs[N:,N:].sum(axis=1)+buv*vdifs[N:,N:].sum(axis=1)+np.diagonal(udifs[:N,N:]*sigma2))/eps
	du[N:,1] = us[N:] + a + bvu*udifs[N:,N:].sum(axis=1)+buu*vdifs[N:,N:].sum(axis=1)
	du[:N,1] = us[:N] + a + bvu*udifs[:N,:N].sum(axis=1)+buu*vdifs[:N,:N].sum(axis=1)
	return du

def Calcnewuv(uv,dt):
	uv += fhnderiv(uv)*dt
	uv[:N,1]+=np.sqrt(2*D*dt)*np.random.standard_normal(N)
	return uv

def Order(u,dz):
	Phis = np.arctan(u[:,:,1]/u[:,:,0])
	dummy = np.zeros(N)
	index  = np.arange(N)
	dummy[index<=dz] = 1
	dummy[len(index)-index<=dz] = 1
	Rangetrix = sp.circulant(dummy)
	del dummy, index
	if Multiplex == True:
		Rangetrix = np.bmat([[Rangetrix,np.zeros((N,N))],[np.zeros((N,N)),Rangetrix]]).A
	Z = np.exp(1j*Phis)
	Z = np.array([np.dot(Rangetrix, Vector) for Vector in Z])
	Z /= (2*dz+1)
	Z = np.abs(Z)
	return Z

u = np.array([[2*np.cos(T),2*np.sin(T)] for T in np.random.rand(N)*2*np.pi])
if Multiplex == True:
	u = np.array([[2*np.cos(T),2*np.sin(T)] for T in np.random.rand(2*N)*2*np.pi])
T = np.arange(0,Tmax+dt,dt)
uv = np.zeros((len(T),u.shape[0],u.shape[1]))
uv[0] = u
buu = np.cos(Phi)
buv = np.sin(Phi)
bvu = np.sin(-Phi)
dummy = np.zeros(N)
index  = np.arange(N)
dummy[index<=R] = 1
dummy[len(index)-index<=R] = 1
Rangetrix = sp.circulant(dummy)
del dummy, index
if Multiplex == True:
	Rangetrix = np.bmat([[sigma/(2*R)*Rangetrix,sigma2*np.eye(N)],[sigma2*np.eye(N),sigma/(2*R)*Rangetrix]]).A
do=True
while do:
	for i in xrange(len(T)-1):
		uv[i+1]=Calcnewuv(uv[i],dt)
		if T[i] == 50:
			do = np.argwhere(uv[int(30/dt):,N+N/2,0]>1.5).any()
			print np.argwhere(uv[int(30/dt):,N+N/2,0]>1.5)
			print np.argwhere(uv[int(30/dt):,N+N/2,0]>1.5).any()
		if T[i]%1.==0:
			print T[i]
	do = False
"""
plt.imshow(uv[-int(20/dt):,:,0], origin = "lower", aspect = "auto", extent = (0,N*(1+Multiplex),Tmax-20,Tmax), cmap = "jet")
plt.title("lastfewu")
plt.colorbar()
plt.show()
plt.imshow(uv[:,:,0], origin = "lower", aspect = "auto", extent = (0,N*(1+Multiplex),0,Tmax), cmap = "jet")
plt.title("uglobal")
plt.colorbar()
plt.show()
Z = Order(uv[-int(20/dt):], dz)
plt.imshow(Z[-int(20/dt):,:], origin = "lower", aspect = "auto", extent = (0,N*(1+Multiplex),Tmax-20,Tmax), cmap = "jet")
plt.title("Z")
plt.colorbar()
plt.show()
"""
with open("UCoupling001.pickle","wb") as f:
	pickle.dump(uv[::10],f)
