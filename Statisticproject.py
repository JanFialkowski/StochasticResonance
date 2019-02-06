import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sp
import pickle

eps = 0.05
a = 1.001
D = 0.0001
sigma = 0.4
sigma2 = 0.4
N = 500
R = int(N*0.12)
Phi = np.pi/2-0.1
Tmax = 100
dt = 0.01
dz = 25
Multiplex = True

def fhnderiv(u):
	us = u[:,0]
	vs = u[:,1]
	udifs = (us-us[:,np.newaxis])
	vdifs = (vs - vs[:,np.newaxis])
	udifs *= Rangetrix
	vdifs *= Rangetrix
	du = np.zeros_like(u)
	du[:,0] = (us-us**3/3.-vs + buu*udifs.sum(axis=1)+buv*vdifs.sum(axis=1))/eps
	du[:,1] = us + a + bvu*udifs.sum(axis=1)+buu*vdifs.sum(axis=1)
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
for i in xrange(len(T)-1):
	uv[i+1]=Calcnewuv(uv[i],dt)
	if T[i]%1.==0:
		print T[i]

Z = Order(uv[-int(20/dt):], dz)
#"""
plt.imshow(Z[-int(20/dt):,:], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax), cmap = "jet")
plt.title("Z")
plt.colorbar()
plt.show()
#plt.imshow(uv[:,:,0], origin = "lower", aspect = "auto", extent = (0,N,0,Tmax), cmap = "jet")
#plt.title("uglobal")
#plt.colorbar()
#plt.show()
plt.imshow(uv[-int(20/dt):,:,0], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax), cmap = "jet")
plt.title("lastfewu")
plt.colorbar()
plt.show()
plt.imshow(uv[:500,:,0], origin = "lower", aspect = "auto", cmap = "jet")
plt.show()
#"""
#with open("Chimeralotsotime.pickle","wb") as f:
#	pickle.dump(uv,f)
"""
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
pic1 = ax1.imshow(uv[-int(20/dt):,:,0], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax))
pic2 = ax2.imshow(Z[-int(20/dt):,:], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax))
"""
#D = 0.0002
#for i in xrange(len(T)-1):
#	uv[i+1]=Calcnewuv(uv[i],dt)
#	if T[i]%1.==0:
#		print T[i]
#Z = Order(uv, dz)
#with open("bullshitlotsotime.pickle","wb") as f:
#	pickle.dump(uv[-int(20/dt):],f)
"""
pic3 = ax3.imshow(uv[-int(20/dt):,:,0], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax))
pic4 = ax4.imshow(Z[-int(20/dt):,:], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax))
fig.colorbar(pic1,ax = ax3, orientation="horizontal", pad = 0.2)
fig.colorbar(pic2, ax = ax4, orientation="horizontal", pad=0.2)
plt.show()
"""
