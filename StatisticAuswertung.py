import matplotlib.pyplot as plt
import pickle
import numpy as np
import scipy.linalg as sp
dt = 0.01
N = 500
Tmax = 500
dz = 25
Multiplex = False

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

with open("Chimera.pickle","rb") as f:
	uv = pickle.load(f)[0]
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
pic1 = ax1.imshow(uv[-int(20/dt):,:,0], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax), cmap = "jet")
uv = Order(uv,dz)
pic2 = ax2.imshow(uv[-int(20/dt):,:], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax), cmap = "jet")
with open("bullshitlotsotime.pickle","rb") as f:
	uv = pickle.load(f)
pic3 = ax3.imshow(uv[-int(200/dt):,:,0], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax), cmap = "jet")
uv = Order(uv,dz)
pic4 = ax4.imshow(uv[-int(200/dt):,:], origin = "lower", aspect = "auto", extent = (0,N,Tmax-20,Tmax), cmap = "jet")
fig.colorbar(pic1,ax = ax3, orientation="horizontal", pad = 0.2)
fig.colorbar(pic2, ax = ax4, orientation="horizontal", pad=0.2)
plt.show()
