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

with open("UCoupling0001.pickle","rb") as f:
	uv = pickle.load(f)

plt.imshow(uv[-int(20/dt):,:,0], origin = "lower", aspect = "auto", extent = (0,2*N,Tmax-20,Tmax), cmap = "jet")
plt.colorbar()
plt.show()
