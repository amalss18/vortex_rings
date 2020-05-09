from mayavi import mlab
import numpy as np

x, y, z = np.mgrid[-10:10:100j, -10:10:100j, -10:10:100j]

x = x.ravel()
y = y.ravel()
z = z.ravel()

theta = np.arctan(y/x)

r = np.sqrt(x**2+y**2)
# r = x

ur, utheta, uz = np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)

v0 = 1
c = 3
a = 1

uz = v0 * c**4 * (r**2 - z**2 - c**2) / (r**2 + z**2 + c**2)**3
ur = -2 * v0 * (c**4) * r * z / (r**2 + z**2 + c**2)**3

ux = ur * np.cos(theta)
uy = ur * np.sin(theta)

idx = ((c - np.sqrt(x**2 + y**2))**2 + z**2 - a**2) < 0

x = x[idx]
y = y[idx]
z = z[idx]

ux = ux[idx]
uy = uy[idx]
uz = uz[idx]

mlab.figure(bgcolor=(0,0,0))
mlab.quiver3d(x,y,z,ux,uy,uz)
mlab.points3d(x,y,z,color = (0.2,0.4,1))