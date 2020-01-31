from mayavi import mlab
import numpy as np

x, z = np.mgrid[-10:10:100j, -10:10:100j]

x = x.ravel()
z = z.ravel()

y = np.zeros_like(x)

theta = np.arctan(y/x)

r = x
# r = x

ur, utheta, uz = np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)

v0 = 1
c = 1.5
a = 1

uz = v0 * c**4 * (r**2 - z**2 - c**2) / (r**2 + z**2 + c**2)**3
ur = -2 * v0 * (c**4) * r * z / (r**2 + z**2 + c**2)**3

ux = ur
uy = np.zeros_like(ux)

idx = ((c - x)**2 + z**2 - a**2) < 0
idxo = ((c + x)**2 + z**2 - a**2) < 0

x1 = x[idx]
y1 = y[idx]
z1 = z[idx]

ux1 = ux[idx]
uy1 = uy[idx]
uz1 = uz[idx]

x2= x[idxo]
y2= y[idxo]
z2= z[idxo]

ux2 = ux[idxo]
uy2 = uy[idxo]
uz2 = uz[idxo]

x = np.concatenate((x1,x2))
y = np.concatenate((y1,y2))
z = np.concatenate((z1,z2))

ux = np.concatenate((ux1,ux2))
uy = np.concatenate((uy1,uy2))
uz = np.concatenate((uz1,uz2))


mlab.figure(bgcolor=(0,0,0))
mlab.quiver3d(x,y,z,ux,uy,uz)
points = mlab.points3d(x,y,z,color = (0.2,0.4,1))
points.glyph.glyph.scale_factor = 0.05
mlab.gcf().scene.y_plus_view()
