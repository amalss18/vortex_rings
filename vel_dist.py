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
c = 5
a = 2
r0 = c
K = 5
# uz = v0 * c**4 * (r**2 - z**2 - c**2) / (r**2 + z**2 + c**2)**3
# ur = -2 * v0 * (c**4) * r * z / (r**2 + z**2 + c**2)**3
pi = np.pi
ur = (K * z /(2* pi * ((r-r0)**2 + z**2))) - (K*(r-r0)*z/(4 * pi * r0 * ((r-r0)**2 + z**2)))
uz = (-K *(r-r0) / (2 * pi * (r-r0)**2 + z**2)) + (K/(4*pi*r0)) * (np.log(8*r0/(np.sqrt((r-r0)**2 + z**2))) - (z**2/((r-r0)**2 + z**2)))
# uz = 1.5 * v0 * (1-((2*r**2+z**2)/a**2))
# ur = 1.5 * v0 * z*r/a**2


ux = ur * np.cos(theta)
uy = ur * np.sin(theta)

idx = ((c - np.sqrt(x**2 + y**2))**2 + z**2 - a**2) < 0
# idx = r**2 + z**2 < a**2

x = x[idx]
y = y[idx]
z = z[idx]

ux = ux[idx]
uy = uy[idx]
uz = uz[idx]

mlab.figure(bgcolor=(0,0,0))
mlab.quiver3d(x,y,z,ux,uy,uz)
# points = mlab.points3d(x,y,z,color = (0.2,0.4,1))
# points.glyph.glyph.scale_factor = 0.05
mlab.gcf().scene.y_plus_view()
