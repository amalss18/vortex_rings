import numpy as np
from mayavi import mlab

x, y, z = np.mgrid[-3:3:30j, -3:3:30j, -3:3:30j]

a = 1.5

# x = np.zeros_like(x)
# x = np.ones_like(x)*-1
# x = x.ravel()
# y = y.ravel()
# z = z.ravel()
v0 = 1.5

r = np.sqrt(x**2 + y**2)
Uz = v0*(a**4)*(r**2 - z**2 -a**2) / (r**2 + z**2 + a**2)**3
Ur = -2*v0*(a**4)*r*z / (r**2 + z**2 + a**2)**3

Uz = Uz + 0.09*v0
theta = np.arctan(y/x)
idx = x<0
theta[idx] = theta[idx] - np.pi
Ux = Ur * np.cos(theta)
Uy = Ur * np.sin(theta)

src = mlab.pipeline.vector_field(x, y, z, Ux, Uy, Uz)
mlab.pipeline.vector_cut_plane(src)
mlab.pipeline.vectors(src)
# mlab.show_pipeline()