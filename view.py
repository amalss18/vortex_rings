import numpy as np
from mayavi import mlab
from poisson import solve

n = 100
sigma = 0.2
yy = np.linspace(0, 1, n)
xx = np.linspace(0, 1, n)
y, x = np.meshgrid(yy, xx)

# the right hand side term in Poisson eq.
phi = np.exp(-(y**2+x**2) / (2*sigma**2))

# the initial guess
u0 = np.zeros_like(phi)

# solve the Poisson equation with fixed boundaries
u = solve(u0, phi, fix="bounds")

print(u)