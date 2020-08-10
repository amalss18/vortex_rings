from pysph.base.utils import get_particle_array
from pysph.solver.application import Application
from pysph.sph.scheme import WCSPHScheme
import numpy as np
from pysph.base.nnps import DomainManager
from pysph.base.kernels import CubicSpline
from pysph.sph.integrator_step import WCSPHStep
from pysph.sph.integrator import PECIntegrator
from pysph.solver.solver import Solver
from pysph.tools.geometry import rotate
from pysph.base.nnps import DomainManager

# Geometry and SPH parameters
R = 3.0 # Radius of ring
a = 1.0  # Radius of core
nu = 1e-3
print("Nu: ", nu)
Re = 7500
Gamma = Re * nu
print("Gamma = ", Gamma)
n = 30

length = R + 2.5
# dx = 2*length/n
dx = 0.35
print("dx: ", dx)

rho = 1000.0
hdx = 1.32
h = hdx * dx
m = rho * dx * dx * dx

vel_scale = 1.0

# vmax = 0.15*vel_scale*Gamma
vmax = 20.651
c0 = 10 * vmax

dt_cfl = 0.25 * h / (c0 + vmax)
dt_viscous = 0.125 * h * h / nu
dt_wcsph = 0.125 * h / c0

print(dt_cfl, dt_viscous, dt_wcsph)

dt = 0.5 * min(dt_cfl, dt_viscous, dt_wcsph)

class VortexRing(Application):
    def create_domain(self):
        # domain for periodicity
        domain = DomainManager(
            xmin=-length-1, xmax=length+1, ymin=-length-1, ymax=length+1,
            zmin=-length+0.75, zmax=length+3.5, periodic_in_x=True,
            periodic_in_y=True, periodic_in_z=True
        )
        return domain

    def create_particles(self):
        x, y, z = np.mgrid[-length-1:length+1:dx, -length-1:length+1:dx, -length+0.75:length+3.5:dx]
        x = x.ravel()
        y = y.ravel()
        z = z.ravel()

        print(len(x))
        idx = (np.sqrt(x**2 + y**2) - R)**2 + z**2 < a**2

        data = np.load("vel_035.npz")
        velocity = data['arr_0']
        print(velocity.shape)
        u, v, w = velocity.T[0], velocity.T[1], velocity.T[2]
        u, v, w = np.array([u, v, w]) * vel_scale

        fluid = get_particle_array(x=x, y=y, z=z, m=m, rho=rho, h=h, u=u, v=v, w=w,
                                   name="fluid")


        self.scheme.setup_properties([fluid])
        return [fluid]

    def create_scheme(self):
        s = WCSPHScheme(
            ['fluid'], [], dim=3, rho0=rho, c0=c0,
            h0=h, hdx=hdx, gamma=7.0, alpha=0.02, beta=0.0, nu=nu
        )
        return s

    def create_solver(self):
        kernel = CubicSpline(dim=3)
        integrator = PECIntegrator(fluid=WCSPHStep())

        # dt = 0.01*0.125*h/c0
        # dt = 0.125*h/c0
        print("Time step: ", dt)
        tf = 1
        solver = Solver(kernel=kernel, dim=3, integrator=integrator,
                        tf=tf, dt=dt, adaptive_timestep=False,
                        fixed_h=False)
        return solver

if __name__=="__main__":
    app = VortexRing()
    app.run()