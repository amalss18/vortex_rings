from pysph.base.utils import get_particle_array
from pysph.solver.application import Application
from pysph.sph.scheme import WCSPHScheme
import numpy as np
from numpy import ones_like, mgrid, sqrt

class VortexRing(Application):
    def initialize(self):
        self.co = 10
        self.ro = 1.0
        self.hdx = 1.3
        self.dx = 0.2
        self.alpha = 0.1

    def create_scheme(self):
        s = WCSPHScheme(
            ['fluid'], [], dim=2, rho0=self.ro, c0=self.co,
            h0=self.dx*self.hdx, hdx=self.hdx, gamma=7.0, alpha=0.1, beta=0.0
        )
        dt = 5e-4
        tf = 1.5
        s.configure_solver(dt=dt, tf=tf)
        return s

    def create_particles(self):
        """Create the circular patch of fluid."""

        dx = self.dx
        ro = self.ro
        hdx = self.hdx

        x, y, z = np.mgrid[-10:10:100j, -10:10:100j, -25:10:250j]

        x = x.ravel()
        y = y.ravel()
        z = z.ravel()

        theta = np.arctan(y/x)

        r = np.sqrt(x**2+y**2)

        ur, utheta, uz = np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)

        v0 = 1
        self.co = 10*v0
        c = 3
        a = 1

        uz = v0 * c**4 * (r**2 - z**2 - c**2) / (r**2 + z**2 + c**2)**3
        ur = -2 * v0 * (c**4) * r * z / (r**2 + z**2 + c**2)**3

        ux = ur * np.cos(theta)
        uy = ur * np.sin(theta)

        name = 'fluid'

        m = ones_like(x)*dx*dx*ro
        h = ones_like(x)*hdx*dx
        rho = ones_like(x) * ro

        # remove particles outside the circle
        indices = []
        for i in range(len(x)):
            if not (((c - r[i])**2 + z[i]**2 - a**2) < 0):
                indices.append(i)

        pa = get_particle_array(x=x, y=y, z=z, m=m, rho=rho, h=h, u=ux, v=uy, w=uz,
                                name=name)
        pa.u[indices] = 0
        pa.v[indices] = 0
        pa.w[indices] = 0

        print("Vortex ring :: %d particles"
              % (pa.get_number_of_particles()))

        self.scheme.setup_properties([pa])
        return [pa]


if __name__ == '__main__':
    app = VortexRing()
    app.run()
