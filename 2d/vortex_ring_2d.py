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
        dt = 1e-3
        tf = 1.5
        s.configure_solver(dt=dt, tf=tf)
        return s

    def create_particles(self):
        """Create the circular patch of fluid."""

        dx = self.dx
        ro = self.ro
        hdx = self.hdx
        x, z = np.mgrid[-10:10:100j, -40:10:500j]
        x = x.ravel()
        z = z.ravel()

        y = np.zeros_like(x)

        theta = 0

        r = x

        ur, utheta, uz = np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)

        v0 = 5
        self.co = v0
        c = 1.5
        a = 1

        uz = v0 * c**4 * (r**2 - z**2 - c**2) / (r**2 + z**2 + c**2)**3
        ur = -2 * v0 * (c**4) * r * z / (r**2 + z**2 + c**2)**3

        ux = ur
        uy = np.zeros_like(ux)

        name = 'fluid'

        m = ones_like(x)*dx*dx*ro
        h = ones_like(x)*hdx*dx
        rho = ones_like(x) * ro

        # remove particles outside the circle
        indices = []
        for i in range(len(x)):
            if not (((c - x[i])**2 + z[i]**2 - a**2) < 0) and (not (((c + x[i])**2 + z[i]**2 - a**2) < 0)) :  #and (not (((c + x[i])**2 + z[i]**2 - a**2) < 0))
                indices.append(i)

        fluid = get_particle_array(x=z, y=x, m=m, rho=rho, h=h, u=uz, v=ux,
                                name=name)
        fluid.u[indices] = 0
        fluid.w[indices] = 0

        print("Vortex ring :: %d particles"
              % (fluid.get_number_of_particles()))

        # boundary
        n = 3
        xb, yb = np.mgrid[-40-n*dx:10+n*dx:dx, -10-n*dx:10+n*dx:dx]

        xb = xb.ravel()
        yb = yb.ravel()


        m = ones_like(xb)*dx*dx*ro
        h = ones_like(xb)*hdx*dx
        rho = ones_like(xb) * ro

        boundary = get_particle_array(name='solid', x=xb, y=yb, m=m, rho=rho, h=h)
        indices = []
        for i in range(len(xb)):
            if xb[i] <= 10.0:
                if yb[i]<=10.0 and yb[i]>=-10.0:
                    indices.append(i)

        boundary.remove_particles(indices)

        self.scheme.setup_properties([fluid, boundary])
        return [fluid, boundary]


if __name__ == '__main__':
    app = VortexRing()
    app.run()
