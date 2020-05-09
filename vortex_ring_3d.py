# from pysph.base.utils import get_particle_array
# from pysph.solver.application import Application
# from pysph.sph.scheme import WCSPHScheme
# import numpy as np
# from numpy import ones_like, mgrid, sqrt
# from vField3d import vField3d
# from pysph.base.nnps import DomainManager

# class VortexRing(Application):
#     def initialize(self):
#         self.co = 10
#         self.ro = 1000.0
#         self.hdx = 1.3
#         self.dx = 0.25
#         self.alpha = 0.1
#         self.nu = 1e-6
#         self.mins = [0,0,0]
#         self.maxs = [0,0,0]

#     def create_scheme(self):
#         # s = WCSPHScheme(
#         #     ['fluid'], ['boundary'], dim=3, rho0=self.ro, c0=self.co,
#         #     h0=self.dx*self.hdx, hdx=self.hdx, gamma=7.0, alpha=0.2, beta=0.0, nu = 1e-6
#         # )
#         s = WCSPHScheme(
#             ['fluid'], [], dim=3, rho0=self.ro, c0=self.co,
#             h0=self.dx*self.hdx, hdx=self.hdx, gamma=7.0, alpha=0.1, beta=0.0, nu = 1e-6
#         )
#         dt = 1e-4
#         tf = 3
#         s.configure_solver(dt=dt, tf=tf)
#         return s

#     def create_domain(self):
#         mins = self.mins
#         maxs = self.maxs
#         return DomainManager(
#             xmin=mins[0], xmax=maxs[0], ymin=mins[1], ymax=maxs[1],
#             zmin=mins[2], zmax=maxs[2], periodic_in_x=True, periodic_in_y=True,
#             periodic_in_z=True
#         )

#     def create_particles(self):
#         print("a")
#         dx = self.dx
#         ro = self.ro
#         hdx = self.hdx

#         nl = 2
#         gamma0 = 7500 * self.nu
#         # x, y, z = np.mgrid[-hx:hx:dx, -hy:hy:dx, -hz1:hz2:dx]
#         # nop = int(2* np.pi * 0.4/dx)
#         nop = 50
#         print(nop)
#         x, y, z, ux, uy, uz = vField3d(C=np.array([[1,0,0],[0,0,0]]), r=1.5, Nvec=np.array([[1,0,0],[1,0,0]]), np=nop, f=2, dx = self.dx, gamma = 0.0075, Rcore=0.3)
#         # x_boundary, y_boundary, z_boundary = np.mgrid[x.min()-nl*dx:x.max()+nl*dx:dx, y.min()-nl*dx:y.max()+nl*dx:dx, z.min()-nl*dx:z.max()+nl*dx:dx]
#         self.mins = [x.min(), y.min(), z.min()]
#         self.maxs = [x.max(), y.max(), z.max()]
#         # x_boundary   = x_boundary.ravel()
#         # y_boundary = y_boundary.ravel()
#         # z_boundary = z_boundary.ravel()
#         theta = np.zeros_like(x)

#         for i in range(len(x)):
#             if not x[i] == 0:
#                 theta[i] = np.arctan(y[i]/x[i])
#             else:
#                 theta[i] = np.pi/2

#         r = np.sqrt(x**2+y**2)

#         ur, utheta, uz = np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)

#         v0 = 15

#         self.co = 5*v0
#         c = 1.5
#         a = 0.3
#         # a = c/5

#         uz = v0 * c**4 * (r**2 - z**2 - c**2) / (r**2 + z**2 + c**2)**3
#         ur = -2 * v0 * (c**4) * r * z / (r**2 + z**2 + c**2)**3
#         # K = gamma0
#         # # r0 = c
#         # # pi = np.pi
#         # # ur = (K * z /(2* pi * ((r-r0)**2 + z**2))) - (K*(r-r0)*z/(4 * pi * r0 * ((r-r0)**2 + z**2)))
#         # # uz = (-K *(r-r0) / (2 * pi * (r-r0)**2 + z**2)) + (K/(4*pi*r0)) * (np.log(8*r0/(np.sqrt((r-r0)**2 + z**2))) - (z**2/((r-r0)**2 + z**2))) 

#         ux = ur * np.cos(theta)
#         uy = ur * np.sin(theta)

#         # print(ux)
#         name = 'fluid'
#         vel = np.sqrt(ux**2+uy**2+uz**2)
#         c0 = 10*vel.max()    
#         print("Maximum initial velocity:",vel.max())
#         m = ones_like(x)*dx*dx*ro
#         h = ones_like(x)*hdx*dx
#         rho = ones_like(x) * ro

#         fluid = get_particle_array(x=x, y=y, z=z, m=m, rho=rho, h=h, u=ux, v=uy, w=uz,
#                                 name=name)

#         # m = ones_like(x_boundary)*dx*dx*ro
#         # h = ones_like(x_boundary)*hdx*dx
#         # rho = ones_like(x_boundary) * ro
        
#         # boundary = get_particle_array(x=x_boundary, y=y_boundary, z=z_boundary, m=m, rho=rho, h=h, name='boundary')
        
#         indices = []
#         indices_boundary = []

#         for i in range(len(x)):
#             if not (((c - r[i])**2 + z[i]**2 - a**2) < 0):
#                 indices.append(i)

#         # for i in range(len(x_boundary)):
#         #     if x_boundary[i]>=-x.min() and x_boundary[i]<=x.max():
#         #         if y_boundary[i]>=-y.min() and y_boundary[i]<=y.max():
#         #             if z_boundary[i]>=-z.min() and z_boundary[i]<=z.max():
#         #                 indices_boundary.append(i)

#         fluid.u[indices] = 0
#         fluid.v[indices] = 0
#         fluid.w[indices] = 0

#         # boundary.remove_particles(indices_boundary)

#         print("Vortex ring :: %d particles"
#               % (fluid.get_number_of_particles()))
#         # print("reynolds number:")

#         # self.scheme.setup_properties([fluid, boundary])
#         # return [fluid, boundary]

#         self.scheme.setup_properties([fluid])
#         return [fluid]

# if __name__ == '__main__':
#     app = VortexRing()
#     app.run()


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
