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

# Geometry and SPH parameters
R = 3.0 # Radius of ring
a = 1.0  # Radius of core
nu = 1e-3
print("Nu: ", nu)
Re = 7500
Gamma = Re * nu
print("Gamma = ", Gamma)
n = 20

sep = 6*R

length = R + 2.5
dx = 2*length/n
print("dx: ", dx)

rho = 1000.0
hdx = 1.32
h = hdx * dx
m = rho * dx * dx * dx

vel_scale = 1.0

vmax = 0.15*vel_scale*Gamma
c0 = 10 * vmax

dt_cfl = 0.25 * h / (c0 + vmax)
dt_viscous = 0.125 * h * h / nu
dt_wcsph = 0.125 * h / c0

print(dt_cfl, dt_viscous, dt_wcsph)

dt = 0.5 * min(dt_cfl, dt_viscous, dt_wcsph)

def min_distance(point, particles):
    diff = particles - point
    difft = diff.T
    dist = np.sqrt(difft[0]**2 + difft[1]**2 + difft[2]**2)
    idx = np.where(dist == dist.min())
    return idx


class VortexRing(Application):
    def create_particles(self):
        x, y, z = np.mgrid[-length:length:dx, -length:length:dx, -length:length+sep:dx]
        x = x.ravel()
        y = y.ravel()
        z = z.ravel()

        idx = (np.sqrt(x**2 + y**2) - R)**2 + z**2 < a**2
        # x = x[idx]
        # y = y[idx]
        # z = z[idx]

        theta = np.linspace(0, 2*np.pi, 2500)
        x_ring, y_ring = R*np.cos(theta), R*np.sin(theta)
        ring = np.vstack((x_ring, y_ring, np.zeros_like(x_ring))).T
        ring_2 = np.vstack((x_ring, y_ring, sep*np.ones_like(x_ring))).T

        velocity = np.zeros((len(x), 3))

        for i in range(len(x)):
            if (np.sqrt(x[i]**2 + y[i]**2) - R)**2 + z[i]**2 < a**2:
                particle = np.array([x[i], y[i], z[i]])
                ring_idx = min_distance(particle, ring)
                ref_point = ring[ring_idx][0]
                diff_vec = ref_point - particle
                dist = np.linalg.norm(diff_vec)
                perp_vec = np.cross(particle, diff_vec)
                xn = np.array([diff_vec[0], 0])
                yn = np.array([diff_vec[1], 0])
                zn = np.array([diff_vec[2], 0])
                if np.linalg.norm(perp_vec) != 0:
                    xn, yn, zn = rotate(xn, yn, zn, perp_vec, 90) # 90 deg rotation about perp vec
                    tangent = np.array([xn[0], yn[0], zn[0]])
                    tangent = tangent/np.linalg.norm(tangent)
                else:
                    tangent = np.array([0, 0, 1])
                velocity[i] = Gamma * (1 - np.e**(-(dist/a)**2)) * tangent / (2 * np.pi * dist)

            # elif (np.sqrt(x[i]**2 + y[i]**2) - R)**2 + (z[i]-sep)**2 < a**2:
            #     particle = np.array([x[i], y[i], z[i]])
            #     ring_idx = min_distance(particle, ring_2)
            #     ref_point = ring_2[ring_idx][0]
            #     diff_vec = ref_point - particle
            #     dist = np.linalg.norm(diff_vec)
            #     perp_vec = np.cross(particle, diff_vec)

            #     xn = np.array([diff_vec[0], 0])
            #     yn = np.array([diff_vec[1], 0])
            #     zn = np.array([diff_vec[2], 0])

            #     if np.linalg.norm(perp_vec) != 0:
            #         xn, yn, zn = rotate(xn, yn, zn, perp_vec, 90) # 90 deg rotation about perp vec
            #         tangent = np.array([xn[0], yn[0], zn[0]])
            #         tangent = tangent/np.linalg.norm(tangent)
            #     else:
            #         tangent = np.array([0, 0, 1])

            #     velocity[i] = -Gamma * (1 - np.e**(-(dist/a)**2)) * tangent / (2 * np.pi * dist)


        idx_flip = z<0
        velocity[idx_flip] = -velocity[idx_flip]
        u, v, w = velocity.T[0], velocity.T[1], velocity.T[2]
        u, v, w = np.array([u, v, w]) * vel_scale


        fluid = get_particle_array(x=x, y=y, z=z, m=m, rho=rho, h=h, u=u, v=v, w=w,
                                   name="fluid")

        # fluid.tag[idx] = 1
        # single_idx = np.where(idx == True)[0][0]
        # fluid.tag[single_idx] = 2

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
        tf = 5
        solver = Solver(kernel=kernel, dim=3, integrator=integrator,
                        tf=tf, dt=dt, adaptive_timestep=False,
                        fixed_h=False)
        return solver

if __name__=="__main__":
    app = VortexRing()
    app.run()