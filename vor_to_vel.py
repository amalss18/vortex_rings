'''
Using the discretized Biot Savart law to get
the initial velocity field from the intial vorticity
field of the ring.

God Speed - Amal
'''

from mayavi import mlab
import numpy as np
from numba import prange, njit, jit, float64
import time
from fasten import vel_calc as vel_calc_2

# Geometry and SPH parameters
R = 3.0 # Radius of ring
a = 1.0  # Radius of core
nu = 1e-3
Re = 7500
Gamma = Re * nu
n = 30

length = R + 2.5
# dx = 2*length/n
dx = 0.35

vmax = 0.15*Gamma

def min_distance(point, particles):
    diff = particles - point
    difft = diff.T
    dist = np.sqrt(difft[0]**2 + difft[1]**2 + difft[2]**2)
    idx = np.where(dist == dist.min())
    return idx

@njit(float64[:, :](float64[:, :], float64[:], float64[:], float64[:], float64[:, :]), parallel=True)
def vel_calc(vorticity, x, y, z, velocity):
    for i in prange(len(x)):
        vel_sum = np.array([0.0, 0.0, 0.0])
        ref_particle = np.array([x[i], y[i], z[i]])
        # print(i, len(x))
        for j in range(len(x)):
            # print(j, i)
            if i!=j:

                diff_vec = ref_particle - np.array([x[j], y[j], z[j]])
                dist = np.linalg.norm(diff_vec)

                vel_sum = vel_sum + np.cross(vorticity[j], diff_vec) / (dist)**3

        velocity[i] = vel_sum / (4 * np.pi)

    return velocity


def vor_to_vel():
    start = time.time()
    x, y, z = np.mgrid[-length-1:length+1:dx, -length-1:length+1:dx, -length+0.75:length+3.5:dx]
    x = x.ravel()
    y = y.ravel()
    z = z.ravel()
    print(len(x), len(y), len(z))

    idx = (np.sqrt(x**2 + y**2) - R)**2 + z**2 < a**2
    # print(idx)
    # print(np.where(True==idx)[0].shape)
    # x = x[idx]
    # y = y[idx]
    # z = z[idx]

    theta = np.linspace(0, 2*np.pi, 2500)
    x_ring, y_ring = R*np.cos(theta), R*np.sin(theta)
    ring = np.vstack((x_ring, y_ring, np.zeros_like(x_ring))).T
    # ring = np.vstack((x_ring, y_ring, R*np.ones_like(x_ring))).T

    vorticity = np.zeros((len(x), 3))
    velocity = np.zeros((len(x), 3))


    # Vorticity computation
    for i in range(len(x)):
        if (np.sqrt(x[i]**2 + y[i]**2) - R)**2 + (z[i])**2 < a**2:
            particle = np.array([x[i], y[i], z[i]])
            ring_idx = min_distance(particle, ring)
            ref_point = ring[ring_idx][0]
            diff_vec = ref_point - particle
            dist = np.linalg.norm(diff_vec)
            perp_vec = np.cross(particle, diff_vec)
            if z[i] < 0:
                perp_vec = -perp_vec
            if np.linalg.norm(perp_vec) != 0:
                perp_vec = perp_vec / np.linalg.norm(perp_vec)

            vorticity[i] = Gamma * (np.e**(-(dist/a)**2)) * perp_vec / (np.pi * a**2)

    # Velocity computation

    # velocity = vel_calc(vorticity, x, y, z, velocity)
    velocity = vel_calc_2(vorticity, x, y, z, velocity)
    # print(velocity)
    velocity = np.array(velocity)

    np.savez("vel_035.npz", velocity)

    # omega_x, omega_y, omega_z = vorticity.T[0], vorticity.T[1], vorticity.T[2]
    # vorticity_vectors = mlab.quiver3d(x, y, z, omega_x, omega_y, omega_z)

    # x = x.reshape(n, n, n)
    # y = y.reshape(n, n, n)
    # z = z.reshape(n, n, n)
    # omega_x = omega_x.reshape(n, n, n)
    # omega_y = omega_y.reshape(n, n, n)
    # omega_z = omega_z.reshape(n, n, n)
    # src = mlab.pipeline.vector_field(x, y, z, omega_x, omega_y, omega_z)
    # # mlab.pipeline.vector_cut_plane(src)
    # mlab.pipeline.vectors(src)


    ##################################3
    end = time.time()
    print(end-start)

    u, v, w = velocity.T[0], velocity.T[1], velocity.T[2]
    # mlab.quiver3d(x, y, z, u, v, w)

    x = x.reshape(n, n, n)
    y = y.reshape(n, n, n)
    z = z.reshape(n, n, n)
    u = u.reshape(n, n, n)
    v = v.reshape(n, n, n)
    w = w.reshape(n, n, n)
    src = mlab.pipeline.vector_field(x, y, z, u, v, w)
    # mlab.pipeline.vector_cut_plane(src)
    mlab.pipeline.vectors(src)

if __name__=="__main__":
    vor_to_vel()
